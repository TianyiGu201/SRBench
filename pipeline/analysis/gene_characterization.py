"""Compute platform-specific gene annotations."""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
from scipy.sparse import csr_matrix

def read_cosmx_to_adata(expr_mat_file, metadata_file, tx_file=None, polygons_file=None, prefix='CHYCH2'):
    """Build an AnnData object from CosMx cell counts and metadata.

    Counts require fov and cell_ID columns; metadata requires cell_ID and the
    global cell identifier. Optional transcript and polygon files are retained
    in the returned object. prefix identifies the dataset.
    """
    print('=' * 70)
    print('Reading CosMx data')
    print('=' * 70)
    print('\nStep 1: Reading the expression matrix...')
    expr_df = pd.read_csv(expr_mat_file, compression='gzip')
    print(f'  - Original shape: {expr_df.shape}')
    print(f'  - Example column names: {list(expr_df.columns[:10])}')
    if 'fov' not in expr_df.columns or 'cell_ID' not in expr_df.columns:
        raise ValueError("The expression matrix must contain 'fov' and 'cell_ID' columns")
    print(f"  - Number of FOVs: {expr_df['fov'].nunique()}")
    print(f"  - FOV range: {expr_df['fov'].min()} - {expr_df['fov'].max()}")
    fov_cellid_map = expr_df[['fov', 'cell_ID']].copy()
    expr_df = expr_df.set_index(['fov', 'cell_ID'])
    print(f'  - Shape after indexing: {expr_df.shape}')
    print('\nStep 2: Reading metadata...')
    metadata = pd.read_csv(metadata_file, compression='gzip')
    print(f'  - Original shape: {metadata.shape}')
    print(f'  - Number of columns: {len(metadata.columns)}')
    if 'cell_ID' not in metadata.columns:
        raise ValueError("Metadata must contain the 'cell_ID' column")
    if 'cell' not in metadata.columns:
        raise ValueError("Metadata must contain the global identifier column 'cell'")
    print(f"  - First cell identifiers: {metadata['cell'].head().tolist()}")
    print(f"  - cell_ID range: {metadata['cell_ID'].min()} - {metadata['cell_ID'].max()}")
    print('\nStep 3: Processing FOV identifiers...')
    if 'fov' not in metadata.columns:
        print('  - No fov column in metadata; parsing cell identifiers...')
        cell_parts = metadata['cell'].str.split('_', expand=True)
        numeric_cols = []
        for i in range(1, cell_parts.shape[1]):
            if cell_parts[i].str.isnumeric().all():
                numeric_cols.append(i)
        if len(numeric_cols) >= 2:
            metadata['fov_from_cell'] = cell_parts[numeric_cols[-2]].astype(int)
            print(f"  - FOV identifiers parsed from cell identifiers: {metadata['fov_from_cell'].unique()[:10]}")
        else:
            print('  - Cannot parse FOV from cell identifiers; using expression-matrix FOV information')
            metadata['fov_from_cell'] = None
    else:
        print('  - Metadata already has a fov column')
        metadata['fov_from_cell'] = metadata['fov']
    print('\nStep 4: Merging FOV identifiers...')
    if metadata['cell_ID'].duplicated().any():
        print('  Warning: Duplicate cell_ID values in metadata require expression-matrix FOV identifiers')
        if metadata['fov_from_cell'].notna().all():
            metadata['fov'] = metadata['fov_from_cell']
        else:
            print('  - Warning: FOV identifiers are ambiguous; using the first cell_ID match')
            metadata = metadata.merge(fov_cellid_map.drop_duplicates(subset='cell_ID', keep='first'), on='cell_ID', how='left', suffixes=('', '_expr'))
            if 'fov_expr' in metadata.columns:
                metadata['fov'] = metadata['fov_expr']
                metadata.drop(columns=['fov_expr'], inplace=True)
    else:
        metadata = metadata.merge(fov_cellid_map, on='cell_ID', how='left', suffixes=('', '_expr'))
        if 'fov_expr' in metadata.columns:
            metadata['fov'] = metadata['fov_expr']
            metadata.drop(columns=['fov_expr'], inplace=True)
        elif 'fov' not in metadata.columns and 'fov_from_cell' in metadata.columns:
            metadata['fov'] = metadata['fov_from_cell']
    if 'fov_from_cell' in metadata.columns:
        metadata.drop(columns=['fov_from_cell'], inplace=True)
    if 'fov' not in metadata.columns:
        raise ValueError('Cannot add the fov column to metadata')
    print(f"  - Metadata now has a fov column: {metadata['fov'].unique()[:10]}")
    print(f"  - FOV summary: {metadata['fov'].value_counts().head()}")
    print('\nStep 5: Aligning data...')
    metadata_indexed = metadata.set_index(['fov', 'cell_ID'])
    common_index = expr_df.index.intersection(metadata_indexed.index)
    print(f'  - Expression-matrix cell count: {len(expr_df)}')
    print(f'  - Metadata cell count: {len(metadata_indexed)}')
    print(f'  - Shared cell count: {len(common_index)}')
    if len(common_index) == 0:
        raise ValueError('The expression matrix and metadata have no shared cells')
    expr_df_aligned = expr_df.loc[common_index]
    metadata_aligned = metadata_indexed.loc[common_index]
    print('\nStep 6: Creating the AnnData object...')
    cell_id_map = metadata.set_index(['fov', 'cell_ID'])['cell']
    unique_cell_ids = cell_id_map.loc[common_index].values
    adata = ad.AnnData(X=csr_matrix(expr_df_aligned.values), obs=metadata_aligned.reset_index(), var=pd.DataFrame(index=expr_df_aligned.columns))
    adata.obs.index = unique_cell_ids
    adata.obs_names_make_unique()
    print('\nStep 7: Adding spatial coordinates...')
    if 'CenterX_local_px' in adata.obs.columns and 'CenterY_local_px' in adata.obs.columns:
        adata.obsm['spatial'] = adata.obs[['CenterX_local_px', 'CenterY_local_px']].values
        print(f'  Added local coordinates')
        print(f"    X range: [{adata.obsm['spatial'][:, 0].min():.2f}, {adata.obsm['spatial'][:, 0].max():.2f}]")
        print(f"    Y range: [{adata.obsm['spatial'][:, 1].min():.2f}, {adata.obsm['spatial'][:, 1].max():.2f}]")
    if 'CenterX_global_px' in adata.obs.columns and 'CenterY_global_px' in adata.obs.columns:
        adata.obsm['spatial_global'] = adata.obs[['CenterX_global_px', 'CenterY_global_px']].values
        print(f'  Added global coordinates')
    if tx_file:
        try:
            print('\nStep 8: Reading transcript information...')
            tx_df = pd.read_csv(tx_file, compression='gzip')
            adata.uns['transcripts'] = tx_df
            print(f'  Transcript count: {len(tx_df):,}')
        except Exception as e:
            print(f'  Warning: Failed to read the transcript file: {e}')
    if polygons_file:
        try:
            print('\nStep 9: Reading cell polygons...')
            polygons_df = pd.read_csv(polygons_file, compression='gzip')
            adata.uns['polygons'] = polygons_df
            print(f'  Polygon data: {polygons_df.shape}')
        except Exception as e:
            print(f'  Warning: Failed to read the polygon file: {e}')
    print('\nStep 10: Calculating summary statistics...')
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    adata.var['n_cells'] = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata.var['total_counts'] = np.array(adata.X.sum(axis=0)).flatten()
    print('\n' + '=' * 70)
    print('Data summary')
    print('=' * 70)
    print(f'Cell count: {adata.n_obs:,}')
    print(f'Gene count: {adata.n_vars:,}')
    print(f"Number of FOVs: {adata.obs['fov'].nunique()}")
    print(f"Total UMI count: {adata.obs['n_counts'].sum():,.0f}")
    print(f"Mean UMI count per cell: {adata.obs['n_counts'].mean():.2f} ± {adata.obs['n_counts'].std():.2f}")
    print(f"Mean gene count per cell: {adata.obs['n_genes'].mean():.2f} ± {adata.obs['n_genes'].std():.2f}")
    print('\nCells per FOV:')
    fov_counts = adata.obs['fov'].value_counts().sort_index()
    for fov, count in fov_counts.head(10).items():
        print(f'  FOV {fov}: {count:,} cells')
    if len(fov_counts) > 10:
        print(f'  ... (total {len(fov_counts)} FOVs)')
    return adata

def characterize_xenium(xenium_dir, output_csv):
    XENIUM_DIR = str(xenium_dir)

    adata = sc.read_10x_h5(
        os.path.join(XENIUM_DIR, "cell_feature_matrix.h5")
    )


    cells_path = os.path.join(XENIUM_DIR, "cells.csv.gz")
    cells = pd.read_csv(cells_path)

    cells["cell_id"] = cells["cell_id"].astype(str)

    adata.obs["x_coord"] = cells.set_index("cell_id").loc[
        adata.obs_names, "x_centroid"
    ].values

    adata.obs["y_coord"] = cells.set_index("cell_id").loc[
        adata.obs_names, "y_centroid"
    ].values

    adata.obsm["spatial"] = adata.obs[["x_coord", "y_coord"]].values

    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)


    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        min_mean=0.01,
        max_mean=3,
        min_disp=0.5
    )


    hvg_df = (
        adata.var[["highly_variable"]]
        .reset_index()
        .rename(columns={"index": "gene"})
    )

    hvg_df["hvg_status"] = np.where(
        hvg_df["highly_variable"],
        "HVG",
        "non-HVG"
    )

    hvg_df = hvg_df.rename(columns={hvg_df.columns[0]: "gene"})
    hvg_df = hvg_df[["gene", "hvg_status"]]


    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        delaunay=True
    )


    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=adata.var_names
    )

    moran_df = (
        adata.uns["moranI"]
        .reset_index()
        .rename(columns={"index": "gene", "I": "morans_I"})
    )

    moran_df = moran_df.rename(columns={moran_df.columns[0]: "gene"})
    moran_df = moran_df[["gene", "morans_I"]]
    gene_stats = hvg_df.merge(moran_df, on='gene', how='inner')
    gene_stats.to_csv(output_csv, index=False)
    return gene_stats

def characterize_cosmx(expr_mat_file, metadata_file, fov_ids, output_csv, prefix="sample"):
    adata = read_cosmx_to_adata(expr_mat_file, metadata_file, prefix=prefix)
    adata = adata[adata.obs['fov'].isin(list(fov_ids))].copy()
    control = adata.var_names.str.startswith('SystemControl') | adata.var_names.str.startswith('Negative')
    adata = adata[:, ~control].copy()
    adata.obsm["spatial"] = adata.obs[["CenterX_global_px", "CenterY_global_px"]].values

    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)


    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        min_mean=0.01,
        max_mean=3,
        min_disp=0.5
    )


    hvg_df = (
        adata.var[["highly_variable"]]
        .reset_index()
        .rename(columns={"index": "gene"})
    )

    hvg_df["hvg_status"] = np.where(
        hvg_df["highly_variable"],
        "HVG",
        "non-HVG"
    )

    hvg_df = hvg_df.rename(columns={hvg_df.columns[0]: "gene"})
    hvg_df = hvg_df[["gene", "hvg_status"]]


    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        delaunay=True
    )


    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=adata.var_names
    )

    moran_df = (
        adata.uns["moranI"]
        .reset_index()
        .rename(columns={"index": "gene", "I": "morans_I"})
    )

    moran_df = moran_df.rename(columns={moran_df.columns[0]: "gene"})
    moran_df = moran_df[["gene", "morans_I"]]
    gene_stats = hvg_df.merge(moran_df, on='gene', how='inner')
    gene_stats.to_csv(output_csv, index=False)
    return gene_stats

def characterize_visiumhd(tenx_h5, gene_names_file, tissue_positions_path, output_csv):
    GENE_NAMES_FILE = str(gene_names_file)
    OUT_CSV = str(output_csv)
    import scanpy as sc
    adata = sc.read_10x_h5(tenx_h5)
    adata.var_names_make_unique()
    with open(GENE_NAMES_FILE, 'r') as f:
        gene_list = [line.strip() for line in f]
    adata = adata[:, adata.var_names.isin(gene_list)]
    positions_path = tissue_positions_path
    positions = pd.read_parquet(positions_path)
    positions['barcode'] = positions['barcode'].astype(str)
    adata.obs['x_coord'] = positions.set_index('barcode').loc[adata.obs_names, 'pxl_col_in_fullres'].values
    adata.obs['y_coord'] = positions.set_index('barcode').loc[adata.obs_names, 'pxl_row_in_fullres'].values
    adata.obsm['spatial'] = adata.obs[['x_coord', 'y_coord']].values
    sc.pp.normalize_total(adata, target_sum=10000.0)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', min_mean=0.01, max_mean=3, min_disp=0.5)
    hvg_df = adata.var[['highly_variable']].reset_index().rename(columns={'index': 'gene'})
    hvg_df['hvg_status'] = np.where(hvg_df['highly_variable'], 'HVG', 'non-HVG')
    hvg_df = hvg_df.rename(columns={hvg_df.columns[0]: 'gene'})
    hvg_df = hvg_df[['gene', 'hvg_status']]
    sq.gr.spatial_neighbors(adata, coord_type='grid', delaunay=False)
    sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names)
    moran_df = adata.uns['moranI'].reset_index().rename(columns={'index': 'gene', 'I': 'morans_I'})
    moran_df = moran_df.rename(columns={moran_df.columns[0]: 'gene'})
    moran_df = moran_df[['gene', 'morans_I']]
    gene_stats = hvg_df.merge(moran_df, on='gene', how='inner')
    gene_stats
    print(gene_stats['hvg_status'].value_counts())
    print(gene_stats.groupby('hvg_status')['morans_I'].describe())
    gene_stats.to_csv(OUT_CSV, index=False)
    return gene_stats

