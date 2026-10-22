"""Package existing counts and coordinates without changing the simulation geometry."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil


def prepare(counts_path, locations_path, image_path, output_dir, pixel_size,
            spot_diameter, dataset, features_path=None, gene_list=None, genome='', radius_pixels=None):
    import h5py
    import numpy as np
    import pandas as pd
    from PIL import Image
    from scipy.sparse import csc_matrix
    if pixel_size <= 0 or spot_diameter <= 0:
        raise ValueError('Physical pixel size and spot diameter must be positive')
    counts = pd.read_csv(counts_path, sep='\t', index_col=0)
    locations = pd.read_csv(locations_path, sep='\t', index_col=0)
    counts.index = counts.index.astype(str)
    locations.index = locations.index.astype(str)
    if not counts.index.is_unique or not locations.index.is_unique or not counts.columns.is_unique:
        raise ValueError('Barcodes and gene names must be unique')
    if set(counts.index) != set(locations.index):
        raise ValueError('Counts and locations contain different barcodes')
    if not {'x', 'y'}.issubset(locations.columns):
        raise ValueError('Location columns must include x and y')
    locations = locations.loc[counts.index, ['x', 'y']]
    if gene_list is not None:
        genes = Path(gene_list).read_text(encoding='utf-8').splitlines()
        if not genes or len(set(genes)) != len(genes):
            raise ValueError('The supplied study gene list must be nonempty and unique')
        missing = set(genes) - set(counts.columns)
        if missing:
            raise ValueError(f'Study gene list contains absent genes: {sorted(missing)[:10]}')
        counts = counts.loc[:, genes]
    values = counts.to_numpy()
    if not counts.shape[0] or not counts.shape[1] or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts must be a nonempty finite nonnegative matrix')
    if not np.isfinite(locations.to_numpy()).all():
        raise ValueError('Locations must be finite')
    if features_path:
        features = pd.read_csv(features_path, sep='\t', dtype=str)
        if not {'id', 'name'}.issubset(features.columns) or not features['name'].is_unique:
            raise ValueError('Feature metadata requires unique name and id columns')
        features = features.set_index('name').reindex(counts.columns).rename_axis('name').reset_index()
        if features['id'].isna().any():
            raise ValueError('Feature metadata is missing expression-matrix genes')
    else:
        features = pd.DataFrame({'id': counts.columns, 'name': counts.columns})
    for field, default in [('feature_type', 'Gene Expression'), ('genome', genome)]:
        if field not in features:
            features[field] = default
        features[field] = features[field].fillna(default)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError('Use an empty destination; existing prepared inputs are not overwritten')
    spatial = output_dir / 'spatial'
    spatial.mkdir()
    counts.to_csv(output_dir / 'cnts.tsv', sep='\t', index_label='spot')
    locations.to_csv(output_dir / 'locs-raw.tsv', sep='\t', index_label='spot')
    (output_dir / 'gene-names.txt').write_text('\n'.join(counts.columns) + '\n', encoding='utf-8')
    (output_dir / 'pixel-size-raw.txt').write_text(str(pixel_size), encoding='utf-8')
    radius = spot_diameter / (2 * pixel_size) if radius_pixels is None else float(radius_pixels)
    if radius <= 0:
        raise ValueError('Radius must be positive')
    (output_dir / 'radius-raw.txt').write_text(str(radius), encoding='utf-8')
    (output_dir / 'radius.txt').write_text(str(round(radius)), encoding='utf-8')
    # Keep the source image unchanged for matrix-format methods.
    shutil.copy2(image_path, output_dir / 'he-raw.tiff')
    with Image.open(image_path) as source:
        source.convert('RGB').save(spatial / 'tissue_hires_image.png')
        source.convert('RGB').save(spatial / 'tissue_lowres_image.png')
    # Coordinates remain in the supplied full-image pixel system, so both scales are 1.
    scales = {'tissue_hires_scalef': 1.0, 'tissue_lowres_scalef': 1.0,
              'spot_diameter_fullres': spot_diameter / pixel_size,
              'microns_per_pixel': pixel_size}
    (spatial / 'scalefactors_json.json').write_text(json.dumps(scales, indent=2) + '\n', encoding='utf-8')
    positions = pd.DataFrame(index=counts.index)
    positions['in_tissue'] = 1
    positions['array_row'] = locations.y.map({v: i for i, v in enumerate(sorted(locations.y.unique()))})
    positions['array_col'] = locations.x.map({v: i for i, v in enumerate(sorted(locations.x.unique()))})
    positions['pxl_row_in_fullres'] = locations.y
    positions['pxl_col_in_fullres'] = locations.x
    positions.index.name = 'barcode'
    positions.to_csv(spatial / 'tissue_positions_list_with_header.csv')
    positions.to_csv(spatial / 'tissue_positions_list.csv', header=False)
    matrix = csc_matrix(values.T)
    string_dtype = h5py.string_dtype('utf-8')
    with h5py.File(output_dir / 'filtered_feature_bc_matrix.h5', 'w') as handle:
        handle.attrs['library_ids'] = [dataset]
        group = handle.create_group('matrix')
        for name, data in [('data', matrix.data), ('indices', matrix.indices),
                           ('indptr', matrix.indptr), ('shape', matrix.shape)]:
            group.create_dataset(name, data=data)
        group.create_dataset('barcodes', data=np.array(counts.index, dtype=object), dtype=string_dtype)
        feature_group = group.create_group('features')
        for name in ['id', 'name', 'feature_type', 'genome']:
            feature_group.create_dataset(name, data=np.array(features[name].astype(str), dtype=object), dtype=string_dtype)
        feature_group.create_dataset('_all_tag_keys', data=np.array(['genome'], dtype=object), dtype=string_dtype)
    return {'spots': len(counts), 'genes': counts.shape[1], 'output': str(output_dir)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['counts', 'locations', 'image', 'output_dir']:
        parser.add_argument('--' + name.replace('_', '-'), required=True, type=Path)
    parser.add_argument('--pixel-size', required=True, type=float, help='Micrometres per image pixel')
    parser.add_argument('--spot-diameter', required=True, type=float, help='Simulated diameter/side length in micrometres')
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--features', type=Path, help='TSV with a header and id/name columns, in any row order')
    parser.add_argument('--gene-list', type=Path, help='Actual study gene list; no automatic gene selection')
    parser.add_argument('--genome', default='', help='Actual reference assembly; optional metadata')
    parser.add_argument('--radius-pixels', type=float, help='Recorded study radius when it differs from diameter/(2*pixel-size)')
    args = parser.parse_args(argv)
    print(json.dumps(prepare(args.counts, args.locations, args.image, args.output_dir,
          args.pixel_size, args.spot_diameter, args.dataset, args.features, args.gene_list, args.genome, args.radius_pixels)))


if __name__ == '__main__':
    main()
