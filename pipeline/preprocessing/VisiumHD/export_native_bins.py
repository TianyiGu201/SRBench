"""Export native 8-micrometre bins and matched locations for HEARTSVG input."""
from __future__ import annotations
import argparse
from pathlib import Path


def export_native_bins(matrix_dir, positions, output_dir, chunk_size=256):
    import h5py
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    from scipy.sparse import csc_matrix

    matrix_dir, positions, output_dir = map(Path, (matrix_dir, positions, output_dir))
    if chunk_size < 1:
        raise ValueError("Chunk size must be positive")
    h5_path = matrix_dir if matrix_dir.is_file() else matrix_dir / 'filtered_feature_bc_matrix.h5'
    if h5_path.is_file():
        with h5py.File(h5_path, 'r') as handle:
            group = handle['matrix']
            matrix = csc_matrix((group['data'][:], group['indices'][:], group['indptr'][:]),
                                shape=tuple(group['shape'][:]))
            barcodes = group['barcodes'].asstr()[:]
            genes = group['features/name'].asstr()[:]
    else:
        folder = matrix_dir / 'filtered_feature_bc_matrix'
        if not folder.is_dir():
            folder = matrix_dir

        def locate(name):
            for suffix in ('.gz', ''):
                path = folder / (name + suffix)
                if path.is_file():
                    return path
            raise FileNotFoundError(folder / name)

        path = locate('matrix.mtx')
        if path.suffix == '.gz':
            import gzip
            with gzip.open(path, 'rb') as handle:
                matrix = mmread(handle).tocsc()
        else:
            with path.open('rb') as handle:
                matrix = mmread(handle).tocsc()
        barcodes = pd.read_csv(locate('barcodes.tsv'), sep='\t', header=None, dtype=str)[0].to_numpy()
        genes = pd.read_csv(locate('features.tsv'), sep='\t', header=None, dtype=str)[1].to_numpy()
    if not len(genes) or not len(barcodes):
        raise ValueError("Native bin matrix must be nonempty")
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError("Matrix dimensions do not match feature and barcode metadata")
    if not pd.Index(barcodes).is_unique or not pd.Index(genes).is_unique:
        raise ValueError("Native bin barcodes and gene names must be unique")
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("Native counts must be finite and nonnegative")
    locs = pd.read_parquet(positions) if positions.suffix == '.parquet' else pd.read_csv(positions)
    required = ['barcode', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
    if not set(required).issubset(locs.columns):
        raise ValueError(f"Positions require columns: {required}")
    locs['barcode'] = locs['barcode'].astype(str)
    if not locs['barcode'].is_unique:
        raise ValueError("Position barcodes must be unique")
    locs = locs.set_index('barcode')
    missing = set(barcodes) - set(locs.index)
    if missing:
        raise ValueError(f"Positions missing for {len(missing)} expression barcodes")
    locs = locs.loc[barcodes, required[1:]].rename(
        columns={'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'})
    if not np.isfinite(locs.to_numpy(dtype=float)).all():
        raise ValueError("Native bin positions must be finite")
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise FileExistsError("Use an empty destination directory")
    locs.to_csv(output_dir / 'locs-raw.tsv', sep='\t', index_label='spot')
    for start in range(0, len(barcodes), chunk_size):
        stop = min(start + chunk_size, len(barcodes))
        block = pd.DataFrame(matrix[:, start:stop].toarray().T,
                             index=barcodes[start:stop], columns=genes)
        block.to_csv(output_dir / 'cnts.tsv', sep='\t', index_label='spot',
                     mode='w' if start == 0 else 'a', header=start == 0)
    return {'bins': len(barcodes), 'genes': len(genes)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix-dir', type=Path, required=True,
                        help='Native square_008um directory, matrix folder or 10x H5 file')
    parser.add_argument('--positions', type=Path, required=True,
                        help='Matching native-bin positions parquet or headered CSV')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--chunk-size', type=int, default=256,
                        help='Number of expression rows materialized per TSV write')
    args = parser.parse_args(argv)
    print(export_native_bins(args.matrix_dir, args.positions, args.output_dir, args.chunk_size))


if __name__ == '__main__':
    main()
