"""Run the external SuperST implementation and export both original outputs."""
import argparse
import gzip
import os
from pathlib import Path
import pickle
import sys
import numpy as np

def find_nonzero_bounds(image):
    rows, cols = np.where(image != 0)
    if rows.size == 0 or cols.size == 0:
        return None
    min_row = rows.min()
    max_row = rows.max() + 1
    min_col = cols.min()
    max_col = cols.max() + 1
    return image[min_row:max_row, min_col:max_col]

def store_results(predictions, gene_names, output_path):
    # Load expr
    for gene_idx in range(len(gene_names)):
        gene_data = predictions[gene_idx]
        with gzip.open(output_path+"/"+gene_names[gene_idx]+".pkl.gz", "wb", compresslevel=9) as f:
            pickle.dump(gene_data, f)

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work_dir', type=Path)
    parser.add_argument('target_s', help='Actual study target_s; passed unchanged to merge')
    parser.add_argument('--method-dir', type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    sys.path.insert(0, str(args.method_dir.resolve()))
    from app import merge
    import h5py
    prefix = str(args.work_dir.resolve())
    target_s = args.target_s
    tissue_dir = prefix
    with h5py.File(f"{prefix}/filtered_feature_bc_matrix.h5","r") as f:
        Tgenes = list(f['matrix']['features']['name'][:])
    Tgenes = [g.decode('utf-8') if isinstance(g, bytes) else g for g in Tgenes]


    demask_image_t, demask_image_t_g, tsimg1_, conv_spatial_images_ = \
        merge(tissue_dir, Tgenes, target_s= target_s, num_iter = 320, is_cut_bg = True)

    demask_image_t = [find_nonzero_bounds(arr) for arr in demask_image_t]
    demask_image_t_g = [find_nonzero_bounds(arr) for arr in demask_image_t_g]

    output_t_dir = f"{prefix}/SuperST_enhanced_t/"
    output_t_g_dir = f"{prefix}/SuperST_enhanced_t_g/"
    os.makedirs(output_t_dir, exist_ok = True)
    os.makedirs(output_t_g_dir, exist_ok = True)

    store_results(demask_image_t, Tgenes, output_t_dir)
    store_results(demask_image_t_g, Tgenes, output_t_g_dir)

if __name__ == '__main__':
    main()
