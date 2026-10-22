#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from scipy.io import mmread
import pyarrow.parquet as pq
try:
    import cupy as cp
except ImportError:
    cp = None
import os
try:
    from numba import cuda
except ImportError:
    class _CPUOnlyCUDA:
        @staticmethod
        def jit(function):
            return function
    cuda = _CPUOnlyCUDA()
import time
import math
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
try:
    if cp is None:
        raise RuntimeError("CuPy is not installed")
    cuda.select_device(0)
    d_array = cuda.device_array(1)
    GPU_AVAILABLE = True
    print("CUDA context initialized successfully")
except Exception as e:
    print(f"CUDA initialization failed: {e}")
    GPU_AVAILABLE = False
    print("Falling back to CPU processing")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Process Visium HD data to merge spots.')
    parser.add_argument('--mode', type=str, required=True, 
                        choices=['dense', 'proportional', 'fixed_spotSize', 'fixed_spotDistance'],
                        help='Mode for spot merging')
    parser.add_argument('--spot_size', type=int, required=True, 
                        help='Spot size in microns')
    parser.add_argument('--tissue_positions', type=str, required=True, 
                        help='Path to tissue_positions_list.csv file')
    parser.add_argument('--matrix_dir', type=str, required=True, 
                        help='Directory containing filtered_feature_bc_matrix')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Output directory')
    parser.add_argument('--gene_names', type=str,
                        help='Existing study gene list; optionally restrict sparse rows before aggregation')
    
    return parser.parse_args()

def calculate_parameters(mode, spot_size):
    """Calculate 8-um-bin aggregation geometry; fixed_spotSize input is spacing."""
    if spot_size < 8 or spot_size % 8:
        raise ValueError("The simulation parameter must be a positive multiple of the 8-um input-bin size")
    original_spot_size = 8
    if mode == "dense":
        spot_distance = spot_size
        spots_per_side = spot_size // original_spot_size
    elif mode == "proportional":
        spot_distance = round(spot_size * (100 / 55) / 8) * 8
        spots_per_side = spot_size // original_spot_size
    elif mode == "fixed_spotSize":
        spot_distance = spot_size
        spot_size = 56
        spots_per_side = spot_size // original_spot_size
    elif mode == "fixed_spotDistance":
        spot_distance = 96
        spots_per_side = spot_size // original_spot_size
    else:
        raise ValueError(f"Unknown simulation mode: {mode}")
    return {
        "spot_size": spot_size,
        "spot_distance": spot_distance,
        "spots_per_side": spots_per_side,
        "step_size": spot_distance // original_spot_size,
    }

@cuda.jit
def process_regions_kernel(position_array, regions, results, spots_per_side, step_size, max_row, max_col):
    """CUDA kernel for processing multiple regions"""
    idx = cuda.grid(1)
    if idx < len(regions):
        row, col = regions[idx]

        base_row = 1 + (row * step_size)
        base_col = 1 + (col * step_size)
        
        row_min = base_row
        row_max = min(base_row + spots_per_side - 1, max_row)
        col_min = base_col
        col_max = min(base_col + spots_per_side - 1, max_col)

        has_tissue = False
        tissue_count = 0

        min_pxl_row = 1e10  # Large initial value
        max_pxl_row = -1e10  # Small initial value
        min_pxl_col = 1e10
        max_pxl_col = -1e10

        valid_points = 0
        for r in range(row_min, row_max + 1):
            for c in range(col_min, col_max + 1):
                if not math.isnan(position_array[r, c, 0]):
                    in_tissue = position_array[r, c, 0]
                    if in_tissue > 0:
                        has_tissue = True
                        tissue_count += 1

                    pxl_row = position_array[r, c, 2]
                    pxl_col = position_array[r, c, 3]

                    if not math.isnan(pxl_row) and not math.isnan(pxl_col):
                        min_pxl_row = min(min_pxl_row, pxl_row)
                        max_pxl_row = max(max_pxl_row, pxl_row)
                        min_pxl_col = min(min_pxl_col, pxl_col)
                        max_pxl_col = max(max_pxl_col, pxl_col)
                        valid_points += 1

        if has_tissue and valid_points > 0:

            mid_pxl_row = (min_pxl_row + max_pxl_row) / 2.0
            mid_pxl_col = (min_pxl_col + max_pxl_col) / 2.0
            
            # Save results: [is_valid, row, col, mid_pxl_row, mid_pxl_col, tissue_count]
            results[idx, 0] = 1.0  # valid
            results[idx, 1] = float(row)
            results[idx, 2] = float(col)
            results[idx, 3] = mid_pxl_row
            results[idx, 4] = mid_pxl_col
            results[idx, 5] = float(tissue_count)
        else:
            results[idx, 0] = 0.0  # invalid

def process_batch_on_gpu(batch, position_array, spots_per_side, step_size, barcode_prefix, barcode_list):
    """Process a batch of regions on GPU"""
    try:

        regions_np = np.array(batch, dtype=np.int32)
        max_row, max_col = position_array.shape[0]-1, position_array.shape[1]-1

        results_np = np.zeros((len(batch), 6), dtype=np.float32)  # [is_valid, row, col, avg_pxl_row, avg_pxl_col, tissue_count]

        d_position_array = cuda.to_device(position_array)
        d_regions = cuda.to_device(regions_np)
        d_results = cuda.to_device(results_np)

        threads_per_block = 128
        min_blocks = 8
        blocks_per_grid = max(min_blocks, (len(batch) + threads_per_block - 1) // threads_per_block)

        process_regions_kernel[blocks_per_grid, threads_per_block](
            d_position_array, d_regions, d_results, 
            spots_per_side, step_size, max_row, max_col
        )

        results_np = d_results.copy_to_host()

        batch_results = []
        for i, result in enumerate(results_np):
            if result[0] > 0:  # Valid result
                row, col = batch[i]
                avg_pxl_row = result[3]
                avg_pxl_col = result[4]

                barcode_new = f"{barcode_prefix}_{int(row)}_{int(col)}"

                base_row = 1 + (row * step_size)
                base_col = 1 + (col * step_size)
                row_min = base_row
                row_max = min(base_row + spots_per_side - 1, max_row)
                col_min = base_col
                col_max = min(base_col + spots_per_side - 1, max_col)
                
                old_barcodes = []
                for r in range(row_min, row_max + 1):
                    for c in range(col_min, col_max + 1):
                        if not np.isnan(position_array[r, c, 0]) and position_array[r, c, 0] > 0:
                            barcode_idx = int(position_array[r, c, 1])
                            if barcode_idx >= 0 and barcode_idx < len(barcode_list):
                                old_barcodes.append(barcode_list[barcode_idx])

                locs_df = pd.DataFrame(
                    [[barcode_new, 1, int(row), int(col), avg_pxl_row, avg_pxl_col]],
                    columns=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
                )

                mapping_df = pd.DataFrame(
                    [[barcode_new, ','.join(old_barcodes)]],
                    columns=['new_barcode', 'old_barcodes']
                )
                
                batch_results.append((locs_df, mapping_df))
        
        return batch_results
            
    except Exception as e:
        raise RuntimeError(f"GPU processing failed for a region batch: {e}") from e

def gpu_process_regions(position_df, tech_df, params, max_row, max_col):
    """Use GPU to accelerate processing of multiple regions"""
    spots_per_side = params['spots_per_side']
    step_size = params['step_size']
    spot_size = params['spot_size']
    barcode_prefix = f"s_{spot_size:03d}um"

    max_new_row = max_row // step_size
    max_new_col = max_col // step_size

    print("Preparing data for GPU processing...")
    position_array = np.zeros((max_row+1, max_col+1, 4), dtype=np.float32)  # [in_tissue, barcode_idx, pxl_row, pxl_col]
    position_array.fill(np.nan)

    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(tech_df.index)}
    barcode_list = list(tech_df.index)

    for _, row in position_df.iterrows():
        r, c = int(row['array_row']), int(row['array_col'])
        if 0 <= r <= max_row and 0 <= c <= max_col:
            barcode = row['barcode']
            if barcode in barcode_to_idx:
                position_array[r, c, 0] = float(row['in_tissue'])
                position_array[r, c, 1] = float(barcode_to_idx.get(barcode, -1))
                position_array[r, c, 2] = float(row['pxl_row_in_fullres'])
                position_array[r, c, 3] = float(row['pxl_col_in_fullres'])

    BATCH_SIZE = 200  # Regions per batch
    
    locs_list = []
    mapping_list = []
    print(f"Processing {max_new_row * max_new_col} regions with GPU acceleration...")

    batch_args = []
    for batch_start in range(0, max_new_row * max_new_col, BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, max_new_row * max_new_col)
        batch = [(i // max_new_col, i % max_new_col) for i in range(batch_start, batch_end)]
        batch_args.append(batch)
    
    results = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(process_batch_on_gpu, 
                                  batch, 
                                  position_array.copy(), 
                                  spots_per_side, 
                                  step_size, 
                                  barcode_prefix, 
                                  barcode_list)
                  for batch in batch_args]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing batches"):
            results.extend(future.result())

    for batch_results in results:
        if batch_results:
            locs, mapping = batch_results
            locs_list.append(locs)
            mapping_list.append(mapping)

    print("Preparing for expression data merging...")
    all_mappings = pd.concat(mapping_list, axis=0) if mapping_list else pd.DataFrame(columns=['new_barcode', 'old_barcodes'])
    
    return locs_list, all_mappings

def gpu_merge_spots(tech_df, mapping_df):
    """Use GPU to accelerate merging of spot expression data"""
    try:
        if not GPU_AVAILABLE:
            raise RuntimeError("Using CPU sparse-matrix summation")

        new_barcode_groups = {}
        for _, row in mapping_df.iterrows():
            new_barcode = row['new_barcode']
            old_barcodes = set(row['old_barcodes'].split(','))

            old_barcodes = {b for b in old_barcodes if b}
            new_barcode_groups[new_barcode] = old_barcodes

        new_barcodes = list(new_barcode_groups.keys())

        result_df = pd.DataFrame(0, index=new_barcodes, columns=tech_df.columns)

        batch_size = 500  # Adjust based on GPU memory
        num_batches = math.ceil(len(new_barcodes) / batch_size)
        
        print(f"Processing expression data in {num_batches} batches...")
        with tqdm(total=len(new_barcodes)) as pbar:
            for i in range(0, len(new_barcodes), batch_size):
                batch = new_barcodes[i:i+batch_size]

                for new_barcode in batch:
                    old_barcodes_set = new_barcode_groups[new_barcode]

                    valid_old_barcodes = [b for b in old_barcodes_set if b in tech_df.index]
                    
                    if valid_old_barcodes:

                        try:

                            counts_subset = tech_df.loc[valid_old_barcodes].values

                            if counts_subset.size > 0:
                                counts_gpu = cp.array(counts_subset)
                                sum_counts = cp.sum(counts_gpu, axis=0)
                                result_df.loc[new_barcode] = cp.asnumpy(sum_counts)

                                del counts_gpu
                                cp.get_default_memory_pool().free_all_blocks()
                        except Exception as gpu_err:
                            print(f"GPU computation failed for {new_barcode}: {gpu_err}")

                            subset = tech_df.loc[valid_old_barcodes]
                            result_df.loc[new_barcode] = np.asarray(subset.sparse.to_coo().sum(axis=0)).ravel() if all(isinstance(dtype, pd.SparseDtype) for dtype in subset.dtypes) else subset.sum(axis=0).to_numpy()
                    
                    pbar.update(1)
        
        result_df.index.name = "spot"
        return result_df
        
    except Exception as e:
        print(f"GPU acceleration failed in merge_spots: {e}")
        print("Falling back to CPU processing...")

        new_barcode_groups = {}
        for _, row in mapping_df.iterrows():
            new_barcode = row['new_barcode']
            old_barcodes = set(row['old_barcodes'].split(','))
            old_barcodes = {b for b in old_barcodes if b}
            new_barcode_groups[new_barcode] = old_barcodes
        
        result_df = pd.DataFrame(0, index=list(new_barcode_groups.keys()), columns=tech_df.columns)
        
        print("Merging spots using CPU...")
        with tqdm(total=len(new_barcode_groups)) as pbar:
            for new_barcode, old_barcodes_set in new_barcode_groups.items():
                valid_old_barcodes = [b for b in old_barcodes_set if b in tech_df.index]
                if valid_old_barcodes:
                    subset = tech_df.loc[valid_old_barcodes]
                    result_df.loc[new_barcode] = np.asarray(subset.sparse.to_coo().sum(axis=0)).ravel() if all(isinstance(dtype, pd.SparseDtype) for dtype in subset.dtypes) else subset.sum(axis=0).to_numpy()
                pbar.update(1)
        
        result_df.index.name = "spot"
        return result_df

def cal_spot_loc_parallel(position, features, tech_df, params, max_dim=None):
    """Calculate new spot locations using GPU acceleration"""
    if max_dim is None:
        max_row = position['array_row'].max()
        max_col = position['array_col'].max()
    else:
        max_row, max_col = max_dim
    
    try:

        if not GPU_AVAILABLE:
            raise RuntimeError("Using the CPU aggregation path")
        print("Using GPU acceleration for spot merging...")
        locs_list, all_mappings = gpu_process_regions(position, tech_df, params, max_row, max_col)

        print("Merging expression data using GPU acceleration...")
        cnts_df = gpu_merge_spots(tech_df, all_mappings)
        
        return locs_list, [cnts_df], all_mappings.to_frame() if isinstance(all_mappings, pd.Series) else [all_mappings]
        
    except Exception as e:
        print(f"GPU acceleration failed: {e}")
        print("Falling back to CPU processing...")

        step_size = params['step_size']
        spots_per_side = params['spots_per_side']
        spot_size = params['spot_size']
        barcode_prefix = f"s_{spot_size:03d}um"

        max_new_row = max_row // step_size
        max_new_col = max_col // step_size

        locs_list = []
        mapping_list = []

        position_dict = {}
        for _, row in position.iterrows():
            r, c = int(row['array_row']), int(row['array_col'])
            if row['barcode'] in tech_df.index:
                position_dict[(r, c)] = row
        
        print(f"Processing {max_new_row * max_new_col} regions with CPU...")
        with tqdm(total=max_new_row * max_new_col) as pbar:
            for i in range(max_new_row):
                for j in range(max_new_col):

                    base_row = 1 + (i * step_size)
                    base_col = 1 + (j * step_size)
                    
                    row_min = base_row
                    row_max = min(base_row + spots_per_side - 1, max_row)
                    col_min = base_col
                    col_max = min(base_col + spots_per_side - 1, max_col)

                    region_spots = []
                    for r in range(row_min, row_max + 1):
                        for c in range(col_min, col_max + 1):
                            if (r, c) in position_dict:
                                region_spots.append(position_dict[(r, c)])
                    
                    if region_spots:
                        region_df = pd.DataFrame(region_spots)

                        if region_df['in_tissue'].any():
                            barcode_new = f"{barcode_prefix}_{i}_{j}"
                            barcode_old = region_df[region_df['in_tissue'] == 1]['barcode'].values.tolist()

                            # Match the GPU region-boundary midpoint calculation.
                            pxl_row_in_fullres = (region_df['pxl_row_in_fullres'].min() + region_df['pxl_row_in_fullres'].max()) / 2
                            pxl_col_in_fullres = (region_df['pxl_col_in_fullres'].min() + region_df['pxl_col_in_fullres'].max()) / 2

                            locs_tmp = pd.DataFrame(
                                [[barcode_new, 1, i, j, pxl_row_in_fullres, pxl_col_in_fullres]],
                                columns=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
                            )
                            locs_list.append(locs_tmp)

                            mapping_tmp = pd.DataFrame(
                                [[barcode_new, ','.join(barcode_old)]],
                                columns=['new_barcode', 'old_barcodes']
                            )
                            mapping_list.append(mapping_tmp)
                    
                    pbar.update(1)

        all_mappings = pd.concat(mapping_list, axis=0) if mapping_list else pd.DataFrame(columns=['new_barcode', 'old_barcodes'])

        print("Merging expression data...")
        cnts_df = gpu_merge_spots(tech_df, all_mappings)
        cnts_list = [cnts_df]
        
        return locs_list, cnts_list, mapping_list

def main():
    """Main function to run the script"""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    params = calculate_parameters(args.mode, args.spot_size)
    print(f"Mode: {args.mode}")
    print(f"Spot size: {params['spot_size']} um")
    print(f"Spot distance: {params['spot_distance']} um")
    print(f"Spots per side: {params['spots_per_side']}")
    print(f"Step size: {params['step_size']}")

    matrix_dir = args.matrix_dir
    print("Loading barcodes and features...")
    barcodes = pd.read_csv(f'{matrix_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz', sep='\t', header=None)
    features = pd.read_csv(f'{matrix_dir}/filtered_feature_bc_matrix/features.tsv.gz', sep='\t', header=None)
    print("Reading matrix file...")
    matrix = mmread(f'{matrix_dir}/filtered_feature_bc_matrix/matrix.mtx.gz')
    if args.gene_names:
        with open(args.gene_names, encoding='utf-8') as handle:
            requested = [line.strip() for line in handle if line.strip()]
        if not requested or len(set(requested)) != len(requested):
            raise ValueError('Study gene list must be nonempty and unique')
        selected = features[1].isin(requested)
        missing = set(requested) - set(features.loc[selected, 1])
        if missing:
            raise ValueError(f'Study genes absent from feature metadata: {sorted(missing)[:10]}')
        matrix = matrix.tocsr()[selected.to_numpy()]
        features = features.loc[selected].reset_index(drop=True)
    
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError("Matrix dimensions do not match feature and barcode files")
    print("Retaining the sparse expression matrix...")
    start_time = time.time()
    tech_matrix = matrix.tocsr()
    print(f"Matrix conversion completed in {time.time() - start_time:.2f} seconds")

    barcodes.set_index(0, inplace=True)
    features.set_index(1, inplace=True)

    print("Creating expression DataFrame...")
    tech_df = pd.DataFrame.sparse.from_spmatrix(
        tech_matrix.T,
        index=barcodes.index,
        columns=features.index
    )
    tech_df.index.name = "spot"

    print("Loading spatial positions...")
    try:

        if args.tissue_positions.endswith('.parquet'):
            parquet_file = pq.ParquetFile(args.tissue_positions)
            position = parquet_file.read().to_pandas()
        else:

            position = pd.read_csv(args.tissue_positions, sep=',')
    except Exception as e:
        raise RuntimeError(f"Error loading tissue positions file: {e}") from e

    print("Processing spot merging...")
    start_time = time.time()
    locs_list, cnts_list, mapping_list = cal_spot_loc_parallel(position, features, tech_df, params)
    print(f"Spot merging completed in {time.time() - start_time:.2f} seconds")

    if locs_list:
        locs_df = pd.concat(locs_list, axis=0)
        mapping_df = pd.concat(mapping_list, axis=0)
        cnts_df = cnts_list[0]

        output_prefix = f"{args.mode}_{params['spot_size']}um"
        
        print("Saving results...")

        locs_df.to_csv(f"{args.output_dir}/locs-raw_{output_prefix}.tsv", sep="\t", index=False)

        position_df = locs_df[['barcode', 'pxl_col_in_fullres', 'pxl_row_in_fullres']].copy()
        position_df.columns = ['spot', 'x', 'y']
        position_df.to_csv(f"{args.output_dir}/locs_{output_prefix}.tsv", sep="\t", index=False)

        cnts_df.to_csv(f"{args.output_dir}/cnts_{output_prefix}.tsv", sep="\t")

        mapping_df.to_csv(f"{args.output_dir}/barcode_mapping_{output_prefix}.tsv", sep="\t", index=False)
        
        print(f"Processing complete. Results saved to {args.output_dir}/")
        print(f"Created {len(locs_df)} merged spots")
    else:
        raise ValueError("No valid spots were found after merging")

if __name__ == "__main__":
    main()
