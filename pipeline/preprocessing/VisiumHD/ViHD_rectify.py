#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Visium HD Data Rectification Tool
This script corrects skewed parallelogram grids to rectangular grids in Visium HD data.
It crops the image before transformation to improve efficiency.
"""

import pandas as pd 
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cv2
import tifffile
import os
import sys
from scipy.stats import linregress
import argparse
import warnings
warnings.filterwarnings('ignore')

def crop_image_and_spots(spots, he_image, margin=200):
    """
    Crop image based on spot distribution with fixed pixel margins
    
    Parameters:
        spots: DataFrame with spot coordinates
        he_image: Original HE image
        margin: Fixed margin in pixels (default: 200)
    
    Returns:
        cropped_image: Cropped image
        spots_cropped: Spots with adjusted coordinates
        crop_info: Dictionary containing crop information
    """

    min_x = spots['pxl_col_in_fullres'].min()
    max_x = spots['pxl_col_in_fullres'].max()
    min_y = spots['pxl_row_in_fullres'].min()
    max_y = spots['pxl_row_in_fullres'].max()

    h, w = he_image.shape[:2]
    crop_x1 = max(0, int(min_x) - margin)
    crop_x2 = min(w, int(max_x) + margin + 1)
    crop_y1 = max(0, int(min_y) - margin)
    crop_y2 = min(h, int(max_y) + margin + 1)
    
    print(f"\nCrop Information:")
    print(f"  Original image size: {w} x {h}")
    print(f"  Spot range: X[{min_x:.0f}, {max_x:.0f}], Y[{min_y:.0f}, {max_y:.0f}]")
    print(f"  Margin: {margin} pixels")
    print(f"  Crop region: X[{crop_x1}, {crop_x2}], Y[{crop_y1}, {crop_y2}]")
    print(f"  Cropped size: {crop_x2-crop_x1} x {crop_y2-crop_y1}")
    print(f"  Area reduction: {(1 - (crop_x2-crop_x1)*(crop_y2-crop_y1)/(w*h))*100:.1f}%")

    cropped_image = he_image[crop_y1:crop_y2, crop_x1:crop_x2].copy()

    spots_cropped = spots.copy()
    spots_cropped['pxl_col_in_fullres'] = spots['pxl_col_in_fullres'] - crop_x1
    spots_cropped['pxl_row_in_fullres'] = spots['pxl_row_in_fullres'] - crop_y1

    crop_info = {
        'crop_x1': crop_x1,
        'crop_y1': crop_y1,
        'crop_x2': crop_x2,
        'crop_y2': crop_y2,
        'original_width': w,
        'original_height': h
    }
    
    return cropped_image, spots_cropped, crop_info

def rectify_visiumhd_data(spot_file, image_file=None, output_dir=None, margin=200):
    """
    Adjust HE image and spot coordinates to convert skewed parallelogram grids to 
    horizontally and vertically aligned rectangular grids.
    
    Parameters:
        spot_file: Path to Parquet file containing spot coordinates
        image_file: Path to TIFF file containing HE image (optional)
        output_dir: Directory path for saving output files
        margin: Fixed margin in pixels for cropping (default: 200)
    
    Returns:
        rectified_spots: DataFrame with corrected spot coordinates
        rectified_image: Corrected HE image (if image was provided)
    """

    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(spot_file))
    
    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, "rectified")

    print("Loading spot coordinate data...")
    try:
        spots = pd.read_parquet(spot_file)
        print(f"Successfully loaded spot data with {len(spots)} rows")
    except Exception as e:
        print(f"Error loading spot file: {e}")
        raise

    required_columns = ['array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    missing_columns = [col for col in required_columns if col not in spots.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in spot data: {missing_columns}")

    he_image_original = None
    he_image_cropped = None
    spots_working = spots.copy()
    crop_info = None
    
    if image_file is not None and os.path.exists(image_file):
        print(f"\nLoading HE image: {image_file}")
        try:
            he_image_original = tifffile.imread(image_file)
            print(f"Original image size: {he_image_original.shape}")

            print("\nPerforming image cropping...")
            he_image_cropped, spots_working, crop_info = crop_image_and_spots(
                spots, he_image_original, margin=margin
            )
            print(f"Cropped image size: {he_image_cropped.shape}")

            cropped_image_file = os.path.join(output_dir, "he_cropped.tiff")
            tifffile.imwrite(cropped_image_file, he_image_cropped)
            print(f"Cropped image saved to: {cropped_image_file}")
            
        except Exception as e:
            print(f"Error processing image: {e}")
            print("Continuing with original spot coordinates")
            spots_working = spots.copy()

    print("\nAnalyzing grid characteristics...")
    min_row, max_row = spots_working['array_row'].min(), spots_working['array_row'].max()
    min_col, max_col = spots_working['array_col'].min(), spots_working['array_col'].max()
    print(f"Grid range: array_row [{min_row}, {max_row}], array_col [{min_col}, {max_col}]")

    print("Detecting coordinate directions...")
    
    row_direction = 1
    col_direction = 1

    middle_col = int((min_col + max_col) / 2)
    col_spots = spots_working[spots_working['array_col'] == middle_col].sort_values('array_row')
    if len(col_spots) >= 2:
        min_row_spot = col_spots.iloc[0]
        max_row_spot = col_spots.iloc[-1]
        if max_row_spot['pxl_row_in_fullres'] < min_row_spot['pxl_row_in_fullres']:
            row_direction = -1
            print("Detected row number direction opposite to Y coordinate direction")

    middle_row = int((min_row + max_row) / 2)
    row_spots = spots_working[spots_working['array_row'] == middle_row].sort_values('array_col')
    if len(row_spots) >= 2:
        min_col_spot = row_spots.iloc[0]
        max_col_spot = row_spots.iloc[-1]
        if max_col_spot['pxl_col_in_fullres'] < min_col_spot['pxl_col_in_fullres']:
            col_direction = -1
            print("Detected column number direction opposite to X coordinate direction")

    print("Calculating grid skew characteristics...")

    row_slopes = []
    for row in sorted(spots_working['array_row'].unique()):
        row_spots = spots_working[spots_working['array_row'] == row]
        if len(row_spots) >= 5:
            slope, _, _, _, _ = linregress(row_spots['array_col'], row_spots['pxl_row_in_fullres'])
            row_slopes.append(slope)

    col_slopes = []
    for col in sorted(spots_working['array_col'].unique()):
        col_spots = spots_working[spots_working['array_col'] == col]
        if len(col_spots) >= 5:
            slope, _, _, _, _ = linregress(col_spots['array_row'], col_spots['pxl_col_in_fullres'])
            col_slopes.append(slope)
    
    avg_row_slope = np.median(row_slopes) if row_slopes else 0
    avg_col_slope = np.median(col_slopes) if col_slopes else 0
    print(f"Row slope (pxl_row vs array_col): {avg_row_slope:.2f} pixels/column")
    print(f"Column slope (pxl_col vs array_row): {avg_col_slope:.2f} pixels/row")

    print("Calculating grid spacing...")
    row_distances = []
    col_distances = []

    for row in spots_working['array_row'].unique():
        row_spots = spots_working[spots_working['array_row'] == row].sort_values('array_col')
        if len(row_spots) > 1:
            for i in range(len(row_spots) - 1):
                if row_spots.iloc[i+1]['array_col'] - row_spots.iloc[i]['array_col'] == 1:
                    dx = row_spots.iloc[i+1]['pxl_col_in_fullres'] - row_spots.iloc[i]['pxl_col_in_fullres']
                    dy = row_spots.iloc[i+1]['pxl_row_in_fullres'] - row_spots.iloc[i]['pxl_row_in_fullres']
                    dist = np.sqrt(dx**2 + dy**2)
                    col_distances.append(dist)

    for col in spots_working['array_col'].unique():
        col_spots = spots_working[spots_working['array_col'] == col].sort_values('array_row')
        if len(col_spots) > 1:
            for i in range(len(col_spots) - 1):
                if col_spots.iloc[i+1]['array_row'] - col_spots.iloc[i]['array_row'] == 1:
                    dx = col_spots.iloc[i+1]['pxl_col_in_fullres'] - col_spots.iloc[i]['pxl_col_in_fullres']
                    dy = col_spots.iloc[i+1]['pxl_row_in_fullres'] - col_spots.iloc[i]['pxl_row_in_fullres']
                    dist = np.sqrt(dx**2 + dy**2)
                    row_distances.append(dist)
    
    median_col_distance = np.median(col_distances)
    median_row_distance = np.median(row_distances)
    print(f"Column spacing (horizontal): {median_col_distance:.2f} pixels")
    print(f"Row spacing (vertical): {median_row_distance:.2f} pixels")

    print("Determining reference point...")
    
    def get_corner_safely(row_val, col_val, corner_name):
        corner_spots = spots_working[(spots_working['array_row'] == row_val) & 
                                     (spots_working['array_col'] == col_val)]
        if len(corner_spots) == 0:
            print(f"Warning: Exact corner point {corner_name} (row={row_val}, col={col_val}) not found")
            nearest_spots = spots_working.copy()
            nearest_spots['distance'] = (nearest_spots['array_row'] - row_val)**2 + \
                                        (nearest_spots['array_col'] - col_val)**2
            nearest_spot = nearest_spots.sort_values('distance').iloc[0]
            print(f"  Using nearest point: (row={nearest_spot['array_row']}, col={nearest_spot['array_col']})")
            return nearest_spot
        return corner_spots.iloc[0]
    
    if row_direction == 1 and col_direction == 1:
        reference_point = get_corner_safely(min_row, min_col, "top-left")
    elif row_direction == 1 and col_direction == -1:
        reference_point = get_corner_safely(min_row, max_col, "top-right")
    elif row_direction == -1 and col_direction == 1:
        reference_point = get_corner_safely(max_row, min_col, "bottom-left")
    else:
        reference_point = get_corner_safely(max_row, max_col, "bottom-right")
    
    print(f"Using reference point: row={reference_point['array_row']}, col={reference_point['array_col']}")

    print("Creating orthogonal grid coordinates...")
    
    rectified_spots = spots.copy()

    pxl_col_per_col = abs(np.median([dx / dc for dx, dc in zip(
        np.diff(row_spots.sort_values('array_col')['pxl_col_in_fullres']),
        np.diff(row_spots.sort_values('array_col')['array_col'])
    ) if dc != 0]))
    
    pxl_row_per_row = abs(np.median([dy / dr for dy, dr in zip(
        np.diff(col_spots.sort_values('array_row')['pxl_row_in_fullres']),
        np.diff(col_spots.sort_values('array_row')['array_row'])
    ) if dr != 0]))
    
    print(f"X coordinate spacing per column: {pxl_col_per_col:.2f} pixels")
    print(f"Y coordinate spacing per row: {pxl_row_per_row:.2f} pixels")
    
    ref_pxl_col = reference_point['pxl_col_in_fullres']
    ref_pxl_row = reference_point['pxl_row_in_fullres']
    ref_array_row = reference_point['array_row']
    ref_array_col = reference_point['array_col']

    for i, spot in spots_working.iterrows():
        delta_row = spot['array_row'] - ref_array_row
        delta_col = spot['array_col'] - ref_array_col
        
        rectified_spots.loc[i, 'rectified_pxl_col_in_fullres'] = ref_pxl_col + col_direction * delta_col * pxl_col_per_col
        rectified_spots.loc[i, 'rectified_pxl_row_in_fullres'] = ref_pxl_row + row_direction * delta_row * pxl_row_per_row

    print("Calculating affine transformation matrix...")
    
    ref_points = []
    n_points = min(20, len(spots_working))
    step_row = max(1, (max_row - min_row) // 5)
    step_col = max(1, (max_col - min_col) // 5)
    
    for r in range(min_row, max_row + 1, step_row):
        for c in range(min_col, max_col + 1, step_col):
            matching_spots = spots_working[(spots_working['array_row'] == r) & 
                                          (spots_working['array_col'] == c)]
            if len(matching_spots) > 0:
                orig_spot = matching_spots.iloc[0]
                rect_spot = rectified_spots[(rectified_spots['array_row'] == r) & 
                                           (rectified_spots['array_col'] == c)].iloc[0]
                
                ref_points.append({
                    'orig_x': orig_spot['pxl_col_in_fullres'],
                    'orig_y': orig_spot['pxl_row_in_fullres'],
                    'rect_x': rect_spot['rectified_pxl_col_in_fullres'],
                    'rect_y': rect_spot['rectified_pxl_row_in_fullres']
                })
                
                if len(ref_points) >= n_points:
                    break
        if len(ref_points) >= n_points:
            break
    
    if len(ref_points) < 3:
        for i in range(min(3, len(spots_working))):
            orig_spot = spots_working.iloc[i]
            rect_spot = rectified_spots.iloc[i]
            ref_points.append({
                'orig_x': orig_spot['pxl_col_in_fullres'],
                'orig_y': orig_spot['pxl_row_in_fullres'],
                'rect_x': rect_spot['rectified_pxl_col_in_fullres'],
                'rect_y': rect_spot['rectified_pxl_row_in_fullres']
            })
    
    src_points = np.array([[p['orig_x'], p['orig_y']] for p in ref_points], dtype=np.float32)
    dst_points = np.array([[p['rect_x'], p['rect_y']] for p in ref_points], dtype=np.float32)
    
    transform_matrix, inliers = cv2.estimateAffinePartial2D(
        src_points, dst_points, method=cv2.RANSAC, ransacReprojThreshold=5.0
    )
    
    print(f"Using {len(ref_points)} point pairs, {np.sum(inliers)} inliers")
    print("Affine transformation matrix:")
    print(transform_matrix)

    rectified_image = None
    if he_image_cropped is not None and transform_matrix is not None:
        print("Applying transformation to cropped image...")
        
        h, w = he_image_cropped.shape[:2]
        corners = np.array([
            [0, 0],
            [w-1, 0],
            [0, h-1],
            [w-1, h-1]
        ], dtype=np.float32).reshape(-1, 1, 2)
        
        transformed_corners = cv2.transform(corners, transform_matrix)
        
        x_min = np.floor(transformed_corners[:, 0, 0].min())
        x_max = np.ceil(transformed_corners[:, 0, 0].max())
        y_min = np.floor(transformed_corners[:, 0, 1].min())
        y_max = np.ceil(transformed_corners[:, 0, 1].max())
        
        translation_matrix = np.array([
            [1, 0, -x_min],
            [0, 1, -y_min],
            [0, 0, 1]
        ])
        
        adjusted_matrix = translation_matrix @ np.vstack([transform_matrix, [0, 0, 1]])
        adjusted_matrix = adjusted_matrix[:2, :]
        
        new_width = int(x_max - x_min + 1)
        new_height = int(y_max - y_min + 1)
        print(f"Rectified image size: {new_width} x {new_height}")
        
        rectified_image = cv2.warpAffine(
            he_image_cropped, adjusted_matrix, (new_width, new_height),
            flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0
        )

        rectified_spots['rectified_pxl_col_in_fullres'] -= x_min
        rectified_spots['rectified_pxl_row_in_fullres'] -= y_min

    print("\nEvaluating rectification effectiveness...")
    row_std_before = spots.groupby('array_row')['pxl_row_in_fullres'].std().mean()
    col_std_before = spots.groupby('array_col')['pxl_col_in_fullres'].std().mean()
    row_std_after = rectified_spots.groupby('array_row')['rectified_pxl_row_in_fullres'].std().mean()
    col_std_after = rectified_spots.groupby('array_col')['rectified_pxl_col_in_fullres'].std().mean()
    
    print(f"Original row coordinate std dev: {row_std_before:.2f} pixels")
    print(f"Original column coordinate std dev: {col_std_before:.2f} pixels")
    print(f"Rectified row coordinate std dev: {row_std_after:.2f} pixels")
    print(f"Rectified column coordinate std dev: {col_std_after:.2f} pixels")
    print(f"Row std dev improvement: {(row_std_before - row_std_after) / row_std_before * 100:.2f}%")
    print(f"Column std dev improvement: {(col_std_before - col_std_after) / col_std_before * 100:.2f}%")

    print("\nGenerating comparison visualizations...")

    plt.figure(figsize=(18, 8))
    
    plt.subplot(1, 2, 1)
    plt.scatter(spots['pxl_col_in_fullres'], spots['pxl_row_in_fullres'], 
                c=spots['array_row'], cmap='viridis', s=10, alpha=0.7)
    plt.title('Original Spot Grid (colored by row)', fontsize=14)
    plt.xlabel('pxl_col_in_fullres', fontsize=12)
    plt.ylabel('pxl_row_in_fullres', fontsize=12)
    plt.colorbar(label='array_row')
    plt.grid(alpha=0.3)
    plt.gca().invert_yaxis()
    
    plt.subplot(1, 2, 2)
    plt.scatter(rectified_spots['rectified_pxl_col_in_fullres'], 
                rectified_spots['rectified_pxl_row_in_fullres'], 
                c=rectified_spots['array_row'], cmap='viridis', s=10, alpha=0.7)
    plt.title('Rectified Spot Grid (colored by row)', fontsize=14)
    plt.xlabel('rectified_pxl_col_in_fullres', fontsize=12)
    plt.ylabel('rectified_pxl_row_in_fullres', fontsize=12)
    plt.colorbar(label='array_row')
    plt.grid(alpha=0.3)
    plt.gca().invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_grid_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(18, 8))
    
    plt.subplot(1, 2, 1)
    plt.scatter(spots['pxl_col_in_fullres'], spots['pxl_row_in_fullres'], 
                c=spots['array_col'], cmap='plasma', s=10, alpha=0.7)
    plt.title('Original Spot Grid (colored by column)', fontsize=14)
    plt.xlabel('pxl_col_in_fullres', fontsize=12)
    plt.ylabel('pxl_row_in_fullres', fontsize=12)
    plt.colorbar(label='array_col')
    plt.grid(alpha=0.3)
    plt.gca().invert_yaxis()
    
    plt.subplot(1, 2, 2)
    plt.scatter(rectified_spots['rectified_pxl_col_in_fullres'], 
                rectified_spots['rectified_pxl_row_in_fullres'], 
                c=rectified_spots['array_col'], cmap='plasma', s=10, alpha=0.7)
    plt.title('Rectified Spot Grid (colored by column)', fontsize=14)
    plt.xlabel('rectified_pxl_col_in_fullres', fontsize=12)
    plt.ylabel('rectified_pxl_row_in_fullres', fontsize=12)
    plt.colorbar(label='array_col')
    plt.grid(alpha=0.3)
    plt.gca().invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_grid_comparison_by_col.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(18, 8))
    
    selected_rows = sorted(spots['array_row'].unique())[::10]
    if len(selected_rows) < 3:
        selected_rows = sorted(spots['array_row'].unique())[:3]
    
    plt.subplot(1, 2, 1)
    for row in selected_rows:
        row_data = spots[spots['array_row'] == row].sort_values('array_col')
        plt.plot(row_data['array_col'], row_data['pxl_row_in_fullres'], 'o-', 
                 label=f'Row {row}', alpha=0.7)
    
    plt.title('Original: Y Coordinate Change per Row', fontsize=14)
    plt.xlabel('array_col', fontsize=12)
    plt.ylabel('pxl_row_in_fullres', fontsize=12)
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 2, 2)
    for row in selected_rows:
        row_data = rectified_spots[rectified_spots['array_row'] == row].sort_values('array_col')
        plt.plot(row_data['array_col'], row_data['rectified_pxl_row_in_fullres'], 'o-', 
                 label=f'Row {row}', alpha=0.7)
    
    plt.title('Rectified: Y Coordinate Change per Row', fontsize=14)
    plt.xlabel('array_col', fontsize=12)
    plt.ylabel('rectified_pxl_row_in_fullres', fontsize=12)
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_row_alignment.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(18, 8))
    
    selected_cols = sorted(spots['array_col'].unique())[::10]
    if len(selected_cols) < 3:
        selected_cols = sorted(spots['array_col'].unique())[:3]
    
    plt.subplot(1, 2, 1)
    for col in selected_cols:
        col_data = spots[spots['array_col'] == col].sort_values('array_row')
        plt.plot(col_data['array_row'], col_data['pxl_col_in_fullres'], 'o-', 
                 label=f'Col {col}', alpha=0.7)
    
    plt.title('Original: X Coordinate Change per Column', fontsize=14)
    plt.xlabel('array_row', fontsize=12)
    plt.ylabel('pxl_col_in_fullres', fontsize=12)
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.subplot(1, 2, 2)
    for col in selected_cols:
        col_data = rectified_spots[rectified_spots['array_col'] == col].sort_values('array_row')
        plt.plot(col_data['array_row'], col_data['rectified_pxl_col_in_fullres'], 'o-', 
                 label=f'Col {col}', alpha=0.7)
    
    plt.title('Rectified: X Coordinate Change per Column', fontsize=14)
    plt.xlabel('array_row', fontsize=12)
    plt.ylabel('rectified_pxl_col_in_fullres', fontsize=12)
    plt.legend()
    plt.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_col_alignment.png", dpi=300, bbox_inches='tight')
    plt.close('all')

    print("\nSaving results...")
    
    output_spots_file = os.path.join(output_dir, "tissue_positions_adj.parquet")
    rectified_spots_use = rectified_spots[["barcode", "in_tissue", "array_row", "array_col",
                                           "rectified_pxl_row_in_fullres", "rectified_pxl_col_in_fullres"]]
    rectified_spots_use.columns = ["barcode", "in_tissue", "array_row", "array_col",
                                    "pxl_row_in_fullres", "pxl_col_in_fullres"]
    rectified_spots_use.to_parquet(output_spots_file, index=False)
    print(f"Rectified spot coordinates saved to: {output_spots_file}")
    
    if rectified_image is not None:
        output_image_file = os.path.join(output_dir, "he_rectified.tiff")
        tifffile.imwrite(output_image_file, rectified_image)
        print(f"Rectified image saved to: {output_image_file}")
    
    print("\nRectification complete!")
    print(f"All output files saved to: {output_dir}")
    
    return

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description='Rectify Visium HD data from parallelogram grid to rectangular grid',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=r"""
Examples:
  # Basic usage
  python ViHD_rectify.py \\
    --spot_file /path/to/tissue_positions.parquet \\
    --image_file /path/to/tissue_image.btf \\
    --output_dir /path/to/output

  # With custom margin
  python ViHD_rectify.py \\
    --spot_file tissue_positions.parquet \\
    --image_file tissue_image.btf \\
    --output_dir data_rectified/ \\
    --margin 200

  # Without image file (only rectify spot coordinates)
  python ViHD_rectify.py \\
    --spot_file /path/to/tissue_positions.parquet \\
    --output_dir /path/to/output
        """
    )

    parser.add_argument('--spot_file', type=str, required=True,
                        help='Path to tissue_positions.parquet file')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for rectified data')

    parser.add_argument('--image_file', type=str, default=None,
                        help='Path to H&E image file (.btf or .tiff). If not provided, only spot coordinates will be rectified')
    parser.add_argument('--margin', type=int, default=200,
                        help='Fixed margin in pixels for image cropping (default: 200)')
    
    args = parser.parse_args()

    if not os.path.exists(args.spot_file):
        parser.error(f"Spot file not found: {args.spot_file}")
    
    if args.image_file is not None and not os.path.exists(args.image_file):
        parser.error(f"Image file not found: {args.image_file}")
    
    if args.margin < 0:
        parser.error("Margin must be non-negative")

    os.makedirs(args.output_dir, exist_ok=True)

    print("="*80)
    print("Visium HD Rectification Pipeline")
    print("="*80)
    print(f"Spot file: {args.spot_file}")
    print(f"Image file: {args.image_file if args.image_file else 'Not provided (spot-only mode)'}")
    print(f"Output directory: {args.output_dir}")
    print(f"Crop margin: {args.margin} pixels")
    print("="*80)
    
    try:

        rectify_visiumhd_data(
            spot_file=args.spot_file,
            image_file=args.image_file,
            output_dir=args.output_dir,
            margin=args.margin
        )
        
        print("\n" + "="*80)
        print("Rectification completed successfully!")
        print(f"Output files saved to: {args.output_dir}")
        print("="*80)
        print("\nGenerated files:")
        print(f"  - tissue_positions_adj.parquet: Rectified spot coordinates")
        if args.image_file:
            print(f"  - he_cropped.tiff: Cropped H&E image")
            print(f"  - he_rectified.tiff: Rectified H&E image")
        print(f"  - rectified_grid_comparison.png: Grid comparison (by row)")
        print(f"  - rectified_grid_comparison_by_col.png: Grid comparison (by column)")
        print(f"  - rectified_row_alignment.png: Row alignment visualization")
        print(f"  - rectified_col_alignment.png: Column alignment visualization")
        
    except Exception as e:
        print(f"\nError during rectification: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

