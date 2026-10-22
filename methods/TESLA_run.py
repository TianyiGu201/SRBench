"""Run TESLA with the supplied study settings and export gene-wise maps.

Third-party TESLA retains its MIT license; see THIRD_PARTY_NOTICES.md.
"""
import argparse
import gzip
import os
from pathlib import Path
import pickle
from multiprocessing import Pool

def _save_one_gene(args):
    import numpy as np
    (gene_raw, gene_idx, Xcsc, Xdense,
     i_coords, j_coords, height, width, out_dir, compresslevel) = args

    if Xcsc is not None:
        expr = Xcsc[:, gene_idx].toarray().ravel()
    else:
        expr = Xdense[:, gene_idx].ravel()
    gene_pd = np.zeros((height, width), dtype=np.float32)
    gene_pd[i_coords, j_coords] = expr.astype(np.float32, copy=False)
    fname = f"{gene_raw}.pkl.gz"
    out_path = os.path.join(out_dir, fname)
    with gzip.open(out_path, "wb", compresslevel=compresslevel) as f:
        pickle.dump(gene_pd, f)

def read_lines(filename):
    with open(filename, 'r') as file:
        lines = [line.rstrip() for line in file]
    return lines

def load_genes_by_name(compressed_dir, gene_name):
    compressed_dir = Path(compressed_dir)
    gene_file = compressed_dir / f"{gene_name}.pkl.gz"
    if gene_file.exists():
        with gzip.open(gene_file, 'rb') as f:
            gene_data = pickle.load(f)
    
        return gene_data

def h5ad2pkl(enhanced_res, data_dir, patch_size, gene_use_raw, n_jobs=8, compresslevel = 9):
    import numpy as np
    from scipy import sparse
    # output dir: data_dir/tesla_enhanced/
    out_dir = os.path.join(data_dir, "tesla_enhanced")
    os.makedirs(out_dir, exist_ok=True)
    # coords(/patch_size - 1)
    coords = enhanced_res.obs
    i_coords = coords['x'].astype(int) // int(patch_size) - 1
    j_coords = coords['y'].astype(int) // int(patch_size) - 1
    # img size
    height = int(i_coords.max()) + 1
    width  = int(j_coords.max()) + 1
    var_index_upper = {str(g).upper(): i for i, g in enumerate(enhanced_res.var_names)}
    X = enhanced_res.X
    Xcsc = X.tocsc() if sparse.issparse(X) else None
    Xdense = None if sparse.issparse(X) else np.asarray(X)
    tasks, missing = [], []
    for gene_raw in gene_use_raw:
        idx = var_index_upper.get(gene_raw.upper(), None)
        if idx is None:
            missing.append(gene_raw)
            continue
        tasks.append((gene_raw, idx, Xcsc, Xdense, i_coords, j_coords, height, width, out_dir, compresslevel))
    with Pool(processes=n_jobs) as pool:
        pool.map(_save_one_gene, tasks)

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_dir')
    parser.add_argument('he_path')
    parser.add_argument('work_dir')
    parser.add_argument('patch_size')
    parser.add_argument('gene_names')
    args = parser.parse_args(argv)
    import cv2
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scanpy import read_10x_h5
    from scipy.sparse import issparse
    import TESLA as tesla
    data_dir, he_dir, work_dir = args.data_dir, args.he_path, args.work_dir
    targetRES, gene_name_path = args.patch_size, args.gene_names
    os.makedirs(work_dir, exist_ok=True)
    img=cv2.imread(he_dir)
    resize_factor=1000/np.min(img.shape[0:2])

    resize_width=int(img.shape[1]*resize_factor)
    resize_height=int(img.shape[0]*resize_factor)
    cnt=tesla.cv2_detect_contour(img, apertureSize=5,L2gradient = True)
    binary=np.zeros((img.shape[0:2]), dtype=np.uint8)
    cv2.drawContours(binary, [cnt], -1, (1), thickness=-1)
    # Enlarged filter
    cnt_enlarged = tesla.scale_contour(cnt, 1.05)
    binary_enlarged = np.zeros(img.shape[0:2])
    cv2.drawContours(binary_enlarged, [cnt_enlarged], -1, (1), thickness=-1)
    img_new = img.copy()
    cv2.drawContours(img_new, [cnt], -1, (255), thickness=50)
    img_new=cv2.resize(img_new, ((resize_width, resize_height)))

    if os.path.exists(f"{data_dir}/feature_bc_matrix.h5"):
      adata = read_10x_h5(f"{data_dir}/feature_bc_matrix.h5")
    else:
      adata = read_10x_h5(f"{data_dir}/filtered_feature_bc_matrix.h5")

    spatial=pd.read_csv(f"{data_dir}/spatial/tissue_positions_list_with_header.csv",sep=",",header=0,na_filter=False,index_col=0)
    adata.obs["in_tissue"]=spatial["in_tissue"]
    adata.obs["array_y"]=spatial["array_row"]
    adata.obs["array_x"]=spatial["array_col"]
    adata.obs["pixel_y"]=spatial["pxl_row_in_fullres"]
    adata.obs["pixel_x"]=spatial["pxl_col_in_fullres"]
    # Select captured samples
    adata=adata[adata.obs["in_tissue"]==1]
    adata.var_names=[i.upper() for i in list(adata.var_names)]
    adata.var["genename"]=adata.var.index.astype("str")
    adata.write_h5ad(f"{work_dir}/tesla_data.h5ad")


    # Read in gene expression and spatial location
    counts=sc.read(f"{work_dir}/tesla_data.h5ad")
    counts.var.index=[i.upper() for i in counts.var.index]
    counts.var_names_make_unique()
    counts.raw=counts
    sc.pp.log1p(counts) # impute on log scale
    if issparse(counts.X):counts.X=counts.X.toarray()
    cv2.imwrite(f"{work_dir}/cnt.jpg", img_new)
    enhanced_exp_adata=tesla.imputation(
        img=img, 
        raw=counts, 
        cnt=cnt, 
        genes=counts.var.index.tolist(), 
        shape="None", 
        res=int(targetRES), 
        s=1, k=2, num_nbs=10
    )

    gene_use_raw = read_lines(gene_name_path)
    h5ad2pkl(enhanced_exp_adata, work_dir, targetRES, gene_use_raw, n_jobs=8, compresslevel = 6)

if __name__ == '__main__':
    main()
