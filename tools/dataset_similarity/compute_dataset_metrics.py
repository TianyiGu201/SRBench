#!/usr/bin/env python3
"""Batch-compute dataset feature tables for spatial transcriptomics similarity search."""

from __future__ import annotations

import argparse
import json
import logging
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import sparse

OUTPUT_COLUMNS = [
    "dataset_name",
    "he_image_feature_json",
    "he_image_feature_dim",
    "tissue_type",
    "spatial_resolution_um_per_pixel",
    "data_sparsity",
    "spatial_autocorrelation",
    "mean_expression_per_100um2",
    "cv_expression_density",
    "svg_number",
]

IMAGE_PATTERNS = ["*.tif", "*.tiff", "*.png", "*.jpg", "*.jpeg"]


def setup_logger(log_path: Optional[str], verbose: bool) -> logging.Logger:
    logger = logging.getLogger("compute_dataset_metrics_batch")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG if verbose else logging.INFO)
    sh.setFormatter(formatter)
    logger.addHandler(sh)

    if log_path:
        log_file = Path(log_path)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    return logger


def parse_mapping_arg(value: Optional[str], arg_name: str) -> Dict[str, object]:
    if not value:
        return {}
    value = value.strip()
    path = Path(value)
    text = path.read_text(encoding="utf-8") if path.exists() else value
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(f"{arg_name} must be a JSON string or a path to a JSON file") from e
    if not isinstance(obj, dict):
        raise ValueError(f"{arg_name} must decode to a JSON object / dictionary")
    return obj


def find_single_existing_file(dataset_dir: Path, patterns: Iterable[str], required: bool = False) -> Optional[Path]:
    pats = list(patterns)
    matches: List[Path] = []
    seen = set()
    for pattern in pats:
        for p in sorted(dataset_dir.glob(pattern)):
            rp = str(p.resolve())
            if rp not in seen:
                matches.append(p)
                seen.add(rp)
    if not matches:
        if required:
            raise FileNotFoundError(f"No file found in {dataset_dir} matching any of: {pats}")
        return None
    if len(matches) > 1:
        raise FileExistsError(f"Multiple files found in {dataset_dir} matching {pats}: {[p.name for p in matches]}")
    return matches[0]


def discover_dataset_dirs(root_dir: Path) -> List[Path]:
    dirs = [p for p in sorted(root_dir.iterdir()) if p.is_dir()]
    if not dirs:
        raise ValueError(f"No dataset directories found under: {root_dir}")
    return dirs


def read_input_as_adata(h5ad_path: Optional[Path], cnts_path: Optional[Path], locs_path: Optional[Path]) -> ad.AnnData:
    if h5ad_path is not None:
        adata = ad.read_h5ad(h5ad_path)
    else:
        if cnts_path is None or locs_path is None:
            raise ValueError("Need either one h5ad file, or both cnts.tsv and locs-raw.tsv")
        cnts = pd.read_csv(cnts_path, sep="\t")
        if cnts.columns[0] != "spot":
            raise ValueError("The first column of cnts.tsv must be named 'spot'")
        cnts = cnts.set_index("spot")

        locs = pd.read_csv(locs_path, sep="\t")
        if list(locs.columns[:3]) != ["spot", "x", "y"]:
            raise ValueError("locs-raw.tsv must have columns: spot, x, y")
        locs = locs[["spot", "x", "y"]].set_index("spot")

        common = cnts.index.intersection(locs.index)
        if len(common) == 0:
            raise ValueError("No overlapping spot IDs between cnts.tsv and locs-raw.tsv")

        cnts = cnts.loc[common]
        locs = locs.loc[common]
        adata = ad.AnnData(X=sparse.csr_matrix(cnts.to_numpy(dtype=np.float32)))
        adata.obs_names = cnts.index.astype(str)
        adata.var_names = cnts.columns.astype(str)
        adata.obs["x"] = locs["x"].to_numpy(dtype=float)
        adata.obs["y"] = locs["y"].to_numpy(dtype=float)
        adata.obsm["spatial"] = locs[["x", "y"]].to_numpy(dtype=float)

    if "spatial" not in adata.obsm:
        if {"x", "y"}.issubset(adata.obs.columns):
            adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy(dtype=float)
        else:
            raise ValueError("Spatial coordinates not found. Expect adata.obsm['spatial'] or adata.obs[['x', 'y']]")

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def to_dense(x) -> np.ndarray:
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def compute_data_sparsity(adata: ad.AnnData) -> float:
    x = adata.X
    if sparse.issparse(x):
        zeros = x.shape[0] * x.shape[1] - x.count_nonzero()
    else:
        zeros = np.count_nonzero(np.asarray(x) == 0)
    total = adata.n_obs * adata.n_vars
    return float(zeros / total)


def estimate_spot_area_um2(adata: ad.AnnData, um_per_pixel: float) -> float:
    if um_per_pixel <= 0:
        raise ValueError("Spatial resolution must be positive (um/pixel)")
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    if coords.shape[0] < 2:
        raise ValueError("At least two spots are required to estimate spot area from coordinates")
    from sklearn.neighbors import NearestNeighbors

    nbrs = NearestNeighbors(n_neighbors=min(2, coords.shape[0]), metric="euclidean")
    nbrs.fit(coords)
    distances, _ = nbrs.kneighbors(coords)
    if distances.shape[1] < 2:
        raise ValueError("Could not estimate nearest-neighbor spot spacing")
    nn_dist_px = distances[:, 1]
    nn_dist_px = nn_dist_px[np.isfinite(nn_dist_px) & (nn_dist_px > 0)]
    if nn_dist_px.size == 0:
        raise ValueError("Could not estimate positive nearest-neighbor spacing from spot coordinates")
    median_spacing_px = float(np.median(nn_dist_px))
    diameter_um = median_spacing_px * um_per_pixel
    return float(math.pi * (diameter_um / 2.0) ** 2)


def compute_expression_density_stats(adata: ad.AnnData, area_spot_um2: float) -> Tuple[float, float]:
    if area_spot_um2 <= 0:
        raise ValueError("Spot area must be positive")
    x = adata.X
    spot_totals = np.asarray(x.sum(axis=1)).reshape(-1).astype(float)
    density = 100.0 * spot_totals / area_spot_um2
    mean_val = float(density.mean()) if density.size else np.nan
    std_val = float(density.std(ddof=0)) if density.size else np.nan
    cv_val = float(std_val / mean_val) if mean_val > 0 else np.nan
    return mean_val, cv_val


def compute_spatial_autocorrelation(adata: ad.AnnData, n_hvg: int) -> float:
    work = adata.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=min(n_hvg, work.n_vars))
    hvgs = work.var_names[work.var["highly_variable"]].tolist()
    if not hvgs:
        return np.nan
    sq.gr.spatial_neighbors(work, coord_type="generic")
    sq.gr.spatial_autocorr(work, mode="moran", genes=hvgs, n_perms=None)
    moran_df = work.uns["moranI"].copy()
    moran_values = pd.to_numeric(moran_df["I"], errors="coerce").dropna()
    return float(moran_values.mean()) if len(moran_values) else np.nan


def build_heartsvg_input(adata: ad.AnnData) -> pd.DataFrame:
    x = to_dense(adata.X).astype(float)
    gene_names = np.asarray(adata.var_names.astype(str))
    spatial = np.asarray(adata.obsm["spatial"], dtype=float)
    df = pd.DataFrame(x, columns=gene_names)
    df.insert(0, "col", spatial[:, 0])
    df.insert(0, "row", spatial[:, 1])
    return df


def compute_svg_number_with_heartsvg(adata: ad.AnnData, rscript_path: str, heartsvg_runner: str, p_adj_threshold: float, pval_threshold: float) -> int:
    if shutil.which(rscript_path) is None:
        raise RuntimeError(f"Rscript executable not found: {rscript_path}")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        input_tsv = tmpdir_path / "heartsvg_input.tsv"
        output_tsv = tmpdir_path / "heartsvg_output.tsv"
        build_heartsvg_input(adata).to_csv(input_tsv, sep="\t", index=False)
        cmd = [rscript_path, heartsvg_runner, str(input_tsv), str(output_tsv), str(p_adj_threshold), str(pval_threshold)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError("HEARTSVG execution failed.\n" f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}")
        res = pd.read_csv(output_tsv, sep="\t")
        if "num_svg" not in res.columns:
            raise RuntimeError("HEARTSVG output is missing the 'num_svg' column")
        return int(res.loc[0, "num_svg"])


def load_uni_model(ckpt_path: str, device: str):
    import timm
    import torch
    model = timm.create_model("vit_large_patch16_224", init_values=1e-5, num_classes=0, dynamic_img_size=True)
    state_dict = torch.load(ckpt_path, map_location="cpu")
    model.load_state_dict(state_dict, strict=True)
    model.eval().to(device)
    return model


def get_uni_transform():
    from torchvision import transforms
    return transforms.Compose([
        transforms.Resize(224, interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])


def load_image_any_format(path: str) -> np.ndarray:
    suffix = Path(path).suffix.lower()
    if suffix in {".tif", ".tiff"}:
        import tifffile
        image = tifffile.imread(path)
    else:
        from PIL import Image
        image = np.asarray(Image.open(path))
    return image


def compute_he_image_feature(he_image_path: str, ckpt_path: str, device: str) -> np.ndarray:
    import torch
    from PIL import Image
    image = load_image_any_format(he_image_path)
    if image.ndim == 2:
        image = np.stack([image] * 3, axis=-1)
    elif image.ndim == 3 and image.shape[-1] >= 4:
        image = image[..., :3]
    if image.dtype != np.uint8:
        image = np.clip(image, 0, 255).astype(np.uint8)
    model = load_uni_model(ckpt_path, device)
    transform = get_uni_transform()
    pil_img = Image.fromarray(image).convert("RGB")
    with torch.no_grad():
        tensor = transform(pil_img).unsqueeze(0).to(device)
        feat = model(tensor).squeeze(0).cpu().numpy().astype(np.float32)
    return feat


def compute_one_dataset(dataset_dir: Path, tissue_map: Dict[str, object], resolution_map: Dict[str, object], uni_ckpt: Optional[str], uni_device: str, n_hvg: int, rscript: str, heartsvg_runner: str, svg_p_adj_threshold: float, svg_pval_threshold: float, logger: logging.Logger) -> Dict[str, object]:
    dataset_name = dataset_dir.name
    logger.info("Processing dataset: %s", dataset_name)
    h5ad_path = find_single_existing_file(dataset_dir, ["*.h5ad"], required=False)
    if h5ad_path is None:
        cnts_path = find_single_existing_file(dataset_dir, ["cnts.tsv"], required=True)
        locs_path = find_single_existing_file(dataset_dir, ["locs-raw.tsv"], required=True)
    else:
        cnts_path = None
        locs_path = None
    image_path = find_single_existing_file(dataset_dir, IMAGE_PATTERNS, required=False)
    adata = read_input_as_adata(h5ad_path, cnts_path, locs_path)

    row: Dict[str, object] = {
        "dataset_name": dataset_name,
        "he_image_feature_json": "",
        "he_image_feature_dim": np.nan,
        "tissue_type": tissue_map.get(dataset_name, ""),
        "spatial_resolution_um_per_pixel": np.nan,
        "data_sparsity": compute_data_sparsity(adata),
        "spatial_autocorrelation": compute_spatial_autocorrelation(adata, n_hvg=n_hvg),
        "mean_expression_per_100um2": np.nan,
        "cv_expression_density": np.nan,
        "svg_number": compute_svg_number_with_heartsvg(adata, rscript, heartsvg_runner, svg_p_adj_threshold, svg_pval_threshold),
    }

    if dataset_name in resolution_map:
        um_per_pixel = float(resolution_map[dataset_name])
        row["spatial_resolution_um_per_pixel"] = um_per_pixel
        area_spot_um2 = estimate_spot_area_um2(adata, um_per_pixel)
        mean_expr, cv_expr = compute_expression_density_stats(adata, area_spot_um2)
        row["mean_expression_per_100um2"] = mean_expr
        row["cv_expression_density"] = cv_expr
    else:
        logger.warning("No spatial resolution provided for dataset '%s'; mean_expression_per_100um2 and cv_expression_density will be empty", dataset_name)

    if image_path is not None:
        if not uni_ckpt:
            logger.warning("HE image found for dataset '%s' but --uni-ckpt not provided; HE feature will be empty", dataset_name)
        else:
            feat = compute_he_image_feature(str(image_path), uni_ckpt, uni_device)
            row["he_image_feature_json"] = json.dumps([float(x) for x in feat])
            row["he_image_feature_dim"] = int(feat.shape[0])
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch-compute feature tables for dataset folders under a root directory")
    parser.add_argument("--dataset-root", required=True, help="Root directory containing one subdirectory per dataset")
    parser.add_argument("--output-csv", required=True, help="Output CSV path. One row per successfully processed dataset")
    parser.add_argument("--tissue-type-json", default=None, help="JSON string or JSON file path mapping dataset_name to tissue type")
    parser.add_argument("--spatial-resolution-json", default=None, help="JSON string or JSON file path mapping dataset_name to spatial resolution in um/pixel")
    parser.add_argument("--uni-ckpt", default=None, help="UNI checkpoint path for HE image feature extraction")
    parser.add_argument("--uni-device", default="cpu", help="Device for UNI inference, e.g. cpu or cuda")
    parser.add_argument("--n-hvg", type=int, default=50, help="Number of HVGs used for Moran's I. Default: 50")
    parser.add_argument("--svg-p-adj-threshold", type=float, default=0.05, help="Adjusted p-value threshold for SVG counting. Default: 0.05")
    parser.add_argument("--svg-pval-threshold", type=float, default=0.05, help="Raw p-value threshold for SVG counting. Default: 0.05")
    parser.add_argument("--heartsvg-runner", default=str(Path(__file__).with_name("run_heartsvg.R")), help="Path to the HEARTSVG R wrapper script")
    parser.add_argument("--rscript", default="Rscript", help="Path to the Rscript executable")
    parser.add_argument("--log-file", default=None, help="Optional log file path")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    args = parser.parse_args()

    logger = setup_logger(args.log_file, args.verbose)
    dataset_dirs = discover_dataset_dirs(Path(args.dataset_root))
    tissue_map = parse_mapping_arg(args.tissue_type_json, "--tissue-type-json")
    resolution_map = parse_mapping_arg(args.spatial_resolution_json, "--spatial-resolution-json")

    rows: List[Dict[str, object]] = []
    failures = 0
    for dataset_dir in dataset_dirs:
        try:
            rows.append(compute_one_dataset(dataset_dir, tissue_map, resolution_map, args.uni_ckpt, args.uni_device, args.n_hvg, args.rscript, args.heartsvg_runner, args.svg_p_adj_threshold, args.svg_pval_threshold, logger))
        except Exception as e:
            failures += 1
            logger.exception("Failed to process dataset '%s': %s", dataset_dir.name, e)

    if not rows:
        raise RuntimeError("No datasets were processed successfully")

    df = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    logger.info("Saved %d successful dataset rows to %s", len(df), output_path)
    if failures:
        logger.warning("Skipped %d datasets due to errors; see logs for details", failures)


if __name__ == "__main__":
    main()
