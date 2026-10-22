#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
extract_top3_uni_features.py

Batch extraction of UNI features for the top-3 FOV regions of each dataset/method.

Functionality:
1. Traverse base_csv_dir/{dataset}/{method}/fov_gene_metrics.csv
2. Select the top-3 FOV regions with the smallest mean PARnorm
3. Crop the corresponding regions from the HE image
4. Extract UNI features
5. Save each feature as .npy
6. Save feature_index.csv
7. Save run logs

Output example:
output_base_dir/
©¸©¤©¤ datasetA/
    ©À©¤©¤ feature_index.csv
    ©À©¤©¤ datasetA_method1_top1_fov1500_region12.npy
    ©À©¤©¤ datasetA_method1_top2_fov1500_region08.npy
    ©À©¤©¤ datasetA_method1_top3_fov1500_region31.npy
    ©¸©¤©¤ run.log
"""

import os
import sys
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import timm
import tifffile
from PIL import Image
from torchvision import transforms


def setup_logger(log_path=None, verbose=False):
    logger = logging.getLogger(f"extract_top3_uni_features_{id(log_path)}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if log_path is not None:
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def safe_mkdir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def load_uni_model(ckpt_path, device, logger):
    logger.info(f"Loading UNI model from: {ckpt_path}")
    model = timm.create_model(
        "vit_large_patch16_224",
        init_values=1e-5,
        num_classes=0,
        dynamic_img_size=True,
    )
    state_dict = torch.load(ckpt_path, map_location="cpu")
    model.load_state_dict(state_dict, strict=True)
    model.eval().to(device)
    logger.info(f"UNI model loaded on device: {device}")
    return model


def get_transform():
    return transforms.Compose([
        transforms.Resize(224, interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.CenterCrop(224),
        transforms.ToTensor(),
        transforms.Normalize(
            mean=[0.485, 0.456, 0.406],
            std=[0.229, 0.224, 0.225],
        ),
    ])


def load_he_image(he_path):
    return tifffile.imread(he_path)


def crop_he_region(he_img, minx, maxx, miny, maxy):
    minx, maxx = int(round(minx)), int(round(maxx))
    miny, maxy = int(round(miny)), int(round(maxy))

    height, width = he_img.shape[:2]
    minx, maxx = max(0, minx), min(width, maxx)
    miny, maxy = max(0, miny), min(height, maxy)

    if minx >= maxx or miny >= maxy:
        return None

    crop = he_img[miny:maxy, minx:maxx]

    if crop.size == 0:
        return None

    if crop.ndim == 2:
        crop = np.stack([crop] * 3, axis=-1)
    elif crop.ndim == 3 and crop.shape[-1] == 4:
        crop = crop[..., :3]

    return crop


@torch.no_grad()
def extract_feature_single(model, pil_img, transform, device):
    tensor = transform(pil_img).unsqueeze(0).to(device)
    feature = model(tensor)
    return feature.squeeze(0).cpu().numpy().astype(np.float32)


def find_top3_fovs(csv_path, fov_size_target="fov1500"):
    df = pd.read_csv(csv_path)

    required_cols = ["FOV_region", "FOV_size", "PARnorm", "minx", "maxx", "miny", "maxy"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in {csv_path}: {missing_cols}")

    def normalize_size(value):
        return float(str(value).strip().lower().replace("fov", "").replace("_", ""))
    size = normalize_size(fov_size_target)
    sizes = df["FOV_size"].map(normalize_size)
    df = df[sizes == size].copy()
    df["PARnorm"] = pd.to_numeric(df["PARnorm"], errors="coerce")
    df.loc[~np.isfinite(df["PARnorm"]), "PARnorm"] = np.nan
    if len(df) == 0:
        return pd.DataFrame()

    fov_agg = df.groupby("FOV_region").agg(
        mean_PARnorm=("PARnorm", "mean"),
        FOV_size=("FOV_size", "first"),
        minx=("minx", "first"),
        maxx=("maxx", "first"),
        miny=("miny", "first"),
        maxy=("maxy", "first"),
    ).reset_index()

    fov_agg = fov_agg.loc[np.isfinite(fov_agg["mean_PARnorm"])].copy()
    fov_agg = fov_agg.sort_values("mean_PARnorm", ascending=True).head(3).reset_index(drop=True)
    return fov_agg


def extract_top3_fov_features(
    base_csv_dir,
    he_base_dir,
    output_base_dir,
    ckpt_path,
    device,
    fov_size_target="fov1500",
    verbose=False,
    csv_filename="fov_gene_metrics.csv",
):
    safe_mkdir(output_base_dir)
    root_log_path = os.path.join(output_base_dir, "run.log")
    logger = setup_logger(root_log_path, verbose=verbose)

    logger.info("Starting top-3 UNI feature extraction")
    logger.info(f"base_csv_dir    : {base_csv_dir}")
    logger.info(f"he_base_dir     : {he_base_dir}")
    logger.info(f"output_base_dir : {output_base_dir}")
    logger.info(f"fov_size_target : {fov_size_target}")

    model = load_uni_model(ckpt_path, device, logger)
    transform = get_transform()

    datasets = sorted([
        d for d in os.listdir(base_csv_dir)
        if os.path.isdir(os.path.join(base_csv_dir, d))
    ])

    logger.info(f"Found datasets: {datasets}")

    if not datasets:
        logger.warning("No datasets found. Exiting.")
        return

    for dataset in datasets:
        dataset_csv_dir = os.path.join(base_csv_dir, dataset)
        he_path = os.path.join(he_base_dir, dataset, "he-raw.tiff")
        out_dir = os.path.join(output_base_dir, dataset)
        safe_mkdir(out_dir)

        dataset_logger = setup_logger(os.path.join(out_dir, "run.log"), verbose=verbose)
        dataset_logger.info("=" * 80)
        dataset_logger.info(f"Dataset: {dataset}")
        dataset_logger.info(f"CSV directory: {dataset_csv_dir}")
        dataset_logger.info(f"HE image path: {he_path}")

        if not os.path.exists(he_path):
            dataset_logger.warning(f"HE image not found: {he_path}")
            continue

        try:
            he_img = load_he_image(he_path)
            dataset_logger.info(f"HE image loaded, shape={he_img.shape}, dtype={he_img.dtype}")
        except Exception as exc:
            dataset_logger.exception(f"Failed to load HE image: {exc}")
            continue

        methods = sorted([
            m for m in os.listdir(dataset_csv_dir)
            if os.path.isdir(os.path.join(dataset_csv_dir, m))
        ])

        if not methods:
            dataset_logger.warning("No method folders found for this dataset.")
            continue

        index_rows = []

        for method in methods:
            csv_path = os.path.join(dataset_csv_dir, method, csv_filename)

            if not os.path.exists(csv_path):
                dataset_logger.warning(f"[{method}] Missing CSV file: {csv_path}")
                continue

            try:
                top3 = find_top3_fovs(csv_path, fov_size_target=fov_size_target)
            except Exception as exc:
                dataset_logger.exception(f"[{method}] Failed to parse CSV: {exc}")
                continue

            if top3.empty:
                dataset_logger.warning(f"[{method}] No valid rows for FOV_size={fov_size_target}")
                continue

            dataset_logger.info(f"[{method}] Selected {len(top3)} valid FOV regions (up to 3):")
            dataset_logger.info("\n" + top3[["FOV_region", "mean_PARnorm"]].to_string(index=False))

            for rank_idx, row in top3.iterrows():
                topk = rank_idx + 1
                fov_region = str(row["FOV_region"])
                fov_size = str(row["FOV_size"])

                try:
                    crop = crop_he_region(
                        he_img,
                        row["minx"], row["maxx"],
                        row["miny"], row["maxy"]
                    )

                    if crop is None:
                        dataset_logger.warning(f"[{method}] Invalid crop for region={fov_region}")
                        continue

                    pil_img = Image.fromarray(crop.astype(np.uint8)).convert("RGB")
                    feature = extract_feature_single(model, pil_img, transform, device)

                    filename = f"{dataset}_{method}_top{topk}_{fov_size}_{fov_region}.npy"
                    feature_path = os.path.join(out_dir, filename)
                    np.save(feature_path, feature)

                    index_rows.append({
                        "dataset": dataset,
                        "method": method,
                        "rank": topk,
                        "feature_file": filename,
                        "feature_path": feature_path,
                        "fov_region": fov_region,
                        "fov_size": fov_size,
                        "mean_PARnorm": float(row["mean_PARnorm"]),
                        "minx": int(round(row["minx"])),
                        "maxx": int(round(row["maxx"])),
                        "miny": int(round(row["miny"])),
                        "maxy": int(round(row["maxy"])),
                        "feature_dim": int(feature.shape[0]),
                    })

                    dataset_logger.info(f"[{method}] Saved feature: {filename}, shape={feature.shape}")

                except Exception as exc:
                    dataset_logger.exception(
                        f"[{method}] Failed on region={fov_region}, rank={topk}: {exc}"
                    )
                    continue

        if index_rows:
            index_df = pd.DataFrame(index_rows)
            index_csv_path = os.path.join(out_dir, "feature_index.csv")
            index_df.to_csv(index_csv_path, index=False)
            dataset_logger.info(f"Saved feature index: {index_csv_path}")
        else:
            dataset_logger.warning("No features were extracted for this dataset.")

        del he_img
        dataset_logger.info(f"Finished dataset: {dataset}")

    logger.info("All datasets processed successfully")


def main():
    parser = argparse.ArgumentParser(
        description="Batch extraction of UNI features for the top-3 FOV regions."
    )
    parser.add_argument("--base_csv_dir", type=str, required=True,
                        help="Path containing {dataset}/{method}/<csv_filename>")
    parser.add_argument("--csv_filename", default="fov_gene_metrics.csv",
                        help="Metric filename in each dataset/method folder; use fov_all.csv for older exports")
    parser.add_argument("--he_base_dir", type=str, required=True,
                        help="Path like he_base_dir/{dataset}/he-raw.tiff")
    parser.add_argument("--output_base_dir", type=str, required=True,
                        help="Directory to save feature files and logs")
    parser.add_argument("--ckpt_path", type=str, required=True,
                        help="Path to the UNI checkpoint")
    parser.add_argument("--fov_size_target", type=str, default="fov1500",
                        help="Target FOV_size value, default: fov1500")
    parser.add_argument("--device", type=str, default=None,
                        help="cuda or cpu, default: auto")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose console logging")

    args = parser.parse_args()

    if args.device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    extract_top3_fov_features(
        base_csv_dir=args.base_csv_dir,
        he_base_dir=args.he_base_dir,
        output_base_dir=args.output_base_dir,
        ckpt_path=args.ckpt_path,
        device=device,
        fov_size_target=args.fov_size_target,
        verbose=args.verbose,
        csv_filename=args.csv_filename,
    )


if __name__ == "__main__":
    main()