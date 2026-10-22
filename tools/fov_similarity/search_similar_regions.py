#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
search_similar_regions.py

Search for the top-N most similar regions in an HE image using a reference UNI feature.

Functionality:
1. Load a reference .npy feature
2. Perform sliding-window feature extraction on the HE image
3. Compute cosine similarity
4. Return the top-N most similar regions
5. Save:
   - top_similar_regions.csv
   - top_similar_regions_preview.png
   - crops/top1_crop.png, top2_crop.png, ...
   - region_features/top1_feature.npy, top2_feature.npy, ...
   - run.log
   - run_args.json
"""

import os
import sys
import json
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
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def setup_logger(log_path=None, verbose=False):
    logger = logging.getLogger(f"search_similar_regions_{id(log_path)}")
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


def load_feature_npy(npy_path):
    feature = np.load(npy_path)
    feature = np.asarray(feature, dtype=np.float32)

    if feature.ndim == 2:
        if feature.shape[0] == 1:
            feature = feature.squeeze(0)
        else:
            raise ValueError(f"Expected a single feature vector, got shape={feature.shape}")

    if feature.ndim != 1:
        raise ValueError(f"Feature must be 1D, got shape={feature.shape}")

    return feature


@torch.no_grad()
def extract_feature_batch(model, pil_imgs, transform, device):
    batch = torch.stack([transform(img) for img in pil_imgs]).to(device)
    features = model(batch)
    return features.cpu().numpy().astype(np.float32)


def cosine_similarity(ref, bank):
    ref = ref.reshape(1, -1).astype(np.float32)
    bank = bank.reshape(bank.shape[0], -1).astype(np.float32)

    ref_norm = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-8)
    bank_norm = bank / (np.linalg.norm(bank, axis=1, keepdims=True) + 1e-8)

    return (bank_norm @ ref_norm.T).squeeze(1)


def generate_positions(image_h, image_w, fov_width, fov_height, step):
    positions = []
    for y in range(0, image_h - fov_height + 1, step):
        for x in range(0, image_w - fov_width + 1, step):
            positions.append((x, y, x + fov_width, y + fov_height))
    return positions


def save_preview_image(he_img, result_df, save_path, preview_max_size=1600):
    if he_img.ndim == 2:
        he_rgb = np.stack([he_img] * 3, axis=-1)
    elif he_img.ndim == 3 and he_img.shape[-1] == 4:
        he_rgb = he_img[..., :3]
    else:
        he_rgb = he_img

    height, width = he_rgb.shape[:2]
    scale = min(preview_max_size / max(height, width), 1.0)
    disp_h, disp_w = int(height * scale), int(width * scale)

    preview = Image.fromarray(he_rgb.astype(np.uint8)).resize(
        (disp_w, disp_h), Image.BILINEAR
    )

    fig_w = 10
    fig_h = max(6, fig_w * disp_h / disp_w)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(preview)
    ax.axis("off")

    for _, row in result_df.iterrows():
        x = row["minx"] * scale
        y = row["miny"] * scale
        w = (row["maxx"] - row["minx"]) * scale
        h = (row["maxy"] - row["miny"]) * scale

        rect = patches.Rectangle(
            (x, y), w, h,
            linewidth=2,
            edgecolor="red",
            facecolor="none"
        )
        ax.add_patch(rect)

        ax.text(
            x,
            max(0, y - 5),
            f"#{int(row['rank'])}: {row['similarity']:.4f}",
            fontsize=9,
            color="yellow",
            bbox=dict(facecolor="black", alpha=0.6, pad=2)
        )

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_topn_crops(he_img, result_df, crop_dir):
    safe_mkdir(crop_dir)
    for _, row in result_df.iterrows():
        rank = int(row["rank"])
        minx, miny, maxx, maxy = int(row["minx"]), int(row["miny"]), int(row["maxx"]), int(row["maxy"])
        crop = crop_he_region(he_img, minx, maxx, miny, maxy)
        if crop is None:
            continue
        out_path = os.path.join(crop_dir, f"top{rank}_crop.png")
        Image.fromarray(crop.astype(np.uint8)).save(out_path)


def save_topn_features(result_df, feature_bank, feature_dir):
    safe_mkdir(feature_dir)
    for _, row in result_df.iterrows():
        rank = int(row["rank"])
        bank_index = int(row["bank_index"])
        feature = feature_bank[bank_index]
        out_path = os.path.join(feature_dir, f"top{rank}_feature.npy")
        np.save(out_path, feature.astype(np.float32))


def save_run_args(args, save_path):
    with open(save_path, "w", encoding="utf-8") as file_obj:
        json.dump(vars(args), file_obj, ensure_ascii=False, indent=2)


def find_similar_regions(
    npy_path,
    he_image_path,
    fov_width,
    fov_height,
    step,
    topn,
    model,
    transform,
    device,
    batch_size,
    logger,
):
    ref_feature = load_feature_npy(npy_path)
    logger.info(f"Reference feature loaded: {npy_path}, shape={ref_feature.shape}")

    he_img = load_he_image(he_image_path)
    image_h, image_w = he_img.shape[:2]
    logger.info(f"HE image loaded: {he_image_path}")
    logger.info(f"HE shape = ({image_h}, {image_w})")
    logger.info(f"Sliding window = ({fov_width}, {fov_height}), step={step}")

    positions = generate_positions(image_h, image_w, fov_width, fov_height, step)
    n_windows = len(positions)
    logger.info(f"Total sliding windows: {n_windows}")

    if n_windows == 0:
        logger.warning("No valid windows found. Check window size versus image size.")
        empty_df = pd.DataFrame(columns=[
            "rank", "similarity", "minx", "miny", "maxx", "maxy", "bank_index"
        ])
        return empty_df, he_img, np.empty((0, ref_feature.shape[0]), dtype=np.float32)

    feature_batches = []

    for start in tqdm(range(0, n_windows, batch_size), desc="Sliding-window feature extraction"):
        batch_positions = positions[start:start + batch_size]
        pil_imgs = []

        for minx, miny, maxx, maxy in batch_positions:
            crop = crop_he_region(he_img, minx, maxx, miny, maxy)
            if crop is None:
                crop = np.zeros((fov_height, fov_width, 3), dtype=np.uint8)
            pil_imgs.append(Image.fromarray(crop.astype(np.uint8)).convert("RGB"))

        batch_features = extract_feature_batch(model, pil_imgs, transform, device)
        feature_batches.append(batch_features)

    feature_bank = np.concatenate(feature_batches, axis=0)
    logger.info(f"Feature bank shape: {feature_bank.shape}")

    if ref_feature.shape[0] != feature_bank.shape[1]:
        raise ValueError(
            f"Feature dimension mismatch: ref={ref_feature.shape[0]}, bank={feature_bank.shape[1]}"
        )

    similarities = cosine_similarity(ref_feature, feature_bank)
    topn = min(topn, n_windows)
    top_indices = np.argsort(similarities)[::-1][:topn]

    rows = []
    for rank, idx in enumerate(top_indices, start=1):
        minx, miny, maxx, maxy = positions[idx]
        rows.append({
            "rank": rank,
            "similarity": float(similarities[idx]),
            "minx": int(minx),
            "miny": int(miny),
            "maxx": int(maxx),
            "maxy": int(maxy),
            "bank_index": int(idx),
        })

    result_df = pd.DataFrame(rows)
    logger.info("Top-N results:\n" + result_df.to_string(index=False))

    return result_df, he_img, feature_bank


def main():
    parser = argparse.ArgumentParser(
        description="Search top-N similar regions in an HE image using a reference UNI feature."
    )
    parser.add_argument("--npy_path", type=str, required=True,
                        help="Reference feature .npy path")
    parser.add_argument("--he_image_path", type=str, required=True,
                        help="HE image path")
    parser.add_argument("--fov_width", type=int, required=True,
                        help="Sliding window width in pixels")
    parser.add_argument("--fov_height", type=int, required=True,
                        help="Sliding window height in pixels")
    parser.add_argument("--step", type=int, default=100,
                        help="Sliding step in pixels")
    parser.add_argument("--topn", type=int, default=5,
                        help="Number of top similar regions to return")
    parser.add_argument("--batch_size", type=int, default=64,
                        help="Batch size for feature extraction")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save outputs")
    parser.add_argument("--ckpt_path", type=str, required=True,
                        help="Path to the UNI checkpoint")
    parser.add_argument("--device", type=str, default=None,
                        help="cuda or cpu, default: auto")
    parser.add_argument("--preview_max_size", type=int, default=1600,
                        help="Maximum edge length for the preview image")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose console logging")

    args = parser.parse_args()

    safe_mkdir(args.output_dir)
    logger = setup_logger(os.path.join(args.output_dir, "run.log"), verbose=args.verbose)
    save_run_args(args, os.path.join(args.output_dir, "run_args.json"))

    if args.device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    try:
        model = load_uni_model(args.ckpt_path, device, logger)
        transform = get_transform()

        result_df, he_img, feature_bank = find_similar_regions(
            npy_path=args.npy_path,
            he_image_path=args.he_image_path,
            fov_width=args.fov_width,
            fov_height=args.fov_height,
            step=args.step,
            topn=args.topn,
            model=model,
            transform=transform,
            device=device,
            batch_size=args.batch_size,
            logger=logger,
        )

        csv_df = result_df.drop(columns=["bank_index"], errors="ignore")
        csv_path = os.path.join(args.output_dir, "top_similar_regions.csv")
        csv_df.to_csv(csv_path, index=False)
        logger.info(f"Saved CSV: {csv_path}")

        preview_path = os.path.join(args.output_dir, "top_similar_regions_preview.png")
        save_preview_image(
            he_img=he_img,
            result_df=result_df,
            save_path=preview_path,
            preview_max_size=args.preview_max_size,
        )
        logger.info(f"Saved preview image: {preview_path}")

        crop_dir = os.path.join(args.output_dir, "crops")
        save_topn_crops(he_img, result_df, crop_dir)
        logger.info(f"Saved crop images to: {crop_dir}")

        feature_dir = os.path.join(args.output_dir, "region_features")
        save_topn_features(result_df, feature_bank, feature_dir)
        logger.info(f"Saved top-N region features to: {feature_dir}")

    except Exception as exc:
        logger.exception(f"Search failed: {exc}")
        raise


if __name__ == "__main__":
    main()