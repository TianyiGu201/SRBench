#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
extract_single_region_feature.py

Extract a UNI feature from one specified region in an HE image.

Functionality:
1. Load an HE image
2. Crop one specified region using minx, maxx, miny, maxy
3. Extract the UNI feature
4. Save:
   - region_feature.npy
   - region_meta.json
   - region_crop.png (optional)
   - run.log
"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path

import numpy as np
import torch
import timm
import tifffile
from PIL import Image
from torchvision import transforms


def setup_logger(log_path=None, verbose=False):
    logger = logging.getLogger(f"extract_single_region_feature_{id(log_path)}")
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


def main():
    parser = argparse.ArgumentParser(
        description="Extract a UNI feature for one specified HE image region."
    )
    parser.add_argument("--he_image_path", type=str, required=True,
                        help="Path to the HE image")
    parser.add_argument("--minx", type=int, required=True,
                        help="Region min x")
    parser.add_argument("--maxx", type=int, required=True,
                        help="Region max x")
    parser.add_argument("--miny", type=int, required=True,
                        help="Region min y")
    parser.add_argument("--maxy", type=int, required=True,
                        help="Region max y")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Directory to save outputs")
    parser.add_argument("--ckpt_path", type=str, required=True,
                        help="Path to the UNI checkpoint")
    parser.add_argument("--device", type=str, default=None,
                        help="cuda or cpu, default: auto")
    parser.add_argument("--save_crop", action="store_true",
                        help="Save the cropped image region")
    parser.add_argument("--feature_name", type=str, default="region_feature.npy",
                        help="Output feature filename")
    parser.add_argument("--verbose", action="store_true",
                        help="Enable verbose console logging")

    args = parser.parse_args()

    safe_mkdir(args.output_dir)
    logger = setup_logger(os.path.join(args.output_dir, "run.log"), verbose=args.verbose)

    if args.device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    try:
        logger.info(f"Loading HE image: {args.he_image_path}")
        he_img = tifffile.imread(args.he_image_path)
        logger.info(f"HE image loaded, shape={he_img.shape}, dtype={he_img.dtype}")

        crop = crop_he_region(he_img, args.minx, args.maxx, args.miny, args.maxy)
        if crop is None:
            raise ValueError("Invalid or empty crop region")

        logger.info(f"Cropped region shape: {crop.shape}")

        model = load_uni_model(args.ckpt_path, device, logger)
        transform = get_transform()

        pil_img = Image.fromarray(crop.astype(np.uint8)).convert("RGB")
        feature = extract_feature_single(model, pil_img, transform, device)

        feature_path = os.path.join(args.output_dir, args.feature_name)
        np.save(feature_path, feature)
        logger.info(f"Saved feature: {feature_path}, shape={feature.shape}")

        metadata = {
            "he_image_path": args.he_image_path,
            "minx": int(args.minx),
            "maxx": int(args.maxx),
            "miny": int(args.miny),
            "maxy": int(args.maxy),
            "feature_path": feature_path,
            "feature_dim": int(feature.shape[0]),
        }

        meta_path = os.path.join(args.output_dir, "region_meta.json")
        with open(meta_path, "w", encoding="utf-8") as file_obj:
            json.dump(metadata, file_obj, ensure_ascii=False, indent=2)
        logger.info(f"Saved metadata: {meta_path}")

        if args.save_crop:
            crop_path = os.path.join(args.output_dir, "region_crop.png")
            Image.fromarray(crop.astype(np.uint8)).save(crop_path)
            logger.info(f"Saved crop image: {crop_path}")

    except Exception as exc:
        logger.exception(f"Feature extraction failed: {exc}")
        raise


if __name__ == "__main__":
    main()