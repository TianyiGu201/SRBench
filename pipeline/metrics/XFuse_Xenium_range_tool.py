#!/usr/bin/env python3

import argparse
import gzip
import os
import pickle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

if __package__:
    from .spatial_eval_common import (
        ensure_dir,
        find_pred_file,
        load_pickle_any,
    )
else:
    from spatial_eval_common import (
        ensure_dir,
        find_pred_file,
        load_pickle_any,
    )


XFUSE_DEFAULT_PRED_SUBDIR = os.path.join(
    "my-run", "analyses", "final", "analysis-gene_maps", "section0"
)


def save_pickle_like(src_path, obj, dst_path):
    ensure_dir(os.path.dirname(dst_path))
    if dst_path.endswith(".pkl.gz"):
        with gzip.open(dst_path, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(dst_path, "wb") as f:
            pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def apply_transform_and_crop(arr, do_transpose=False, do_fliplr=False, top=None, bottom=None, left=None, right=None):
    out = np.asarray(arr)

    if do_transpose:
        out = out.T
    if do_fliplr:
        out = np.fliplr(out)

    h, w = out.shape[:2]

    t = 0 if top is None else int(top)
    b = h if bottom is None else int(bottom)
    l = 0 if left is None else int(left)
    r = w if right is None else int(right)

    if not (0 <= t <= b <= h):
        raise ValueError(f"Invalid crop rows: top={t}, bottom={b}, valid range is [0, {h}]")
    if not (0 <= l <= r <= w):
        raise ValueError(f"Invalid crop cols: left={l}, right={r}, valid range is [0, {w}]")

    out = out[t:b, l:r]
    return np.asarray(out)


def build_param_text(args_or_dict):
    if isinstance(args_or_dict, dict):
        transpose = bool(args_or_dict.get("transpose", False))
        fliplr = bool(args_or_dict.get("fliplr", False))
        top = args_or_dict.get("top", None)
        bottom = args_or_dict.get("bottom", None)
        left = args_or_dict.get("left", None)
        right = args_or_dict.get("right", None)
    else:
        transpose = bool(args_or_dict.transpose)
        fliplr = bool(args_or_dict.fliplr)
        top = args_or_dict.top
        bottom = args_or_dict.bottom
        left = args_or_dict.left
        right = args_or_dict.right

    return (
        f"transpose={transpose}\n"
        f"fliplr={fliplr}\n"
        f"crop(top,bottom,left,right)=({top},{bottom},{left},{right})"
    )


def load_he_image_rgb(path):
    with Image.open(path) as im:
        return np.asarray(im.convert("RGB"))


def preview_gene(args):
    pred_dir = os.path.abspath(args.pred_dir)
    he_path = os.path.abspath(args.he_path)
    pred_path = find_pred_file(pred_dir, args.gene)
    if pred_path is None:
        raise FileNotFoundError(f"Missing gene file for {args.gene} in {pred_dir}")

    arr = load_pickle_any(pred_path)
    arr = apply_transform_and_crop(
        arr,
        do_transpose=args.transpose,
        do_fliplr=args.fliplr,
        top=args.top,
        bottom=args.bottom,
        left=args.left,
        right=args.right,
    )

    he_img = load_he_image_rgb(he_path)

    out_png = os.path.abspath(args.out_png)
    ensure_dir(os.path.dirname(out_png))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].imshow(he_img)
    axes[0].set_title("Reference HE image")
    axes[0].axis("off")

    axes[1].imshow(arr, origin="upper")
    axes[1].set_title(f"{args.gene}\n{build_param_text(args)}")
    axes[1].axis("off")

    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"[DONE] preview saved: {out_png}")


def is_prediction_dir(path):
    if not os.path.isdir(path):
        return False
    for name in os.listdir(path):
        if name.endswith((".pkl.gz", ".pkl", ".pickle")):
            return True
    return False


def is_sample_dir_for_xfuse(path, pred_subdir):
    if not os.path.isdir(path):
        return False
    pred_dir = os.path.join(path, pred_subdir)
    return os.path.isdir(pred_dir) and is_prediction_dir(pred_dir)


def collect_target_pred_dirs(prediction_path, pred_subdir=XFUSE_DEFAULT_PRED_SUBDIR, recursive=True):
    prediction_path = os.path.abspath(prediction_path)

    if not os.path.exists(prediction_path):
        raise FileNotFoundError(f"Path not found: {prediction_path}")

    if os.path.isfile(prediction_path):
        raise ValueError(f"prediction-path must be a directory, got file: {prediction_path}")


    if is_prediction_dir(prediction_path):
        return [prediction_path], True, prediction_path


    if is_sample_dir_for_xfuse(prediction_path, pred_subdir):
        return [os.path.join(prediction_path, pred_subdir)], True, prediction_path


    pred_dirs = []
    seen = set()

    walker = os.walk(prediction_path) if recursive else [(prediction_path, [], [])]
    for root, _, _ in walker:
        if is_sample_dir_for_xfuse(root, pred_subdir):
            pred_dir = os.path.join(root, pred_subdir)
            if pred_dir not in seen:
                pred_dirs.append(pred_dir)
                seen.add(pred_dir)

    pred_dirs.sort()
    return pred_dirs, False, prediction_path


def iter_prediction_files(pred_dir):
    files = []
    for name in sorted(os.listdir(pred_dir)):
        if name.endswith((".pkl.gz", ".pkl", ".pickle")):
            files.append(os.path.join(pred_dir, name))
    return files


def build_output_pred_dir(src_pred_dir, search_root, output_root):
    src_pred_dir = os.path.abspath(src_pred_dir)
    search_root = os.path.abspath(search_root)
    output_root = os.path.abspath(output_root)

    rel = os.path.relpath(src_pred_dir, search_root)
    return os.path.join(output_root, rel)


def process_one_pred_dir(pred_dir, dst_pred_dir, args):
    files = iter_prediction_files(pred_dir)
    if not files:
        print(f"[WARN] no prediction files found: {pred_dir}")
        return 0

    ensure_dir(dst_pred_dir)

    count = 0
    for src_path in files:
        arr = load_pickle_any(src_path)
        out = apply_transform_and_crop(
            arr,
            do_transpose=args.transpose,
            do_fliplr=args.fliplr,
            top=args.top,
            bottom=args.bottom,
            left=args.left,
            right=args.right,
        )

        if args.replace_src:
            dst_path = src_path
        else:
            dst_path = os.path.join(dst_pred_dir, os.path.basename(src_path))

        save_pickle_like(src_path, out, dst_path)
        count += 1

    return count


def batch_process(args):
    pred_dirs, is_single, search_root = collect_target_pred_dirs(
        prediction_path=args.prediction_path,
        pred_subdir=args.pred_subdir,
        recursive=not args.no_recursive,
    )

    if not pred_dirs:
        raise RuntimeError("No valid XFuse prediction directories found")

    if not args.replace_src and not args.output_root:
        raise ValueError("When --replace-src is not set, --output-root is required")

    print("[INFO] scheme:")
    print(build_param_text(args))
    print(f"[INFO] found prediction dirs: {len(pred_dirs)}")

    total_files = 0

    for idx, pred_dir in enumerate(pred_dirs, start=1):
        if args.replace_src:
            dst_pred_dir = pred_dir
        else:
            dst_pred_dir = build_output_pred_dir(pred_dir, search_root, args.output_root)

        n = process_one_pred_dir(pred_dir, dst_pred_dir, args)
        total_files += n

        print(f"[{idx}/{len(pred_dirs)}] src={pred_dir}")
        if args.replace_src:
            print(f"          mode=replace-src")
        else:
            print(f"          dst={dst_pred_dir}")
        print(f"          processed_files={n}")

    print(f"[DONE] total prediction dirs: {len(pred_dirs)}")
    print(f"[DONE] total processed files: {total_files}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="XFuse Xenium tool: preview or batch transform prediction files"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p1 = sub.add_parser("preview", help="Preview one gene after transpose/fliplr/crop")
    p1.add_argument("--pred-dir", required=True, help="Prediction directory containing gene pickle files")
    p1.add_argument("--he-path", required=True, help="Reference HE image path")
    p1.add_argument("--gene", required=True, help="Gene name")
    p1.add_argument("--out-png", required=True, help="Output preview png path")
    p1.add_argument("--transpose", action="store_true", help="Apply transpose before fliplr/crop")
    p1.add_argument("--fliplr", action="store_true", help="Apply horizontal flip after transpose")
    p1.add_argument("--top", type=int, default=None, help="Crop top boundary after transform")
    p1.add_argument("--bottom", type=int, default=None, help="Crop bottom boundary after transform")
    p1.add_argument("--left", type=int, default=None, help="Crop left boundary after transform")
    p1.add_argument("--right", type=int, default=None, help="Crop right boundary after transform")

    p2 = sub.add_parser("batch-process", help="Batch apply transpose/fliplr/crop to prediction files")
    p2.add_argument(
        "--prediction-path",
        required=True,
        help=(
            "Can be: "
            "(1) one prediction dir containing .pkl/.pkl.gz/.pickle files, or "
            "(2) one sample dir, or "
            "(3) one root dir to recursively search sample dirs"
        ),
    )
    p2.add_argument(
        "--pred-subdir",
        default=XFUSE_DEFAULT_PRED_SUBDIR,
        help=f"Relative prediction subdir under each sample dir (default: {XFUSE_DEFAULT_PRED_SUBDIR})",
    )
    p2.add_argument("--no-recursive", action="store_true", help="Do not recursively search under prediction-path")
    p2.add_argument("--transpose", action="store_true", help="Apply transpose before fliplr/crop")
    p2.add_argument("--fliplr", action="store_true", help="Apply horizontal flip after transpose")
    p2.add_argument("--top", type=int, default=None, help="Crop top boundary after transform")
    p2.add_argument("--bottom", type=int, default=None, help="Crop bottom boundary after transform")
    p2.add_argument("--left", type=int, default=None, help="Crop left boundary after transform")
    p2.add_argument("--right", type=int, default=None, help="Crop right boundary after transform")

    mode = p2.add_mutually_exclusive_group(required=True)
    mode.add_argument("--replace-src", action="store_true", help="Replace source files in place")
    mode.add_argument("--output-root", default=None, help="Save processed files to a new root directory")

    return parser


def main():
    args = build_parser().parse_args()

    if args.command == "preview":
        preview_gene(args)
    elif args.command == "batch-process":
        batch_process(args)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()