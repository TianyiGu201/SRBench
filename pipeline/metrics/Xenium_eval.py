#!/usr/bin/env python3

import argparse
import csv
import os

import numpy as np

if __package__:
    from .spatial_eval_common import (
        COLUMN_NAMES,
        EPS,
        algorithm_key,
        build_gt_gene_coords,
        build_gt_hist,
        build_unique_names,
        collect_sample_dirs,
        dump_pkl_gz,
        ensure_dir,
        evaluate_pair,
        find_pred_file,
        get_algorithm_pred_subdirs,
        is_xenium_sample_dir,
        istar_scale_pred,
        load_he_hw,
        load_pickle_any,
        load_transcripts_tsv,
        make_gt_cache_map,
        make_output_dir_map,
        read_lines,
        resolve_default_gt_cache_dir,
        should_filter_gene,
    )
else:
    from spatial_eval_common import (
        COLUMN_NAMES,
        EPS,
        algorithm_key,
        build_gt_gene_coords,
        build_gt_hist,
        build_unique_names,
        collect_sample_dirs,
        dump_pkl_gz,
        ensure_dir,
        evaluate_pair,
        find_pred_file,
        get_algorithm_pred_subdirs,
        is_xenium_sample_dir,
        istar_scale_pred,
        load_he_hw,
        load_pickle_any,
        load_transcripts_tsv,
        make_gt_cache_map,
        make_output_dir_map,
        read_lines,
        resolve_default_gt_cache_dir,
        should_filter_gene,
    )


def write_csv(path, rows):
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(COLUMN_NAMES)
        writer.writerows(rows)


def evaluate_one_head(sample_dir, head_name, pred_subdir, out_dir, gt_cache_dir, gt_by_gene, he_h, he_w, args):
    pred_dir = os.path.join(sample_dir, pred_subdir)
    if not os.path.isdir(pred_dir):
        print(f"[WARN] missing pred dir: {pred_dir}")
        return

    out_csv = os.path.join(out_dir, f"{args.name_map[sample_dir]}.csv")
    if os.path.exists(out_csv) and not args.force:
        print(f"[SKIP] output exists: {out_csv}")
        return

    gene_file = os.path.join(sample_dir, "gene-names.txt")
    if not os.path.isfile(gene_file):
        print(f"[WARN] missing: {gene_file}")
        return

    rows = []
    missing = 0

    for gene in read_lines(gene_file):
        if should_filter_gene(gene):
            continue

        pred_path = find_pred_file(pred_dir, gene)
        if pred_path is None:
            missing += 1
            continue

        pred = np.asarray(load_pickle_any(pred_path), dtype=np.float64)

        if args.algorithm == "istar":
            pred = istar_scale_pred(pred, sample_dir=sample_dir)

        if pred.ndim != 2:
            print(f"[WARN] prediction is not 2D for {gene}: shape={getattr(pred, 'shape', None)}")
            continue

        pd_h, pd_w = pred.shape
        gt_shape_dir = os.path.join(gt_cache_dir, f"{pd_h}x{pd_w}")
        gt_cache_path = os.path.join(gt_shape_dir, f"{gene}.pkl.gz")

        if os.path.exists(gt_cache_path):
            y_true = load_pickle_any(gt_cache_path)
        else:
            coords = gt_by_gene.get(gene)
            y_true = build_gt_hist(coords, pd_h, pd_w, he_h, he_w)
            dump_pkl_gz(y_true, gt_cache_path)

        y_true = np.asarray(y_true, dtype=np.float64)

        if y_true.shape != pred.shape:
            print(f"[WARN] shape mismatch {gene}: gt={y_true.shape}, pred={pred.shape}")
            continue

        rows.append([gene] + evaluate_pair(y_true, pred, eps=args.eps))

    if not rows:
        print(f"[WARN] no valid genes: {sample_dir} ({head_name}), missing={missing}")
        return

    write_csv(out_csv, rows)
    print(f"[DONE] saved: {out_csv}")


def build_parser():
    parser = argparse.ArgumentParser(description="Xenium evaluation")
    parser.add_argument("--he-path", required=True, help="Full-res HE image path")
    parser.add_argument("--transcript-use-path", required=True, help="transcripts_use.tsv path")
    parser.add_argument("--algorithm", required=True, help="iSTAR, scstGCN, TESLA, XFuse, or SuperST")
    parser.add_argument("--output-path", required=True, help="Output directory")
    parser.add_argument("--prediction-path", required=True, help="One sample dir or a parent dir")
    parser.add_argument("--gt-cache-dir", default=None, help="GT cache root")
    parser.add_argument("--feature-col", default="feature_name")
    parser.add_argument("--y-col-idx", type=int, default=2)
    parser.add_argument("--x-col-idx", type=int, default=1)
    parser.add_argument("--no-recursive", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--eps", type=float, default=EPS)
    return parser


def main():
    args = build_parser().parse_args()
    args.algorithm = algorithm_key(args.algorithm)

    he_path = os.path.abspath(args.he_path)
    transcript_use_path = os.path.abspath(args.transcript_use_path)
    prediction_path = os.path.abspath(args.prediction_path)

    gt_rows = load_transcripts_tsv(transcript_use_path, feature_col=args.feature_col)
    gt_by_gene = build_gt_gene_coords(
        gt_rows,
        feature_col=args.feature_col,
        y_col_idx=args.y_col_idx,
        x_col_idx=args.x_col_idx,
    )
    he_h, he_w = load_he_hw(he_path)

    pred_subdirs = list(get_algorithm_pred_subdirs(args.algorithm).values())
    sample_dirs, is_single_sample = collect_sample_dirs(
        prediction_path,
        sample_check=lambda p: is_xenium_sample_dir(p, pred_subdirs),
        recursive=not args.no_recursive,
    )
    if not sample_dirs:
        raise RuntimeError("No valid sample directories found")

    gt_cache_root = args.gt_cache_dir or resolve_default_gt_cache_dir(
        prediction_path,
        is_single_sample,
        suffix="gt_cache",
    )
    gt_cache_map = make_gt_cache_map(gt_cache_root, args.algorithm)
    out_dir_map = make_output_dir_map(os.path.abspath(args.output_path), args.algorithm)
    args.name_map = build_unique_names(sample_dirs)

    for head_name, pred_subdir in get_algorithm_pred_subdirs(args.algorithm).items():
        for sample_dir in sample_dirs:
            evaluate_one_head(
                sample_dir=sample_dir,
                head_name=head_name,
                pred_subdir=pred_subdir,
                out_dir=out_dir_map[head_name],
                gt_cache_dir=gt_cache_map[head_name],
                gt_by_gene=gt_by_gene,
                he_h=he_h,
                he_w=he_w,
                args=args,
            )


if __name__ == "__main__":
    main()