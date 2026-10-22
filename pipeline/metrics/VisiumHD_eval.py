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
        build_unique_names,
        build_visiumhd_gt_cache_name,
        collect_sample_dirs,
        crop_array_by_range,
        dump_pkl_gz,
        ensure_dir,
        evaluate_pair,
        find_pred_file,
        gene_grid_from_anndata,
        get_algorithm_pred_subdirs,
        infer_visiumhd_range_info,
        is_visiumhd_sample_dir,
        load_json,
        load_pickle_any,
        load_pred_array_for_range_detection,
        make_gt_cache_map,
        make_output_dir_map,
        read_lines,
        resolve_default_gt_cache_dir,
        save_json,
        should_filter_gene,
    )
else:
    from spatial_eval_common import (
        COLUMN_NAMES,
        EPS,
        algorithm_key,
        build_unique_names,
        build_visiumhd_gt_cache_name,
        collect_sample_dirs,
        crop_array_by_range,
        dump_pkl_gz,
        ensure_dir,
        evaluate_pair,
        find_pred_file,
        gene_grid_from_anndata,
        get_algorithm_pred_subdirs,
        infer_visiumhd_range_info,
        is_visiumhd_sample_dir,
        load_json,
        load_pickle_any,
        load_pred_array_for_range_detection,
        make_gt_cache_map,
        make_output_dir_map,
        read_lines,
        resolve_default_gt_cache_dir,
        save_json,
        should_filter_gene,
    )


def write_csv(path, rows):
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(COLUMN_NAMES)
        writer.writerows(rows)


def detect_range(args):
    import pyarrow.parquet as pq
    import scanpy as sc
    gt_data = sc.read_10x_h5(os.path.abspath(args.gt_data))
    tissue_positions = pq.ParquetFile(os.path.abspath(args.tissue_positions)).read().to_pandas()

    if args.pd_file:
        pd_path = os.path.abspath(args.pd_file)
    else:
        pred_dir = os.path.abspath(args.pred_dir)
        pd_path = find_pred_file(pred_dir, args.gene)
        if pd_path is None:
            raise FileNotFoundError(f"Cannot find prediction file for gene={args.gene} in {pred_dir}")

    pd_array = load_pred_array_for_range_detection(pd_path)
    info = infer_visiumhd_range_info(
        pd_array=pd_array,
        gt_data=gt_data,
        tissue_positions=tissue_positions,
        gene=args.gene,
    )
    info["pd_source"] = pd_path
    info["auxiliary_check_only"] = True
    info["manual_review_required"] = True
    save_json(info, os.path.abspath(args.out_range_json))
    print(f"[DONE] saved range json: {os.path.abspath(args.out_range_json)}")
    print(info)
    print("[CHECK] Auxiliary ranges only. Inspect spatial alignment and finalize ranges manually before evaluation.")


def evaluate_one_head(sample_dir, head_name, pred_subdir, out_dir, gt_cache_dir, gt_data, tissue_positions, range_info, args):
    pred_dir = os.path.join(sample_dir, pred_subdir)
    if not os.path.isdir(pred_dir):
        print(f"[WARN] missing pred dir: {pred_dir}")
        return

    out_csv = os.path.join(out_dir, f"{args.name_map[sample_dir]}.csv")
    if os.path.exists(out_csv) and not args.force:
        print(f"[SKIP] output exists: {out_csv}")
        return

    gt_row_range = tuple(range_info["gt_row_range"])
    gt_col_range = tuple(range_info["gt_col_range"])
    pd_row_range = tuple(range_info["pd_row_range"])
    pd_col_range = tuple(range_info["pd_col_range"])
    flip_rows = bool(range_info.get("flip_rows", False))
    flip_cols = bool(range_info.get("flip_cols", False))

    gt_shape_dir = os.path.join(gt_cache_dir, build_visiumhd_gt_cache_name(gt_row_range, gt_col_range))

    rows = []
    gene_file = os.path.join(sample_dir, "gene-names.txt")

    for gene in read_lines(gene_file):
        if should_filter_gene(gene):
            continue

        pred_path = find_pred_file(pred_dir, gene)
        if pred_path is None:
            continue

        pred = np.asarray(load_pickle_any(pred_path), dtype=np.float64)
        pred = crop_array_by_range(pred, row_range=pd_row_range, col_range=pd_col_range)

        gt_cache_path = os.path.join(gt_shape_dir, f"{gene}.pkl.gz")
        if os.path.exists(gt_cache_path):
            y_true = load_pickle_any(gt_cache_path)
        else:
            try:
                y_true = gene_grid_from_anndata(
                    gene=gene,
                    gt_data=gt_data,
                    tissue_positions=tissue_positions,
                    gt_row_range=gt_row_range,
                    gt_col_range=gt_col_range,
                    flip_rows=flip_rows,
                    flip_cols=flip_cols,
                )
            except KeyError:
                continue
            dump_pkl_gz(y_true, gt_cache_path)

        y_true = np.asarray(y_true, dtype=np.float64)

        if y_true.shape != pred.shape:
            print(f"[WARN] shape mismatch {gene}: gt={y_true.shape}, pred={pred.shape}")
            continue

        rows.append([gene] + evaluate_pair(y_true, pred, eps=args.eps))

    if not rows:
        print(f"[WARN] no valid genes: {sample_dir} ({head_name})")
        return

    write_csv(out_csv, rows)
    print(f"[DONE] saved: {out_csv}")


def build_parser():
    parser = argparse.ArgumentParser(description="VisiumHD evaluation")
    sub = parser.add_subparsers(dest="command", required=True)

    p1 = sub.add_parser("detect-range", help="Generate auxiliary candidate ranges for subsequent manual inspection")
    src_group = p1.add_mutually_exclusive_group(required=True)
    src_group.add_argument("--pd-file", default=None, help="One representative prediction file")
    src_group.add_argument("--pred-dir", default=None, help="Prediction directory; must be used together with --gene")
    p1.add_argument("--gene", default=None, help="Gene used for range detection; required when --pred-dir is used")
    p1.add_argument("--tissue-positions", required=True, help="tissue_positions parquet path")
    p1.add_argument("--gt-data", required=True, help="10x h5 path")
    p1.add_argument("--out-range-json", required=True, help="Output range json path")

    p2 = sub.add_parser("evaluate", help="Evaluate with a precomputed range json")
    p2.add_argument("--algorithm", required=True, help="iSTAR, scstGCN, TESLA, XFuse, or SuperST")
    p2.add_argument("--prediction-path", required=True, help="One sample dir or a parent dir")
    p2.add_argument("--output-path", required=True, help="Output directory")
    p2.add_argument("--gt-cache-dir", default=None, help="GT cache root")
    p2.add_argument("--range-info", required=True, help="Manually finalized prediction and GT crop-range JSON")
    p2.add_argument("--tissue-positions", required=True, help="tissue_positions parquet path")
    p2.add_argument("--gt-data", required=True, help="10x h5 path")
    p2.add_argument("--no-recursive", action="store_true")
    p2.add_argument("--force", action="store_true")
    p2.add_argument("--eps", type=float, default=EPS)
    p2.add_argument("--manually-reviewed", action="store_true", help="Confirm auxiliary ranges have been inspected and finalized manually")

    return parser


def run_evaluate(args):
    import pyarrow.parquet as pq
    import scanpy as sc
    args.algorithm = algorithm_key(args.algorithm)

    prediction_path = os.path.abspath(args.prediction_path)
    range_info = load_json(os.path.abspath(args.range_info))
    if range_info.get("auxiliary_check_only", False) and not args.manually_reviewed:
        raise ValueError("Review and finalize auxiliary ranges manually, then use --manually-reviewed")
    tissue_positions = pq.ParquetFile(os.path.abspath(args.tissue_positions)).read().to_pandas()
    gt_data = sc.read_10x_h5(os.path.abspath(args.gt_data))

    pred_subdirs = list(get_algorithm_pred_subdirs(args.algorithm).values())
    sample_dirs, is_single_sample = collect_sample_dirs(
        prediction_path,
        sample_check=lambda p: is_visiumhd_sample_dir(p, pred_subdirs),
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
                gt_data=gt_data,
                tissue_positions=tissue_positions,
                range_info=range_info,
                args=args,
            )


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "detect-range":
        if args.pred_dir and not args.gene:
            raise ValueError("--gene is required when using --pred-dir")
        detect_range(args)
    elif args.command == "evaluate":
        run_evaluate(args)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()