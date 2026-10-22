#!/usr/bin/env python3
"""Rank candidate datasets by similarity to query datasets using shared metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


SUPPORTED_COLUMNS = [
    "he_image_feature_json",
    "tissue_type",
    "spatial_resolution_um_per_pixel",
    "data_sparsity",
    "spatial_autocorrelation",
    "mean_expression_per_100um2",
    "cv_expression_density",
    "svg_number",
]

CATEGORY_MAP = {
    "he_image_feature_json": "morphological",
    "tissue_type": "morphological",
    "spatial_resolution_um_per_pixel": "technical",
    "data_sparsity": "technical",
    "spatial_autocorrelation": "technical",
    "mean_expression_per_100um2": "gene_expression",
    "cv_expression_density": "gene_expression",
    "svg_number": "gene_expression",
}

ALIASES = {
    "he_image_feature": "he_image_feature_json",
    "uni_embedding_json": "he_image_feature_json",
    "spatial_resolution_um": "spatial_resolution_um_per_pixel",
    "moranI": "spatial_autocorrelation",
    "spatial_autocorrelation_top50hvg_moranI": "spatial_autocorrelation",
    "mean_gene_expression": "mean_expression_per_100um2",
    "cv_gene_expression": "cv_expression_density",
    "svg_count": "svg_number",
}


def canonicalize_column_name(name: str) -> str:
    return ALIASES.get(name.strip(), name.strip())


def parse_user_columns(columns_str: Optional[str]) -> Optional[List[str]]:
    if not columns_str:
        return None
    columns = [canonicalize_column_name(x) for x in columns_str.split(",") if x.strip()]
    invalid = [c for c in columns if c not in SUPPORTED_COLUMNS]
    if invalid:
        raise ValueError(f"Unsupported columns: {invalid}")
    return columns


def parse_user_weights(weights_str: Optional[str], columns: Optional[Sequence[str]]) -> Optional[Dict[str, float]]:
    if not weights_str:
        return None
    if not columns:
        raise ValueError("--weights requires --columns")
    weights = {c: 1.0 for c in columns}
    for item in weights_str.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"Invalid weight specification: {item}")
        key, value = item.split("=", 1)
        key = canonicalize_column_name(key)
        if key not in columns:
            raise ValueError(f"Weight provided for a column not present in --columns: {key}")
        weights[key] = float(value.strip())
    return weights


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        return df
    df = df.copy()
    df.columns = [canonicalize_column_name(c) for c in df.columns]
    return df


def parse_embedding(value) -> Optional[np.ndarray]:
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        arr = np.asarray(json.loads(value), dtype=np.float32)
    except Exception:
        return None
    if arr.ndim != 1 or arr.size == 0:
        return None
    return arr


def has_value(row: pd.Series, column: str) -> bool:
    if column not in row.index:
        return False
    value = row[column]
    if isinstance(value, str):
        return bool(value.strip())
    return not pd.isna(value)


def safe_float(value) -> Optional[float]:
    try:
        if pd.isna(value):
            return None
        return float(value)
    except Exception:
        return None


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> Optional[float]:
    if a.shape != b.shape:
        return None
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return None
    return float(np.dot(a, b) / denom)


def bounded_relative_similarity(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    denom = max(a, b, 1e-12)
    score = 1.0 - abs(a - b) / denom
    return float(max(0.0, min(1.0, score)))


def tissue_similarity(a: str, b: str) -> Optional[float]:
    a = (a or "").strip().lower()
    b = (b or "").strip().lower()
    if not a or not b:
        return None
    return 1.0 if a == b else 0.0


def is_comparable(query_row: pd.Series, candidate_row: pd.Series, column: str) -> bool:
    if column == "he_image_feature_json":
        return parse_embedding(query_row.get(column, "")) is not None and parse_embedding(candidate_row.get(column, "")) is not None
    return has_value(query_row, column) and has_value(candidate_row, column)


def compute_metric_similarity(query_row: pd.Series, candidate_row: pd.Series, column: str) -> Optional[float]:
    if column == "he_image_feature_json":
        q = parse_embedding(query_row.get(column, ""))
        t = parse_embedding(candidate_row.get(column, ""))
        return cosine_similarity(q, t) if q is not None and t is not None else None
    if column == "tissue_type":
        return tissue_similarity(str(query_row.get(column, "")), str(candidate_row.get(column, "")))
    if column == "spatial_resolution_um_per_pixel":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    if column == "data_sparsity":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    if column == "spatial_autocorrelation":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    if column == "mean_expression_per_100um2":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    if column == "cv_expression_density":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    if column == "svg_number":
        return bounded_relative_similarity(safe_float(query_row.get(column)), safe_float(candidate_row.get(column)))
    raise ValueError(f"Unsupported column: {column}")


def weighted_average(pairs: Sequence[Tuple[Optional[float], float]]) -> Optional[float]:
    valid = [(score, weight) for score, weight in pairs if score is not None and weight > 0]
    if not valid:
        return None
    total_weight = sum(weight for _, weight in valid)
    return float(sum(score * weight for score, weight in valid) / total_weight)


def get_default_common_columns(query_row: pd.Series, candidate_row: pd.Series) -> List[str]:
    return [c for c in SUPPORTED_COLUMNS if is_comparable(query_row, candidate_row, c)]


def get_active_scheme(
    query_row: pd.Series,
    candidate_row: pd.Series,
    user_columns: Optional[Sequence[str]],
    user_weights: Optional[Dict[str, float]],
) -> Tuple[List[str], Dict[str, float], str]:
    if user_columns and all(is_comparable(query_row, candidate_row, c) for c in user_columns):
        return list(user_columns), (user_weights or {c: 1.0 for c in user_columns}), "user_defined"
    active_columns = get_default_common_columns(query_row, candidate_row)
    return active_columns, {c: 1.0 for c in active_columns}, "default_common"


def compare_one(
    query_row: pd.Series,
    candidate_row: pd.Series,
    active_columns: Sequence[str],
    active_weights: Dict[str, float],
) -> Dict[str, object]:
    metric_scores = {c: compute_metric_similarity(query_row, candidate_row, c) for c in active_columns}

    category_scores = {}
    for category in ["morphological", "technical", "gene_expression"]:
        cols = [c for c in active_columns if CATEGORY_MAP[c] == category]
        category_scores[category] = weighted_average([(metric_scores[c], active_weights[c]) for c in cols])

    overall_similarity = weighted_average([(metric_scores[c], active_weights[c]) for c in active_columns])

    row = {
        "overall_similarity": overall_similarity,
        "morphological_similarity": category_scores["morphological"],
        "technical_similarity": category_scores["technical"],
        "gene_expression_similarity": category_scores["gene_expression"],
    }
    for col in SUPPORTED_COLUMNS:
        row[f"{col}_similarity"] = metric_scores.get(col, np.nan)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank candidate datasets from one candidate CSV against one query CSV")
    parser.add_argument("--query-csv", required=True, help="CSV containing one or more query rows")
    parser.add_argument("--candidate-csv", required=True, help="CSV containing one or more candidate rows")
    parser.add_argument("--columns", default=None, help="Optional comma-separated columns to force")
    parser.add_argument("--weights", default=None, help="Optional weights for --columns")
    parser.add_argument("--topk", type=int, default=3, help="Top-k candidates to keep for each query row. Default: 3")
    parser.add_argument("--output-csv", default=None, help="Optional output CSV path")
    args = parser.parse_args()

    user_columns = parse_user_columns(args.columns)
    user_weights = parse_user_weights(args.weights, user_columns)

    query_df = load_csv(args.query_csv)
    candidate_df = load_csv(args.candidate_csv)
    if query_df.empty:
        raise ValueError("Query CSV is empty")
    if candidate_df.empty:
        raise ValueError("Candidate CSV is empty")

    results = []
    for _, query_row in query_df.iterrows():
        query_name = query_row.get("dataset_name", "")
        pair_rows = []
        for _, candidate_row in candidate_df.iterrows():
            active_columns, active_weights, scheme_used = get_active_scheme(
                query_row=query_row,
                candidate_row=candidate_row,
                user_columns=user_columns,
                user_weights=user_weights,
            )
            base = {
                "query_dataset_name": query_name,
                "candidate_dataset_name": candidate_row.get("dataset_name", ""),
                "columns_used": json.dumps(active_columns),
                "weights_used_json": json.dumps(active_weights),
                "scheme_used": scheme_used,
            }
            if active_columns:
                sim = compare_one(query_row, candidate_row, active_columns, active_weights)
            else:
                sim = {
                    "overall_similarity": np.nan,
                    "morphological_similarity": np.nan,
                    "technical_similarity": np.nan,
                    "gene_expression_similarity": np.nan,
                }
                for col in SUPPORTED_COLUMNS:
                    sim[f"{col}_similarity"] = np.nan
            pair_rows.append({**base, **sim})

        pair_df = pd.DataFrame(pair_rows)
        pair_df = pair_df.sort_values("overall_similarity", ascending=False, na_position="last").reset_index(drop=True)
        pair_df.insert(0, "rank_within_query", np.arange(1, len(pair_df) + 1))
        results.append(pair_df.head(args.topk))

    result_df = pd.concat(results, ignore_index=True)
    if args.output_csv:
        output_path = Path(args.output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_path, index=False)

    print(result_df.to_string(index=False))


if __name__ == "__main__":
    main()
