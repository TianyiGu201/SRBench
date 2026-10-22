#!/usr/bin/env python3

import csv
import gzip
import json
import os
import pickle
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.stats import entropy, pearsonr, spearmanr

EPS = 1e-12
COLUMN_NAMES = ["gene", "gt_counts", "PCC", "SCC", "KLD", "JSD", "SSIM", "PSNR"]
VALID_ALGORITHMS = {"istar", "scstgcn", "tesla", "xfuse", "superst"}


def ensure_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)


def read_lines(path):
    with open(path, "r", encoding="utf-8") as f:
        return [line.rstrip("\n") for line in f]


def save_json(obj, path):
    ensure_dir(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_he_hw(path):
    with Image.open(path) as im:
        w, h = im.size
    return h, w


def load_image_gray(path, out_shape=None):
    with Image.open(path) as im:
        gray = im.convert("L")
        if out_shape is not None:
            gray = gray.resize((int(out_shape[1]), int(out_shape[0])), Image.Resampling.BILINEAR)
        return np.asarray(gray, dtype=np.float64)


def load_pickle_any(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    with open(path, "rb") as f:
        return pickle.load(f)


def dump_pkl_gz(obj, path, compresslevel=6):
    ensure_dir(os.path.dirname(path))
    with gzip.open(path, "wb", compresslevel=compresslevel) as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)


def find_pred_file(pred_dir, gene):
    for ext in (".pkl.gz", ".pkl", ".pickle"):
        path = os.path.join(pred_dir, f"{gene}{ext}")
        if os.path.exists(path):
            return path
    return None


def algorithm_key(name):
    key = str(name).strip().lower()
    if key not in VALID_ALGORITHMS:
        raise ValueError(f"Unsupported algorithm: {name}. Use one of: iSTAR, scstGCN, TESLA, XFuse, SuperST")
    return key


def get_algorithm_pred_subdirs(algorithm):
    algorithm = algorithm_key(algorithm)
    if algorithm in {"istar", "scstgcn"}:
        return {"main": "cnts-super"}
    if algorithm == "tesla":
        return {"main": "tesla_enhanced"}
    if algorithm == "xfuse":
        return {"main": os.path.join("my-run", "analyses", "final", "analysis-gene_maps", "section0")}
    if algorithm == "superst":
        return {"SuperST_t": "SuperST_enhanced_t", "SuperST_t_g": "SuperST_enhanced_t_g"}
    raise AssertionError("unreachable")


def make_gt_cache_map(gt_cache_root, algorithm):
    algorithm = algorithm_key(algorithm)
    if algorithm == "superst":
        return {
            "SuperST_t": os.path.join(gt_cache_root, "SuperST_t"),
            "SuperST_t_g": os.path.join(gt_cache_root, "SuperST_t_g"),
        }
    return {"main": gt_cache_root}


def make_output_dir_map(output_root, algorithm):
    algorithm = algorithm_key(algorithm)
    if algorithm == "superst":
        return {
            "SuperST_t": os.path.join(output_root, "SuperST_t"),
            "SuperST_t_g": os.path.join(output_root, "SuperST_t_g"),
        }
    return {"main": output_root}


def resolve_default_gt_cache_dir(prediction_path, is_single_sample, suffix="gt_cache"):
    base = os.path.abspath(prediction_path)
    parent = os.path.dirname(base)
    return os.path.join(parent, suffix)


def _has_gene_names(path):
    return os.path.isfile(os.path.join(path, "gene-names.txt"))


def is_xenium_sample_dir(path, pred_subdirs):
    if not os.path.isdir(path) or not _has_gene_names(path):
        return False
    for subdir in pred_subdirs:
        if os.path.isdir(os.path.join(path, subdir)):
            return True
    return False


def is_visiumhd_sample_dir(path, pred_subdirs):
    return is_xenium_sample_dir(path, pred_subdirs)


def collect_sample_dirs(prediction_path, sample_check, recursive=True):
    prediction_path = os.path.abspath(prediction_path)
    if not os.path.isdir(prediction_path):
        raise FileNotFoundError(f"Directory not found: {prediction_path}")

    if sample_check(prediction_path):
        return [prediction_path], True

    sample_dirs = []
    seen = set()

    if recursive:
        walker = os.walk(prediction_path)
    else:
        try:
            child_dirs = next(os.walk(prediction_path))[1]
        except StopIteration:
            child_dirs = []
        walker = [(prediction_path, child_dirs, [])]

    for root, _, _ in walker:
        if root in seen:
            continue
        if sample_check(root):
            sample_dirs.append(root)
            seen.add(root)

    sample_dirs.sort()
    return sample_dirs, False


def build_unique_names(sample_dirs):
    def parts_of(path):
        p = Path(path).resolve()
        out = [p.name]
        for parent in p.parents:
            if not parent.name:
                break
            out.append(parent.name)
        return out

    parts_map = {d: parts_of(d) for d in sample_dirs}
    depth_map = {d: 1 for d in sample_dirs}

    while True:
        buckets = {}
        for d in sample_dirs:
            parts = parts_map[d]
            depth = min(depth_map[d], len(parts))
            name = "_".join(reversed(parts[:depth]))
            buckets.setdefault(name, []).append(d)
        conflicts = [v for v in buckets.values() if len(v) > 1]
        if not conflicts:
            break
        progressed = False
        for group in conflicts:
            for d in group:
                if depth_map[d] < len(parts_map[d]):
                    depth_map[d] += 1
                    progressed = True
        if not progressed:
            break

    out = {}
    used = set()
    for d in sample_dirs:
        parts = parts_map[d]
        depth = min(depth_map[d], len(parts))
        name = "_".join(reversed(parts[:depth]))
        base = name
        k = 2
        while name in used:
            name = f"{base}_{k}"
            k += 1
        out[d] = name
        used.add(name)
    return out


def should_filter_gene(gene):
    if gene is None:
        return True
    gene = str(gene).strip()
    return gene == ""


def safe_minmax(a):
    a = np.asarray(a, dtype=np.float64)
    if a.size == 0:
        return a
    amin = float(np.nanmin(a))
    amax = float(np.nanmax(a))
    if not np.isfinite(amin) or not np.isfinite(amax) or amax - amin <= EPS:
        return np.zeros_like(a, dtype=np.float64)
    return (a - amin) / (amax - amin)


def safe_corr(y_true, y_pred):
    yt = np.asarray(y_true, dtype=np.float64).ravel()
    yp = np.asarray(y_pred, dtype=np.float64).ravel()
    mask = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[mask]
    yp = yp[mask]
    if yt.size < 2 or np.std(yt) <= EPS or np.std(yp) <= EPS:
        return np.nan, np.nan
    return float(pearsonr(yt, yp)[0]), float(spearmanr(yt, yp)[0])


def to_prob(x, eps=EPS):
    x = np.asarray(x, dtype=np.float64).ravel()
    if x.size == 0:
        return np.array([1.0], dtype=np.float64)
    x = np.maximum(x, eps)
    s = float(x.sum())
    if s <= 0:
        return np.full_like(x, 1.0 / x.size, dtype=np.float64)
    return x / s


def kl_divergence(p, q, eps=EPS):
    return float(entropy(to_prob(p, eps), to_prob(q, eps), base=2))


def js_divergence(p, q, eps=EPS):
    p = to_prob(p, eps)
    q = to_prob(q, eps)
    m = 0.5 * (p + q)
    return float(0.5 * entropy(p, m, base=2) + 0.5 * entropy(q, m, base=2))


def safe_psnr(a, b):
    a = safe_minmax(a)
    b = safe_minmax(b)
    mse = float(np.mean((a - b) ** 2))
    if mse <= EPS:
        return np.inf
    return float(10.0 * np.log10(1.0 / mse))


def safe_ssim(a, b):
    a = safe_minmax(a)
    b = safe_minmax(b)
    c1 = 0.01 ** 2
    c2 = 0.03 ** 2
    mu_a = float(np.mean(a))
    mu_b = float(np.mean(b))
    var_a = float(np.var(a))
    var_b = float(np.var(b))
    cov = float(np.mean((a - mu_a) * (b - mu_b)))
    num = (2 * mu_a * mu_b + c1) * (2 * cov + c2)
    den = (mu_a ** 2 + mu_b ** 2 + c1) * (var_a + var_b + c2)
    if abs(den) <= EPS:
        return np.nan
    return float(num / den)


def evaluate_pair(y_true, y_pred, eps=EPS):
    yt = np.asarray(y_true, dtype=np.float64)
    yp = np.asarray(y_pred, dtype=np.float64)
    yt1 = yt.ravel()
    yp1 = yp.ravel()
    pcc, scc = safe_corr(yt1, yp1)
    return [
        float(np.sum(yt1)),
        pcc,
        scc,
        kl_divergence(yt1, yp1, eps=eps),
        js_divergence(yt1, yp1, eps=eps),
        safe_ssim(yt, yp),
        safe_psnr(yt, yp),
    ]


def apply_xfuse_transform(arr, do_transpose=False, do_fliplr=False):
    arr = np.asarray(arr)
    if do_transpose:
        arr = arr.T
    if do_fliplr:
        arr = np.fliplr(arr)
    return np.asarray(arr, dtype=np.float64)


def crop_array_by_range(arr, row_range=None, col_range=None):
    arr = np.asarray(arr)
    h, w = arr.shape[:2]

    if row_range is None:
        r0, r1 = 0, h
    else:
        r0, r1 = int(row_range[0]), int(row_range[1])

    if col_range is None:
        c0, c1 = 0, w
    else:
        c0, c1 = int(col_range[0]), int(col_range[1])

    if not (0 <= r0 <= r1 <= h):
        raise ValueError(f"Invalid row range: {(r0, r1)} for shape {arr.shape}")
    if not (0 <= c0 <= c1 <= w):
        raise ValueError(f"Invalid col range: {(c0, c1)} for shape {arr.shape}")

    return arr[r0:r1, c0:c1]


def istar_scale_pred(arr, sample_dir):
    """
    Minimal iSTAR preprocessing:
    only crop prediction to the top-left area whose shape is
    (he_scaled_h // 16, he_scaled_w // 16).
    """
    out = np.asarray(arr, dtype=np.float64)

    he_scaled = os.path.join(sample_dir, "he-scaled.jpg")
    if not os.path.exists(he_scaled):
        return out

    he_h, he_w = load_he_hw(he_scaled)
    target_h = max(int(he_h // 16), 1)
    target_w = max(int(he_w // 16), 1)

    crop_h = min(out.shape[0], target_h)
    crop_w = min(out.shape[1], target_w)

    return out[:crop_h, :crop_w]


def load_transcripts_tsv(path, feature_col="feature_name"):
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def build_gt_gene_coords(rows, feature_col="feature_name", y_col_idx=2, x_col_idx=1):
    if not rows:
        return {}

    headers = list(rows[0].keys())
    y_col = headers[int(y_col_idx)]
    x_col = headers[int(x_col_idx)]

    out = {}
    for row in rows:
        gene = row[feature_col]
        if gene not in out:
            out[gene] = []
        out[gene].append([float(row[y_col]), float(row[x_col])])

    for gene in out:
        out[gene] = np.asarray(out[gene], dtype=np.float64)
    return out


def build_gt_hist(coords, pd_h, pd_w, he_h, he_w):
    hist = np.zeros((int(pd_h), int(pd_w)), dtype=np.float64)
    if coords is None or len(coords) == 0:
        return hist

    coords = np.asarray(coords, dtype=np.float64)
    y = np.clip(np.floor(coords[:, 0] / max(he_h / pd_h, EPS)).astype(int), 0, pd_h - 1)
    x = np.clip(np.floor(coords[:, 1] / max(he_w / pd_w, EPS)).astype(int), 0, pd_w - 1)
    np.add.at(hist, (y, x), 1.0)
    return hist


def build_visiumhd_gt_cache_name(gt_row_range, gt_col_range):
    return f"r{int(gt_row_range[0])}_{int(gt_row_range[1])}__c{int(gt_col_range[0])}_{int(gt_col_range[1])}"


def _find_col(df, candidates):
    lower = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    raise KeyError(f"Missing columns: {candidates}")


def align_positions_to_adata(gt_data, tissue_positions):
    """Align positions to expression rows by unique barcode; never truncate silently."""
    positions = tissue_positions.copy()
    if "barcode" in positions.columns:
        positions["barcode"] = positions["barcode"].astype(str)
        positions = positions.set_index("barcode", verify_integrity=True)
    else:
        if positions.index.name != "barcode":
            raise ValueError("Positions must contain a barcode column or barcode-named index")
        positions.index = positions.index.astype(str)
        if not positions.index.is_unique:
            raise ValueError("Duplicate barcodes in tissue positions")
    obs_names = [str(v) for v in gt_data.obs_names]
    if len(obs_names) != len(set(obs_names)):
        raise ValueError("Duplicate barcodes in expression matrix")
    missing = set(obs_names).difference(positions.index)
    if missing:
        raise ValueError(f"Missing coordinates for {len(missing)} expression barcodes")
    return positions.loc[obs_names]


def build_full_visiumhd_gene_grid(gene, gt_data, tissue_positions):
    if gene not in gt_data.var_names:
        raise KeyError(gene)

    tissue_positions = align_positions_to_adata(gt_data, tissue_positions)

    row_col = _find_col(tissue_positions, ["array_row", "row", "pxl_row_in_fullres"])
    col_col = _find_col(tissue_positions, ["array_col", "col", "pxl_col_in_fullres"])

    expr = gt_data[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().ravel()
    else:
        expr = np.asarray(expr).ravel()

    rows = tissue_positions[row_col].to_numpy().astype(int)
    cols = tissue_positions[col_col].to_numpy().astype(int)

    if not (len(expr) == len(rows) == len(cols)):
        raise ValueError("Expression and coordinate row counts differ after barcode alignment")

    row_min = int(rows.min())
    row_max = int(rows.max())
    col_min = int(cols.min())
    col_max = int(cols.max())

    grid = np.zeros((row_max - row_min + 1, col_max - col_min + 1), dtype=np.float64)
    rr = rows - row_min
    cc = cols - col_min
    np.add.at(grid, (rr, cc), expr)

    meta = {
        "row_min": row_min,
        "row_max": row_max,
        "col_min": col_min,
        "col_max": col_max,
    }
    return grid, meta


def gene_grid_from_anndata(gene, gt_data, tissue_positions, gt_row_range, gt_col_range, flip_rows=False, flip_cols=False):
    if gene not in gt_data.var_names:
        raise KeyError(gene)

    tissue_positions = align_positions_to_adata(gt_data, tissue_positions)

    row_col = _find_col(tissue_positions, ["array_row", "row", "pxl_row_in_fullres"])
    col_col = _find_col(tissue_positions, ["array_col", "col", "pxl_col_in_fullres"])

    expr = gt_data[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().ravel()
    else:
        expr = np.asarray(expr).ravel()

    rows = tissue_positions[row_col].to_numpy().astype(int)
    cols = tissue_positions[col_col].to_numpy().astype(int)

    if not (len(expr) == len(rows) == len(cols)):
        raise ValueError("Expression and coordinate row counts differ after barcode alignment")

    gt_row_range = tuple(int(v) for v in gt_row_range)
    gt_col_range = tuple(int(v) for v in gt_col_range)

    n_rows = int(gt_row_range[1] - gt_row_range[0])
    n_cols = int(gt_col_range[1] - gt_col_range[0])

    grid = np.zeros((n_rows, n_cols), dtype=np.float64)

    rr = rows - int(gt_row_range[0])
    cc = cols - int(gt_col_range[0])

    mask = (rr >= 0) & (rr < n_rows) & (cc >= 0) & (cc < n_cols)
    np.add.at(grid, (rr[mask], cc[mask]), expr[mask])

    if flip_rows:
        grid = np.flipud(grid)
    if flip_cols:
        grid = np.fliplr(grid)

    return grid


def _bbox_from_mask(mask):
    coords = np.argwhere(mask)
    if coords.size == 0:
        return (0, int(mask.shape[0]), 0, int(mask.shape[1]))
    r0, c0 = coords.min(axis=0)
    r1, c1 = coords.max(axis=0) + 1
    return (int(r0), int(r1), int(c0), int(c1))


def _fit_range_with_size(center, size, low, high):
    size = int(size)
    low = int(low)
    high = int(high)

    if size <= 0:
        raise ValueError(f"size must be positive, got {size}")

    full = high - low
    if size > full:
        return (low, high)

    start = int(round(center - size / 2.0))
    end = start + size

    if start < low:
        start = low
        end = start + size
    if end > high:
        end = high
        start = end - size

    return (int(start), int(end))


def infer_visiumhd_range_info(pd_array, gt_data, tissue_positions, gene=None):
    pd_array = np.asarray(pd_array, dtype=np.float64)
    if pd_array.ndim != 2:
        raise ValueError(f"pd_array must be 2D, got shape={pd_array.shape}")

    if gene is None:
        gene = str(gt_data.var_names[0])

    gt_grid_full, meta = build_full_visiumhd_gene_grid(gene, gt_data, tissue_positions)

    pd_mask = np.abs(pd_array) > 0
    gt_mask = np.abs(gt_grid_full) > 0

    pd_bbox = _bbox_from_mask(pd_mask)
    gt_bbox_local = _bbox_from_mask(gt_mask)

    pd_h = int(pd_bbox[1] - pd_bbox[0])
    pd_w = int(pd_bbox[3] - pd_bbox[2])

    gt_center_r = (gt_bbox_local[0] + gt_bbox_local[1]) / 2.0
    gt_center_c = (gt_bbox_local[2] + gt_bbox_local[3]) / 2.0

    gt_row_local = _fit_range_with_size(
        center=gt_center_r,
        size=pd_h,
        low=0,
        high=gt_grid_full.shape[0],
    )
    gt_col_local = _fit_range_with_size(
        center=gt_center_c,
        size=pd_w,
        low=0,
        high=gt_grid_full.shape[1],
    )

    gt_row_range = [
        int(meta["row_min"] + gt_row_local[0]),
        int(meta["row_min"] + gt_row_local[1]),
    ]
    gt_col_range = [
        int(meta["col_min"] + gt_col_local[0]),
        int(meta["col_min"] + gt_col_local[1]),
    ]

    info = {
        "gene_for_range_detection": str(gene),
        "pd_shape": [int(pd_array.shape[0]), int(pd_array.shape[1])],
        "pd_row_range": [int(pd_bbox[0]), int(pd_bbox[1])],
        "pd_col_range": [int(pd_bbox[2]), int(pd_bbox[3])],
        "gt_row_range": gt_row_range,
        "gt_col_range": gt_col_range,
        "flip_rows": False,
        "flip_cols": False,
        "gt_full_meta": {
            "row_min": int(meta["row_min"]),
            "row_max": int(meta["row_max"]),
            "col_min": int(meta["col_min"]),
            "col_max": int(meta["col_max"]),
        },
    }
    return info


def load_pred_array_for_range_detection(pd_path):
    return np.asarray(load_pickle_any(pd_path), dtype=np.float64)