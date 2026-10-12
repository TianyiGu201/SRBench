"""
FastAPI backend for Spatial Transcriptomics Visualization Platform
Uses .npy format exclusively for gene expression data
"""
import os
import json
import re
import math
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Optional
from io import BytesIO

from fastapi import FastAPI, HTTPException, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dotenv import load_dotenv

load_dotenv(override=True)

app = FastAPI(title="SpatialVis API", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost",
        "http://localhost:3000",
        "http://127.0.0.1",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
    allow_headers=["Accept", "Authorization", "Content-Type", "Origin", "User-Agent"],
)

DATA_ROOT = Path(os.getenv("DATA_ROOT", "./data_root"))


# ============================================================================
# Core IO
# ============================================================================

def load_array(file_path: Path) -> np.ndarray:
    if not file_path.exists():
        raise HTTPException(status_code=404, detail=f"Data file not found: {file_path.name}")
    if file_path.suffix != ".npy":
        raise HTTPException(status_code=400, detail=f"Expected .npy file, got: {file_path.name}")
    try:
        data = np.load(file_path, allow_pickle=False)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading {file_path.name}: {e}")
    data = np.asarray(data, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def list_gene_array_files(dir_path: Path) -> List[Path]:
    if not dir_path.exists():
        return []
    return [p for p in dir_path.iterdir() if p.is_file() and p.suffix == ".npy"]


def build_gene_index_from_dir(dir_path: Path) -> Dict:
    return {
        "genes": [{"name": f.stem, "file": f.name} for f in list_gene_array_files(dir_path)]
    }


def load_gene_index(data_dir: Path) -> dict:
    gene_index_path = data_dir / "gene_index.json"
    if gene_index_path.exists():
        with open(gene_index_path, "r", encoding="utf-8") as f:
            idx = json.load(f)
        genes = [g for g in idx.get("genes", []) if g["file"].endswith(".npy")]
        return {"genes": genes}
    return build_gene_index_from_dir(data_dir)


def get_gene_spatial_data(data_dir: Path, gene: str, gene_index: dict) -> Dict:
    gene_info = next((g for g in gene_index["genes"] if g["name"] == gene), None)
    if not gene_info:
        raise HTTPException(status_code=404, detail=f"Gene '{gene}' not found")
    data = load_array(data_dir / gene_info["file"])
    n_rows, n_cols = data.shape
    points = [
        {"x": float(j), "y": float(i), "expression": float(data[i, j])}
        for i in range(n_rows)
        for j in range(n_cols)
        if data[i, j] > 0
    ]
    return {"gene": gene, "points": points, "shape": [n_cols, n_rows]}


def render_array_image(data: np.ndarray, cmap: str, vmin: Optional[float], vmax: Optional[float]) -> StreamingResponse:
    fig, ax = plt.subplots(figsize=(data.shape[1] / 100, data.shape[0] / 100), dpi=100)
    ax.imshow(data, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.axis("off")
    plt.tight_layout(pad=0)
    buf = BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    buf.seek(0)
    return StreamingResponse(buf, media_type="image/png")


# ============================================================================
# Helpers
# ============================================================================

def get_dataset_path(dataset_name: str) -> Path:
    path = DATA_ROOT / dataset_name
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Dataset '{dataset_name}' not found")
    return path


def list_datasets() -> List[str]:
    if not DATA_ROOT.exists():
        return []
    return [d.name for d in DATA_ROOT.iterdir() if d.is_dir()]


def list_csv_files(dir_path: Path) -> List[Path]:
    if not dir_path.exists():
        return []
    return [p for p in dir_path.iterdir() if p.is_file() and p.suffix.lower() == ".csv"]


def parse_FOV_region(region_str: str) -> Optional[Dict]:
    m = re.match(r"minx([\d.]+)_maxx([\d.]+)_miny([\d.]+)_maxy([\d.]+)", region_str)
    if m:
        return {
            "xStart": float(m.group(1)),
            "xEnd": float(m.group(2)),
            "yStart": float(m.group(3)),
            "yEnd": float(m.group(4)),
        }
    return None


def normalize_region_id(region_id: str) -> str:
    return str(region_id).strip()


def extract_numeric_value(filename: str) -> float:
    m = re.search(r"(\d+\.?\d*)", filename)
    return float(m.group(1)) if m else 0.0


def read_genes_from_resolution_csv(csv_path: Path) -> List[str]:
    df = pd.read_csv(csv_path)
    if df.shape[1] == 0:
        return []
    gene_col = df.columns[0]
    genes = []
    for v in df[gene_col].tolist():
        s = v.strip() if isinstance(v, str) else str(v)
        if s:
            genes.append(s)
    return sorted(set(genes))


# Directory helpers
def metric_vis_root(dataset_name: str) -> Path:
    return get_dataset_path(dataset_name) / "metric_visualization"

def metric_vis_reference_dir(dataset_name: str) -> Path:
    return metric_vis_root(dataset_name) / "reference"

def metric_vis_prediction_dir(dataset_name: str, algorithm: str = "", resolution: str = "") -> Path:
    return metric_vis_root(dataset_name) / "prediction"

def diff_res_root(dataset_name: str, algorithm: str) -> Path:
    return DATA_ROOT / dataset_name / "prediction" / algorithm / "raw_resolution"

def diff_res_rep_arr_dir(dataset_name: str, algorithm: str, rep: str, arrangement: str) -> Path:
    return diff_res_root(dataset_name, algorithm) / rep / arrangement

# metric
def find_metric_value_entry(metric_map: Dict, metric_value: str):
    """
    Safely find the exact or uniquely equivalent key for metric_value.

    Matching order:
    1. exact string match
    2. stripped string match
    3. numeric-equivalent unique match, e.g. '40' == '40.0'

    If multiple numeric-equivalent keys exist, return None to avoid ambiguity.
    """
    # 1) exact match
    if metric_value in metric_map:
        return metric_value, metric_map[metric_value]

    # 2) stripped exact match
    stripped_target = metric_value.strip()
    if stripped_target in metric_map:
        return stripped_target, metric_map[stripped_target]

    for k, v in metric_map.items():
        if str(k).strip() == stripped_target:
            return k, v

    # 3) numeric-equivalent unique match
    try:
        target_num = float(stripped_target)
    except Exception:
        return None, None

    matched = []
    for k, v in metric_map.items():
        try:
            if float(str(k).strip()) == target_num:
                matched.append((k, v))
        except Exception:
            continue

    if len(matched) == 1:
        return matched[0]

    return None, None


# ============================================================================
# 4.1 Dataset Overview
# ============================================================================

@app.get("/api/datasets")
async def get_datasets():
    datasets = []
    for dataset_name in list_datasets():
        try:
            dataset_path = get_dataset_path(dataset_name)
            metadata_path = dataset_path / "dataset_metaINFO.json"
            if metadata_path.exists():
                with open(metadata_path, "r", encoding="utf-8") as f:
                    metadata = json.load(f)
                datasets.append({
                    "id": dataset_name,
                    "name": metadata.get("dataset_name", dataset_name),
                    "species": metadata.get("species", "Unknown"),
                    "technology": metadata.get("sequencing_technology", "Unknown"),
                    "tissue": metadata.get("tissue_type", "Unknown"),
                    "description": f"{metadata.get('sequencing_technology', 'Unknown')} data from {metadata.get('tissue_type', 'Unknown')} tissue",
                    "heImageUrl": f"/api/datasets/{dataset_name}/he-thumbnail",
                })
        except Exception as e:
            print(f"Error loading dataset {dataset_name}: {e}")
    tech_order = {
        "Xenium": 0,
        "Visium HD": 1,
        "CosMx SMI": 2,
    }

    datasets.sort(key=lambda d: (
        tech_order.get(d["technology"], 999),  # technology priority
        d["name"]                              # secondary alphabetical
    ))

    return datasets


@app.get("/api/datasets/{dataset_name}/metadata")
async def get_dataset_metadata(dataset_name: str):
    metadata_path = get_dataset_path(dataset_name) / "dataset_metaINFO.json"
    if not metadata_path.exists():
        raise HTTPException(status_code=404, detail="Metadata not found")
    with open(metadata_path, "r", encoding="utf-8") as f:
        return json.load(f)


@app.get("/api/datasets/{dataset_name}/he-thumbnail")
async def get_he_thumbnail(dataset_name: str):
    p = get_dataset_path(dataset_name) / "he-raw.jpg"
    if not p.exists():
        raise HTTPException(status_code=404, detail="Thumbnail not found")
    return FileResponse(p, media_type="image/jpeg")


@app.get("/api/datasets/{dataset_name}/he-image")
async def get_he_image(dataset_name: str):
    p = get_dataset_path(dataset_name) / "he-raw.tiff"
    if not p.exists():
        raise HTTPException(status_code=404, detail="HE image not found")
    return FileResponse(p, media_type="image/tiff", filename=f"{dataset_name}_he.tiff")


# ============================================================================
# 4.2 Metric Visualization
# ============================================================================

@app.get("/api/metric-vis/{dataset_name}/options")
async def get_metric_options(dataset_name: str):
    p = get_dataset_path(dataset_name) / "metric_vis.json"
    if not p.exists():
        raise HTTPException(status_code=404, detail="metric_vis.json not found")
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


@app.get("/api/metric-vis/{dataset_name}/{metric_name}/{metric_value}/genes")
async def get_similar_genes(dataset_name: str, metric_name: str, metric_value: str):
    p = get_dataset_path(dataset_name) / "metric_vis.json"
    if not p.exists():
        raise HTTPException(status_code=404, detail="metric_vis.json not found")

    with open(p, "r", encoding="utf-8") as f:
        metric_vis = json.load(f)

    if metric_name not in metric_vis:
        raise HTTPException(status_code=404, detail=f"Metric '{metric_name}' not found")

    metric_map = metric_vis[metric_name]
    matched_key, genes = find_metric_value_entry(metric_map, metric_value)

    if matched_key is None:
        raise HTTPException(
            status_code=404,
            detail=f"Value '{metric_value}' not found for metric '{metric_name}'"
        )

    return {
        "metric": metric_name,
        "value": matched_key,
        "genes": genes,
    }


@app.get("/api/spatial-data/{dataset_name}/compute-global-expression-domain/{gene}")
async def compute_global_expression_domain(
    dataset_name: str, gene: str,
    algorithm: str = "iSTAR", resolution: str = "spotDistance_100",
):
    ref_dir = metric_vis_reference_dir(dataset_name)
    pred_dir = metric_vis_prediction_dir(dataset_name)
    if not ref_dir.exists():
        raise HTTPException(status_code=404, detail="reference directory not found")
    if not pred_dir.exists():
        raise HTTPException(status_code=404, detail="prediction directory not found")

    ref_pts = get_gene_spatial_data(ref_dir, gene, load_gene_index(ref_dir))["points"]
    pred_pts = get_gene_spatial_data(pred_dir, gene, load_gene_index(pred_dir))["points"]
    all_expr = [p["expression"] for p in ref_pts + pred_pts]
    if not all_expr:
        return 0, 1
    return min(all_expr), max(all_expr)


@app.get("/api/spatial-data/{dataset_name}/reference/{gene}")
async def get_reference_spatial_data(dataset_name: str, gene: str):
    ref_dir = metric_vis_reference_dir(dataset_name)
    if not ref_dir.exists():
        raise HTTPException(status_code=404, detail="reference directory not found")
    return get_gene_spatial_data(ref_dir, gene, load_gene_index(ref_dir))


@app.get("/api/spatial-data/{dataset_name}/prediction/{gene}")
async def get_prediction_spatial_data(
    dataset_name: str, gene: str,
    algorithm: str = "iSTAR", resolution: str = "spotDistance_100",
):
    pred_dir = metric_vis_prediction_dir(dataset_name)
    if not pred_dir.exists():
        raise HTTPException(status_code=404, detail="prediction directory not found")
    return get_gene_spatial_data(pred_dir, gene, load_gene_index(pred_dir))


@app.get("/api/spatial-image/{dataset_name}/reference/{gene}")
async def get_reference_spatial_image(
    dataset_name: str, gene: str,
    vmin: Optional[float] = None, vmax: Optional[float] = None,
):
    ref_dir = metric_vis_reference_dir(dataset_name)
    if not ref_dir.exists():
        raise HTTPException(status_code=404, detail="reference directory not found")
    gene_index = load_gene_index(ref_dir)
    gene_info = next((g for g in gene_index["genes"] if g["name"] == gene), None)
    if not gene_info:
        raise HTTPException(status_code=404, detail=f"Gene '{gene}' not found")
    data = load_array(ref_dir / gene_info["file"])
    return render_array_image(data, "RdBu_r", vmin, vmax)


@app.get("/api/spatial-image/{dataset_name}/prediction/{gene}")
async def get_prediction_spatial_image(
    dataset_name: str, gene: str,
    algorithm: str = "iSTAR", resolution: str = "spotDistance_100",
    vmin: Optional[float] = None, vmax: Optional[float] = None,
):
    pred_dir = metric_vis_prediction_dir(dataset_name)
    if not pred_dir.exists():
        raise HTTPException(status_code=404, detail="prediction directory not found")
    gene_index = load_gene_index(pred_dir)
    gene_info = next((g for g in gene_index["genes"] if g["name"] == gene), None)
    if not gene_info:
        raise HTTPException(status_code=404, detail=f"Gene '{gene}' not found")
    data = load_array(pred_dir / gene_info["file"])
    return render_array_image(data, "RdBu_r", vmin, vmax)


# ============================================================================
# 4.3 Out-of-sample Prediction
# ============================================================================

@app.get("/api/out-of-sample/{dataset_name}/{algorithm}/regions")
async def get_out_of_sample_regions(dataset_name: str, algorithm: str):
    oos_dir = get_dataset_path(dataset_name) / "prediction" / algorithm / "out_of_sample_prediction"
    if not oos_dir.exists():
        raise HTTPException(status_code=404, detail="Out-of-sample prediction directory not found")

    region_to_parnorms: Dict[str, List[float]] = defaultdict(list)
    region_to_coords: Dict[str, Dict] = {}

    for csv_file in list_csv_files(oos_dir):
        try:
            df = pd.read_csv(csv_file)
            if "FOV_region" not in df.columns or "PARnorm" not in df.columns:
                continue

            df["FOV_region"] = df["FOV_region"].astype(str).str.strip()

            for region_str, sub_df in df.groupby("FOV_region"):
                coords = parse_FOV_region(region_str)
                if not coords:
                    continue

                parnorm_vals = pd.to_numeric(sub_df["PARnorm"], errors="coerce").dropna().tolist()
                if parnorm_vals:
                    region_to_parnorms[region_str].extend([float(v) for v in parnorm_vals])

                if region_str not in region_to_coords:
                    region_to_coords[region_str] = coords

        except Exception as e:
            print(f"Error processing {csv_file}: {e}")

    regions = []
    for region_id, coords in region_to_coords.items():
        vals = region_to_parnorms.get(region_id, [])
        parnorm_mean = float(np.mean(vals)) if vals else float("nan")
        regions.append({
            "id": region_id,
            **coords,
            "parNormMean": parnorm_mean,
        })

    regions.sort(key=lambda r: (r["yStart"], r["xStart"]))
    return {"regions": regions}


@app.get("/api/out-of-sample/{dataset_name}/{algorithm}/region/{region_id}")
async def get_region_data(dataset_name: str, algorithm: str, region_id: str):
    oos_dir = get_dataset_path(dataset_name) / "prediction" / algorithm / "out_of_sample_prediction"
    if not oos_dir.exists():
        raise HTTPException(status_code=404, detail="Out-of-sample prediction directory not found")

    normalized_region_id = normalize_region_id(region_id)
    matched_frames = []

    for csv_file in list_csv_files(oos_dir):
        try:
            df = pd.read_csv(csv_file)
            if "FOV_region" not in df.columns:
                continue

            df["FOV_region"] = df["FOV_region"].astype(str).str.strip()
            region_data = df[df["FOV_region"] == normalized_region_id]

            if not region_data.empty:
                matched_frames.append(region_data)

        except Exception as e:
            print(f"Error reading {csv_file}: {e}")

    if not matched_frames:
        raise HTTPException(status_code=404, detail=f"Region '{region_id}' not found")

    merged = pd.concat(matched_frames, ignore_index=True)

    if "gene" in merged.columns:
        merged = merged.drop_duplicates(subset=["gene"], keep="first")

    merged = merged.drop(columns=["FOV_region"], errors="ignore")
    merged = merged.where(pd.notnull(merged), None)

    return {
        "region_id": normalized_region_id,
        "data": merged.to_dict(orient="records"),
    }


@app.get("/api/out-of-sample/{dataset_name}/{algorithm}/region/{region_id}/download")
async def download_region_csv(dataset_name: str, algorithm: str, region_id: str):
    region_data = await get_region_data(dataset_name, algorithm, region_id)
    csv_content = pd.DataFrame(region_data["data"]).to_csv(index=False)
    return Response(
        content=csv_content,
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={dataset_name}_{algorithm}_{region_id}.csv"},
    )


@app.get("/api/out-of-sample/{dataset_name}/{algorithm}/region/{region_id}/image")
async def get_region_parnorm_image(
    dataset_name: str, algorithm: str, region_id: str,
    vmin: Optional[float] = None, vmax: Optional[float] = None,
):
    region_data = await get_region_data(dataset_name, algorithm, region_id)
    df = pd.DataFrame(region_data["data"])
    if "PARnorm" not in df.columns:
        raise HTTPException(status_code=404, detail="PARnorm column not found")
    vals = df["PARnorm"].to_numpy(dtype=float)
    data_2d = vals.reshape(1, -1)
    _vmin = vmin if vmin is not None else float(np.nanmin(vals))
    _vmax = vmax if vmax is not None else float(np.nanmax(vals))
    fig, ax = plt.subplots(figsize=(max(data_2d.shape[1] / 50, 4), 1), dpi=100)
    ax.imshow(data_2d, cmap="RdBu_r", aspect="auto", vmin=_vmin, vmax=_vmax, interpolation="nearest")
    ax.axis("off")
    plt.tight_layout(pad=0)
    buf = BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    buf.seek(0)
    return StreamingResponse(buf, media_type="image/png")


# ============================================================================
# 4.4 Different Initial Resolution Prediction
# ============================================================================

@app.get("/api/diff-res/{dataset_name}/{algorithm}/arrangements")
async def get_spot_arrangements(dataset_name: str, algorithm: str):
    diff_res_dir = diff_res_root(dataset_name, algorithm)
    if not diff_res_dir.exists():
        raise HTTPException(status_code=404, detail="Different resolution directory not found")
    rep1_dir = diff_res_dir / "rep1"
    if not rep1_dir.exists():
        raise HTTPException(status_code=404, detail="Replicate directory not found")
    return {"arrangements": sorted([d.name for d in rep1_dir.iterdir() if d.is_dir()])}


@app.get("/api/diff-res/{dataset_name}/{algorithm}/{arrangement}/resolutions")
async def get_resolutions(dataset_name: str, algorithm: str, arrangement: str):
    rep1_arr_dir = diff_res_rep_arr_dir(dataset_name, algorithm, "rep1", arrangement)
    if not rep1_arr_dir.exists():
        raise HTTPException(status_code=404, detail="Replicate directory not found")
    resolutions = sorted([f.stem for f in list_csv_files(rep1_arr_dir)], key=extract_numeric_value)
    return {"resolutions": resolutions}


@app.get("/api/diff-res/{dataset_name}/{algorithm}/{metric}/plot-data")
async def get_metrics_plot_data(dataset_name: str, algorithm: str, metric: str):
    """
    Returns plot data with mean and stdev across rep1/rep2/rep3.
    Each data point: { axisValue, metricValue (mean), stdev }
    
    Axis value logic:
      - proportional → spotSize chart: raw numeric from filename (= spot size in μm)
      - proportional → distance chart: spot_size / (55/100) = spot_size * 100/55
      - dense → both charts: raw numeric
      - fixed_spotDistance → spotSize chart: raw numeric
      - fixed_spotSize → distance chart: raw numeric
    """
    diff_res_dir = diff_res_root(dataset_name, algorithm)
    if not diff_res_dir.exists():
        raise HTTPException(status_code=404, detail="Different resolution directory not found")

    def collect_with_std(arrangement: str, x_transform=None) -> List[Dict]:
        rep_data: Dict[float, List[float]] = defaultdict(list)

        for rep in ["rep1", "rep2", "rep3"]:
            arr_dir = diff_res_rep_arr_dir(dataset_name, algorithm, rep, arrangement)
            if not arr_dir.exists():
                continue

            for csv_file in list_csv_files(arr_dir):
                try:
                    df = pd.read_csv(csv_file)
                    if metric not in df.columns:
                        continue

                    metric_series = pd.to_numeric(df[metric], errors="coerce").dropna()
                    if metric_series.empty:
                        continue

                    raw_av = extract_numeric_value(csv_file.stem)
                    av = x_transform(raw_av) if x_transform else raw_av
                    av = round(av, 6)

                    val = float(metric_series.mean())
                    if math.isfinite(val):
                        rep_data[av].append(val)

                except Exception as e:
                    print(f"Error reading {csv_file}: {e}")

        pts = []
        for av, vals in rep_data.items():
            finite_vals = [v for v in vals if math.isfinite(v)]
            if not finite_vals:
                continue

            mean_val = float(np.mean(finite_vals))
            std_val = float(np.std(finite_vals, ddof=1)) if len(finite_vals) > 1 else 0.0

            if not math.isfinite(mean_val):
                continue
            if not math.isfinite(std_val):
                std_val = 0.0

            pts.append({
                "axisValue": float(av),
                "metricValue": mean_val,
                "stdev": std_val,
            })

        pts.sort(key=lambda x: x["axisValue"])
        return pts

    # proportional distance axis: spot_size / (55/100) = spot_size * 100/55
    proportional_to_distance = lambda av: av / (55.0 / 100.0)

    return {
        "proportional_spotSize":       collect_with_std("proportional"),
        "proportional_distance":        collect_with_std("proportional", proportional_to_distance),
        "dense_spotSize":               collect_with_std("dense"),
        "dense_distance":               collect_with_std("dense"),
        "fixed_spotDistance_spotSize":  collect_with_std("fixed_spotDistance"),
        "fixed_spotSize_distance":      collect_with_std("fixed_spotSize"),
    }



@app.get("/api/diff-res/{dataset_name}/{algorithm}/{arrangement}/{resolution}/replicates")
async def get_replicate_data(dataset_name: str, algorithm: str, arrangement: str, resolution: str):
    diff_res_dir = diff_res_root(dataset_name, algorithm)
    if not diff_res_dir.exists():
        raise HTTPException(status_code=404, detail="Different resolution directory not found")

    result: Dict[str, List[Dict]] = {}
    for rep in ["rep1", "rep2", "rep3"]:
        csv_path = diff_res_rep_arr_dir(dataset_name, algorithm, rep, arrangement) / f"{resolution}.csv"
        if not csv_path.exists():
            result[rep] = []
            continue
        try:
            df = pd.read_csv(csv_path)
            df = df.where(pd.notnull(df), None)
            result[rep] = df.to_dict(orient="records")
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")
            result[rep] = []
    return result


@app.get("/api/diff-res/{dataset_name}/{algorithm}/{arrangement}/{resolution}/genes")
async def get_genes_for_resolution(dataset_name: str, algorithm: str, arrangement: str, resolution: str):
    csv_path = diff_res_rep_arr_dir(dataset_name, algorithm, "rep1", arrangement) / f"{resolution}.csv"
    if not csv_path.exists():
        raise HTTPException(status_code=404, detail="CSV file not found")
    return {"genes": read_genes_from_resolution_csv(csv_path)}


@app.get("/api/diff-res/{dataset_name}/{algorithm}/{arrangement}/{resolution}/gene/{gene}")
async def get_gene_metrics_and_spatial(dataset_name: str, algorithm: str, arrangement: str, resolution: str, gene: str):
    csv_path = diff_res_rep_arr_dir(dataset_name, algorithm, "rep1", arrangement) / f"{resolution}.csv"
    if not csv_path.exists():
        raise HTTPException(status_code=404, detail="CSV file not found")
    df = pd.read_csv(csv_path)
    gene_col = df.columns[0]
    gene_data = df[df[gene_col] == gene]
    if gene_data.empty:
        raise HTTPException(status_code=404, detail=f"Gene '{gene}' not found in CSV")
    metrics = gene_data.iloc[0].to_dict()

    pred_dir = metric_vis_prediction_dir(dataset_name)
    if not pred_dir.exists():
        raise HTTPException(status_code=404, detail="prediction directory not found")
    spatial_data = get_gene_spatial_data(pred_dir, gene, load_gene_index(pred_dir))

    ref_dir = metric_vis_reference_dir(dataset_name)
    reference_spatial = None
    if ref_dir.exists():
        try:
            reference_spatial = get_gene_spatial_data(ref_dir, gene, load_gene_index(ref_dir))
        except Exception:
            pass

    return {"gene": gene, "metrics": metrics, "prediction": spatial_data, "reference": reference_spatial}


# ============================================================================
# Health Check
# ============================================================================

@app.get("/api/health")
async def health_check():
    return {
        "status": "healthy",
        "data_root": str(DATA_ROOT),
        "data_root_exists": DATA_ROOT.exists(),
        "datasets_count": len(list_datasets()),
    }


@app.get("/")
async def root():
    return {"message": "SpatialVis API", "version": "1.0.0", "docs": "/docs"}


if __name__ == "__main__":
    import uvicorn
    port = int(os.getenv("PORT", 7001))
    host = os.getenv("HOST", "0.0.0.0")
    uvicorn.run(app, host=host, port=port)
