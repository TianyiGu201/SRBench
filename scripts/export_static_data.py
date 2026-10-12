"""Export the small API-derived files needed by the Pages + object storage build.

The original NPY, CSV, JPEG, TIFF, and metadata JSON files stay in data_root.
This script writes only the derived JSON under web_data/site.
"""

import argparse
import asyncio
import json
import math
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from backend import main  # noqa: E402


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")))


def export_oos(dataset: str, algorithm: str, source: Path, output: Path) -> None:
    csv_files = sorted(source.glob("*.csv"))
    if not csv_files:
        return
    frames = []
    for csv_file in csv_files:
        frame = pd.read_csv(csv_file)
        if "FOV_region" in frame.columns:
            frame["FOV_region"] = frame["FOV_region"].astype(str).str.strip()
            frames.append(frame)
    if not frames:
        return

    merged = pd.concat(frames, ignore_index=True)
    root = output / "oos" / dataset / algorithm
    regions = []
    for region_id, frame in merged.groupby("FOV_region", sort=False):
        coords = main.parse_FOV_region(region_id)
        if coords is None:
            continue
        parnorm = pd.to_numeric(frame.get("PARnorm"), errors="coerce").dropna()
        mean = float(parnorm.mean()) if not parnorm.empty else None
        if mean is not None and not math.isfinite(mean):
            mean = None
        regions.append({"id": region_id, **coords, "parNormMean": mean})
        if "gene" in frame.columns:
            frame = frame.drop_duplicates(subset=["gene"], keep="first")
        frame = frame.drop(columns=["FOV_region"], errors="ignore")
        records = json.loads(frame.to_json(orient="records", double_precision=15))
        write_json(root / "region" / f"{region_id}.json", {"region_id": region_id, "data": records})
    regions.sort(key=lambda region: (region["yStart"], region["xStart"]))
    write_json(root / "regions.json", {"regions": regions})
    print(f"OOS  {dataset}/{algorithm}: {len(regions)} regions", flush=True)


def export_diff(dataset: str, algorithm: str, source: Path, output: Path) -> None:
    if not (source / "rep1").exists():
        return
    arrangements = sorted(path.name for path in (source / "rep1").iterdir() if path.is_dir())
    resolutions = {}
    for arrangement in arrangements:
        files = (source / "rep1" / arrangement).glob("*.csv")
        resolutions[arrangement] = sorted((path.stem for path in files), key=main.extract_numeric_value)
    root = output / "diff" / dataset / algorithm
    write_json(root / "manifest.json", {"arrangements": arrangements, "resolutions": resolutions})
    for metric in ("PCC", "SCC", "KLD", "JSD"):
        plot = asyncio.run(main.get_metrics_plot_data(dataset, algorithm, metric))
        write_json(root / "plots" / f"{metric}.json", plot)
    print(f"DIFF {dataset}/{algorithm}: {len(arrangements)} arrangements", flush=True)


def main_cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--datasets", nargs="*")
    parser.add_argument("--algorithms", nargs="*")
    args = parser.parse_args()
    main.DATA_ROOT = args.data_root.resolve()
    output = args.output.resolve()
    datasets = asyncio.run(main.get_datasets())
    for item in datasets:
        item["heImageUrl"] = f"/data_root/{item['id']}/he-raw.jpg"
    write_json(output / "datasets.json", datasets)
    selected = set(args.datasets or [item["id"] for item in datasets])
    algorithms = set(args.algorithms or ("XFuse", "TESLA", "iSTAR", "scstGCN", "SuperST_t", "SuperST_t_g"))
    for item in datasets:
        dataset = item["id"]
        if dataset not in selected:
            continue
        prediction = main.DATA_ROOT / dataset / "prediction"
        if not prediction.exists():
            continue
        for algorithm_dir in sorted(prediction.iterdir()):
            if not algorithm_dir.is_dir() or algorithm_dir.name not in algorithms:
                continue
            algorithm = algorithm_dir.name
            export_oos(dataset, algorithm, algorithm_dir / "out_of_sample_prediction", output)
            export_diff(dataset, algorithm, algorithm_dir / "raw_resolution", output)
    print(f"Exported derived JSON to {output}", flush=True)


if __name__ == "__main__":
    main_cli()
