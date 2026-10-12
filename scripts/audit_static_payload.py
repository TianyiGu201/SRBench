"""Check that every file referenced by the static website exists before R2 upload."""

import argparse
import json
from pathlib import Path


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def collect_references(root: Path) -> tuple[set[Path], set[Path], int, int, int]:
    """Return required and optional source paths used by the static frontend."""
    data_root = root / "data_root"
    site_root = root / "site"
    optional_missing: set[Path] = set()
    referenced: set[Path] = set()

    def require(path: Path) -> None:
        referenced.add(path)

    def optional(path: Path) -> None:
        referenced.add(path)
        if not path.is_file():
            optional_missing.add(path)

    catalog = site_root / "datasets.json"
    require(catalog)
    if not catalog.is_file():
        raise SystemExit(f"Missing dataset catalog: {catalog}")

    datasets = read_json(catalog)
    for dataset in datasets:
        dataset_id = dataset["id"]
        folder = data_root / dataset_id
        for name in ("dataset_metaINFO.json", "metric_vis.json", "he-raw.jpg", "he-raw.tiff"):
            require(folder / name)
        metric_file = folder / "metric_vis.json"
        if metric_file.is_file():
            for levels in read_json(metric_file).values():
                for genes in levels.values():
                    for gene in genes:
                        for kind in ("reference", "prediction"):
                            require(folder / "metric_visualization" / kind / f"{gene}.npy")

    oos_count = 0
    for regions_file in site_root.glob("oos/*/*/regions.json"):
        require(regions_file)
        for region in read_json(regions_file)["regions"]:
            require(regions_file.parent / "region" / f"{region['id']}.json")
            oos_count += 1

    diff_count = 0
    for manifest_file in site_root.glob("diff/*/*/manifest.json"):
        require(manifest_file)
        dataset_id = manifest_file.parent.parent.name
        algorithm = manifest_file.parent.name
        manifest = read_json(manifest_file)
        for metric in ("PCC", "SCC", "KLD", "JSD"):
            require(manifest_file.parent / "plots" / f"{metric}.json")
        for arrangement, resolutions in manifest["resolutions"].items():
            for resolution in resolutions:
                for rep in ("rep1", "rep2", "rep3"):
                    path = data_root / dataset_id / "prediction" / algorithm / "raw_resolution" / rep / arrangement / f"{resolution}.csv"
                    if rep == "rep1":
                        require(path)
                    else:
                        optional(path)
                    diff_count += 1

    return referenced - optional_missing, optional_missing, len(datasets), oos_count, diff_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--web-data", type=Path, required=True)
    args = parser.parse_args()
    root = args.web_data.resolve()
    data_root = root / "data_root"
    site_root = root / "site"
    referenced, optional_missing, dataset_count, oos_count, diff_count = collect_references(root)
    missing = {path for path in referenced if not path.is_file()}
    upload_files = [p for folder in (data_root, site_root) for p in folder.rglob("*") if p.is_file() and p.suffix != ".bak"]
    upload_bytes = sum(p.stat().st_size for p in upload_files)
    print(f"Datasets: {dataset_count}")
    print(f"Referenced files: {len(referenced) + len(optional_missing)}; OOS regions: {oos_count}; replicate CSVs: {diff_count}")
    print(f"Upload scope excluding .bak: {len(upload_files)} files, {upload_bytes / 1_000_000_000:.2f} GB")
    print(f"Optional repeat CSVs absent in source: {len(optional_missing)} (shown as empty results)")
    if missing:
        print(f"Missing referenced files: {len(missing)}")
        for path in sorted(missing)[:20]:
            print(f"  {path.relative_to(root)}")
        raise SystemExit(1)
    print("All referenced files exist")


if __name__ == "__main__":
    main()
