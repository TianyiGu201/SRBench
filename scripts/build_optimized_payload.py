"""Build and verify the lossless, browser-readable R2 payload in a new directory."""

import argparse
import gzip
import hashlib
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import tifffile

from audit_static_payload import collect_references


COMPRESS_SUFFIXES = {".npy", ".csv", ".json"}
COPY_SUFFIXES = {".jpg"}
TIFF_SUFFIXES = {".tiff"}


def destination(source: Path, web_data: Path, output: Path) -> Path:
    relative = source.relative_to(web_data)
    path = output / relative
    return path.with_name(path.name + ".gz") if source.suffix in COMPRESS_SUFFIXES else path


def hash_file(path: Path, compressed: bool = False) -> str:
    digest = hashlib.sha256()
    opener = gzip.open if compressed else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_chunk(files: list[str], root: str, output: str) -> tuple[int, int, int]:
    web_data = Path(root)
    stage = Path(output)
    original_bytes = staged_bytes = 0
    for name in files:
        source = Path(name)
        target = destination(source, web_data, stage)
        target.parent.mkdir(parents=True, exist_ok=True)
        raw = source.read_bytes()
        packed = gzip.compress(raw, compresslevel=6, mtime=0)
        target.write_bytes(packed)
        original_bytes += len(raw)
        staged_bytes += len(packed)
    return len(files), original_bytes, staged_bytes


def chunks_by_size(paths: list[Path], max_bytes: int = 100_000_000) -> list[list[str]]:
    chunks: list[list[str]] = []
    current: list[str] = []
    size = 0
    for path in paths:
        file_size = path.stat().st_size
        if current and size + file_size > max_bytes:
            chunks.append(current)
            current = []
            size = 0
        current.append(str(path))
        size += file_size
    if current:
        chunks.append(current)
    return chunks


def validate_source(web_data: Path) -> list[Path]:
    references, missing_optional, _, _, _ = collect_references(web_data)
    missing = [path for path in references if not path.is_file()]
    if missing:
        raise RuntimeError(f"{len(missing)} required files are missing; first: {missing[0]}")
    unsupported = [path for path in references if path.suffix not in COMPRESS_SUFFIXES | COPY_SUFFIXES | TIFF_SUFFIXES]
    if unsupported:
        raise RuntimeError(f"Unsupported source type: {unsupported[0]}")
    print(f"Referenced files to stage: {len(references)}; optional CSVs absent: {len(missing_optional)}", flush=True)
    return sorted(references)


def build(web_data: Path, output: Path, workers: int) -> None:
    files = validate_source(web_data)
    if output.exists():
        raise FileExistsError(f"Output already exists; choose a new directory: {output}")
    output.mkdir(parents=True)
    compressed = [path for path in files if path.suffix in COMPRESS_SUFFIXES]
    chunks = chunks_by_size(compressed)
    completed = source_bytes = staged_bytes = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(gzip_chunk, group, str(web_data), str(output)) for group in chunks]
        for future in as_completed(futures):
            count, before, after = future.result()
            completed += count
            source_bytes += before
            staged_bytes += after
            if completed % 2000 < count:
                print(f"gzip: {completed}/{len(compressed)} files", flush=True)
    print(f"gzip complete: {source_bytes} -> {staged_bytes} bytes", flush=True)

    for source in files:
        if source.suffix in TIFF_SUFFIXES:
            target = destination(source, web_data, output)
            target.parent.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(["tiffcp", "-c", "zip:2:p6", str(source), str(target)], capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(f"TIFF compression failed for {source}: {result.stderr}")
            print(f"TIFF {source.parent.name}: {source.stat().st_size} -> {target.stat().st_size} bytes", flush=True)
        elif source.suffix in COPY_SUFFIXES:
            target = destination(source, web_data, output)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)

    stage_files = [path for path in output.rglob("*") if path.is_file()]
    if len(stage_files) != len(files):
        raise RuntimeError(f"Staged file count mismatch: {len(stage_files)} vs {len(files)}")
    total = sum(path.stat().st_size for path in stage_files)
    print(f"Staging complete: {len(stage_files)} files, {total} bytes ({total / 1e9:.3f} GB)", flush=True)


def verify(web_data: Path, output: Path) -> None:
    files = validate_source(web_data)
    expected = {destination(path, web_data, output) for path in files}
    actual = {path for path in output.rglob("*") if path.is_file()}
    if actual != expected:
        raise RuntimeError(f"Staged inventory differs: missing={len(expected - actual)}, extra={len(actual - expected)}")
    for index, source in enumerate(files, 1):
        target = destination(source, web_data, output)
        if source.suffix in COMPRESS_SUFFIXES:
            if hash_file(source) != hash_file(target, compressed=True):
                raise RuntimeError(f"gzip round trip differs: {source}")
        elif source.suffix in TIFF_SUFFIXES:
            original = tifffile.memmap(source)
            packed = tifffile.imread(target)
            if not np.array_equal(original, packed):
                raise RuntimeError(f"TIFF pixels differ: {source}")
            with tifffile.TiffFile(source) as before, tifffile.TiffFile(target) as after:
                original_tags = {tag.name: tag.value for tag in before.pages[0].tags.values()}
                packed_tags = {tag.name: tag.value for tag in after.pages[0].tags.values()}
                exempt = {"Compression", "StripOffsets", "StripByteCounts", "RowsPerStrip", "Predictor"}
                if any(name not in packed_tags or packed_tags[name] != value for name, value in original_tags.items() if name not in exempt):
                    raise RuntimeError(f"TIFF metadata differs: {source}")
        elif hash_file(source) != hash_file(target):
            raise RuntimeError(f"Copied file differs: {source}")
        if index % 4000 == 0:
            print(f"Verified {index}/{len(files)} files", flush=True)
    total = sum(path.stat().st_size for path in actual)
    print(f"PASS: {len(files)} files validated, {total} bytes ({total / 1e9:.3f} GB)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "verify"))
    parser.add_argument("--web-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    root = args.web_data.resolve()
    output = args.output.resolve()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.mode == "build":
        build(root, output, args.workers)
    else:
        verify(root, output)


if __name__ == "__main__":
    main()
