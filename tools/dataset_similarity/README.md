# Dataset similarity tool

This standalone tool characterizes user-supplied SRT data and compares the resulting features with reference dataset features. It does not predict super-resolved expression or guarantee performance on a new dataset.

Install this directory's requirements in a suitable Scanpy/Squidpy environment:

```bash
python -m pip install -r tools/dataset_similarity/requirements.txt
python tools/dataset_similarity/compute_dataset_metrics.py --help
python tools/dataset_similarity/find_similar_dataset.py --help
```

The computation input is a directory with one supported dataset per subdirectory. The comparison input is a user feature CSV plus a reference feature CSV. Use `--help` for the explicit comparison schemes and weights supported by the script. The reference feature CSV must be provided or regenerated; it is not embedded in this code archive.

Example feature extraction:

```bash
python tools/dataset_similarity/compute_dataset_metrics.py --dataset-root DATA/user_datasets --output-csv RESULTS/user_dataset_features.csv --heartsvg-runner tools/dataset_similarity/run_heartsvg.R
```

Supply `--uni-ckpt` only when image-feature extraction is required. SVG computation uses the separately installed R/HEARTSVG environment through `Rscript`. `install_heartsvg.R` and `run_heartsvg.R` provide the R installation/entry scripts. Record the checkpoint and feature-table versions used for each comparison.
