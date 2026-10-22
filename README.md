# SRBench

Code accompanying **Cost-Effective High-Definition Spatial Transcriptomics via Histology-Based Super-Resolution**.

SRBench provides data preparation, method execution, evaluation and analysis code for studying histology-guided spatial transcriptomics super-resolution. The study generates controlled low-resolution (LR) inputs from high-resolution (HR) data and evaluates predictions against matched HR references. It examines in-sample prediction, prediction in unmeasured regions and adjacent sections, and changes in input resolution, spot geometry and transcript detection.

The package runs **XFuse, TESLA, iSTAR, scstGCN and SuperST** through independent scripts. SuperST produces two outputs, `SuperST_t` and `SuperST_t_g`. A batch runner calls method scripts sequentially. Two standalone tools compare dataset characteristics and histological image features.

## What is included?

| Component | Purpose | Inputs and outputs | Instructions |
| --- | --- | --- | --- |
| `pipeline/preprocessing/` | Image/coordinate alignment and LR, FOV, resolution and dropout input generation for Xenium, Visium HD and CosMx | Platform data and registration settings → counts, coordinates, images and method input files | [Preprocessing](docs/PREPROCESSING.md) |
| `methods/` | Five method launchers, study implementation modifications and prediction export | Prepared inputs, upstream installations and model settings → gene-expression maps | [Method commands](docs/METHODS.md) |
| `run_methods.py` | Sequential execution of configured method jobs | JSON command list → predictions, per-step logs and a run summary | Batch usage below |
| `pipeline/metrics/` | PCC, SCC, KLD, JSD, SSIM and PSNR calculation | Prediction maps and matched HR references → per-gene metric CSVs | Metric usage below |
| `pipeline/analysis/` | Replicate summaries, rankings, robustness, gene-level effects, FOVs, adjacent sections, resolution, usability, gene selection and supplied post hoc tables | Metric CSVs or published numerical workbooks → summary CSVs, figures and workbooks | [Analysis guide](pipeline/analysis/README.md) |
| `tools/dataset_similarity/` | Dataset characterization and comparison with reference datasets | SRT data or feature tables → dataset features and ranked references | [Dataset tool](tools/dataset_similarity/README.md) |
| `tools/fov_similarity/` | UNI image-feature extraction and comparison with reference training FOVs | H&E images, UNI weights and reference features → ranked regions, features and previews | [FOV tool](tools/fov_similarity/README.md) |

The archive includes code, example configurations and documentation. Obtain data, Source Data workbooks, upstream method installations, model weights and reference feature collections separately. Data access follows the manuscript's Data Availability statement. [Input requirements](docs/DATA_REQUIREMENTS.md) lists the resources needed by each workflow.

## Where should I start?

| What you already have | Starting point |
| --- | --- |
| Raw platform data and histology | Register the data and generate LR inputs using platform preprocessing |
| Aligned counts, coordinates and an H&E image | Package method inputs with `prepare_inputs` |
| Prepared method inputs | Install the upstream method and study modifications, then use its launcher |
| Predicted expression maps and matched HR data | Calculate metrics, then run the relevant analysis |
| Published Source Data workbooks | Export numerical tables for supported downstream analyses without retraining |
| A new dataset or tissue image to compare with reference data | Use a standalone similarity tool |

The model workflow is **prepare inputs → run methods → calculate metrics → analyze results**. These are separate steps. The batch runner handles method execution; it does not automatically register images, generate all simulated conditions or run all downstream analyses.

## Installation and paths

Run examples from the extracted package root unless a command changes directory. `DATA/` denotes user-supplied data and `RESULTS/` denotes generated outputs. Replace placeholder paths and dataset-specific settings with actual values. Shell examples use Linux/Bash. Activate environments using the appropriate command on other systems.

Use **Python 3.12** for numerical analyses:

```bash
python -m venv .venv-analysis
source .venv-analysis/bin/activate
python -m pip install -r requirements-analysis.txt
python -m pipeline.analysis.cli --help
```

| Task | Environment/dependencies |
| --- | --- |
| Numerical analyses and figures | `requirements-analysis.txt`; locally checked versions are in `requirements-analysis-tested.txt` |
| Python preprocessing | `requirements-preprocessing.txt` |
| Notebook execution | `requirements-notebooks.txt` and platform-specific dependencies |
| Prediction metrics | `requirements-metrics.txt`, including Scanpy/PyArrow for Visium HD |
| Prediction methods | Each original method's own environment and GPU dependencies |
| Standalone tools | Each tool's `requirements.txt`, with R/HEARTSVG or UNI weights where required |

The analysis environment does not install all five models. iSTAR and XFuse launchers require Bash. Visium HD SVG selection uses R/IRkernel and HEARTSVG; CosMx registration uses an interactive napari environment. See [environment notes](docs/ENVIRONMENTS.md) for installation boundaries and the dependency requirements.

## Prepare method inputs

The platform notebooks/scripts in `pipeline/preprocessing/` align histology and expression coordinates and generate simulated inputs. Configure paths, physical pixel sizes, selected genes and tissue/FOV bounds before execution. Visually check the registration before training.

If aligned simulated counts are already available, package them as follows:

```bash
python -m pip install -r requirements-preprocessing.txt
python -m pipeline.preprocessing.prepare_inputs \
  --counts DATA/cnts.tsv \
  --locations DATA/locs-raw.tsv \
  --image DATA/he-raw.tiff \
  --output-dir DATA/prepared \
  --pixel-size 1.7 --spot-diameter 55 \
  --dataset Xe_hsa_breast_c \
  --gene-list DATA/gene-names.txt
```

- `cnts.tsv`: headered **spot-by-gene** matrix, with spot identifiers in column 1.
- `locs-raw.tsv`: the same spot identifiers and `x`, `y` image-pixel coordinates.
- `gene-names.txt`: the intended genes, one per line.
- `--pixel-size`: micrometres per image pixel; `--spot-diameter`: simulated spot diameter or square side length in micrometres. The values above illustrate one configuration and must match the actual data.

The packager aligns barcodes and gene metadata and exports matrix-format files, histology, `filtered_feature_bc_matrix.h5`, position tables and a `spatial/` directory. It packages existing counts; it does not perform registration or generate the LR simulation itself. Use an empty output directory and separate inputs for each method, condition and repeat. XFuse needs its own H5 conversion; SuperST also has an image-size preparation step. See [preprocessing](docs/PREPROCESSING.md) and [input schemas](docs/DATA_REQUIREMENTS.md).

## Run one method

Install the chosen upstream method in its own environment. Before running iSTAR or scstGCN, or converting XFuse inputs, install the corresponding study files:

```bash
python methods/install_patches.py --method iSTAR --method-dir /path/to/istar --apply
python methods/install_patches.py --method scstGCN --method-dir /path/to/scstGCN --apply
python methods/install_patches.py --method XFuse --method-dir /path/to/xfuse --apply
```

Run the installation command for the method you need. Original files are backed up as `*.srbench-original`; use `--verify` to check installation. See the [study modification guide](methods/patches/README.md).

For example, activate the iSTAR environment, then run from its upstream checkout:

```bash
cd /path/to/istar
bash /path/to/SRBench/methods/iSTAR_run.sh /path/to/DATA/iSTAR/rep1
```

The launcher runs image preparation, feature extraction and prediction. It keeps the supplied gene list and writes gene maps under the input's `cnts-super/`.

| Method | Launcher | Outputs/settings |
| --- | --- | --- |
| iSTAR | `methods/iSTAR_run.sh` | `cnts-super/`; 400 epochs |
| scstGCN | `methods/scstGCN_run.py` | `cnts-super/`; 500 epochs; set `UNI_CHECKPOINT` or `--uni-checkpoint` |
| TESLA | `methods/TESLA_run.py` | `tesla_enhanced/`; supply the actual patch size |
| XFuse | `methods/XFuse_run.sh`, then `methods/xfuse_export.py` | XFuse analysis output, then mean expression maps; use the study TOML |
| SuperST | `methods/SuperST_run.py` | `SuperST_enhanced_t/` and `SuperST_enhanced_t_g/`; supply the actual target size |

[Method commands](docs/METHODS.md) provides complete commands, required files and working directories for all five methods. Preserve the actual gene lists, conversion scales, model settings and evaluation domains when reproducing results.

## Run several methods

Return to the package root, copy the configuration and edit its paths/settings:

```bash
cp configs/method_runs.example.json configs/my_runs.json
# Edit project_root, installations, environment names, inputs and outputs.
python run_methods.py --config configs/my_runs.json --dry-run
python run_methods.py --config configs/my_runs.json --check
python run_methods.py --config configs/my_runs.json
```

`--dry-run` prints commands without executing models or requiring example data paths to exist. `--check` checks working directories, declared inputs and output destinations; it does not load models or verify GPU compatibility. With neither flag, jobs run sequentially. Each command selects its own environment; the supplied configuration uses `conda run`.

A job represents a method/dataset/condition/repeat. Use unique IDs and separate run directories. XFuse uses two steps, prediction and export. Replace example TESLA patch sizes and SuperST target sizes with the intended values.

```bash
python run_methods.py --config configs/my_runs.json --methods iSTAR TESLA
python run_methods.py --config configs/my_runs.json --keep-going
```

`--methods` selects configured methods. A failure stops that job's remaining steps; `--keep-going` permits subsequent independent jobs. Existing declared outputs are protected. The configured `log_dir` (`RUN_LOGS/` by default) contains per-step logs and `summary.json` with commands, exit codes, elapsed times and errors.

In the JSON, `project_root` is relative to the configuration file. Relative `cwd`, `inputs`, `outputs` and `log_dir` resolve against that project root. Command arguments use `{root}` for the project root and `{package}` for the code package. Analysis configurations may use different path conventions; consult their guides.

## Calculate prediction metrics

Invoke the platform evaluation scripts **as file commands** from the package root:

```bash
python -m pip install -r requirements-metrics.txt
python pipeline/metrics/Xenium_eval.py --help
python pipeline/metrics/VisiumHD_eval.py --help
```

For aligned transcript-coordinate data, evaluate iSTAR predictions with:

```bash
python pipeline/metrics/Xenium_eval.py \
  --he-path DATA/he-raw.tiff \
  --transcript-use-path DATA/transcripts_use.tsv \
  --algorithm iSTAR \
  --prediction-path DATA/iSTAR/rep1 \
  --output-path RESULTS/metrics
```

This route also handles prepared CosMx transcript-coordinate data. Set the gene and coordinate column arguments if their order differs from `feature_name, x, y`. The image and transcript coordinates must describe the same registered region. For dropout evaluation, use the unchanged HR reference, rather than the subset used to generate the LR input.

Visium HD uses the 8 µm reference H5, positions parquet and a finalized range JSON:

```bash
python pipeline/metrics/VisiumHD_eval.py evaluate \
  --algorithm iSTAR \
  --prediction-path DATA/iSTAR/rep1 \
  --gt-data DATA/filtered_feature_bc_matrix.h5 \
  --tissue-positions DATA/tissue_positions.parquet \
  --range-info DATA/finalized_ranges.json \
  --output-path RESULTS/metrics
```

`detect-range` generates candidate ranges for manual inspection; it does not replace alignment review. The XFuse range tool supports the explicit study transform/cropping step. Review [evaluation inputs](docs/DATA_REQUIREMENTS.md) before evaluating cropped predictions or reusing ground-truth caches.

Use the metric-to-analysis file layout in [the analysis guide](pipeline/analysis/README.md#metric-file-layout) before running `insample`.

Outputs are per-gene CSVs with `PCC`, `SCC`, `KLD`, `JSD`, `SSIM`, `PSNR` and `gt_counts`. Keep separate per-run files to analyze computational variability.

## Analyze metrics or published Source Data

```bash
python -m pipeline.analysis.cli --help
```

For raw metric CSVs, use `insample` to summarize three computational runs, then `rank-insample`. For published gene-level Source Data, start from the workbook:

```bash
python -m pipeline.analysis.cli export-source \
  --workbook DATA/Gu_Supplementary_Data_SuppFig6_95.xlsx \
  --output-dir RESULTS/source_tables

python -m pipeline.analysis.cli hvg \
  --input-dir RESULTS/source_tables/classified \
  --dataset Xe_hsa_breast_c --methods iSTAR scstGCN \
  --output-dir RESULTS/hvg
```

The exporter writes `merged/` metric tables and `classified/` tables where annotations are available. Published values are already averaged across runs; exporting them does not reconstruct the original repeats or their variability.

| Analysis | Entry point/module | Outputs |
| --- | --- | --- |
| In-sample accuracy, rankings and robustness | CLI `insample`, `rank-insample` | Replicate summaries, rankings and figures |
| HVG/non-HVG and gene-property effects | CLI `hvg`, `strata`, `gene-effects` | Statistical CSVs, heatmaps and expression/spatial-variability contrasts |
| FOVs and adjacent sections | CLI `fov`, `adjacent` | Ranked FOV tables, maps, comparisons and paired tests |
| Resolution and spot geometry | CLI `resolution` | Replicate means, ranks and mean performance curves |
| Gene selection | CLI `top-genes`; `pipeline.analysis.gene_selection` | Selected genes and Supplementary Tables 3–5 with complete inputs |
| Post hoc absolute in-sample summaries | `pipeline.analysis.posthoc.build_insample_absolute_summaries` | Supplementary Table 2 |
| Dropout summaries | `pipeline.analysis.posthoc.build_dropout_replicate_means`, `build_dropout_paired_statistics`, `build_dropout_supplement_a4` | Seed-level summaries, Supplementary Table 6 and figures |
| Usability | CLI `usability` | Normalized observations and a comparison figure |

The [analysis guide](pipeline/analysis/README.md) supplies full commands, input layouts, dataset/method identifiers and output definitions. Full-study table builders require the stated dataset/method collections. Individual analysis tasks can use the relevant subsets where supported.

## Use the two standalone tools

Install the chosen tool's dependencies in its own environment. Supply reference data separately; neither tool requires running all prediction methods first.

**Dataset similarity:** extract characteristics with `compute_dataset_metrics.py`, then compare query and reference feature tables:

```bash
python tools/dataset_similarity/find_similar_dataset.py \
  --query-csv DATA/query_features.csv \
  --candidate-csv DATA/reference_features.csv \
  --topk 3 --output-csv RESULTS/similar_datasets.csv
```

The [dataset tool guide](tools/dataset_similarity/README.md) covers extraction, R/HEARTSVG requirements and comparison settings.

**Histological FOV similarity:** build reference features with `extract_top3_uni_features.py` or extract a region with `extract_single_region_feature.py`, then search a user image:

```bash
python tools/fov_similarity/search_similar_regions.py \
  --npy_path DATA/reference_features.npy \
  --he_image_path DATA/query_he.tiff \
  --ckpt_path CHECKPOINTS/UNI/pytorch_model.bin \
  --fov_width 1500 --fov_height 1500 \
  --output_dir RESULTS/similar_regions
```

Window dimensions above are examples in pixels. Query and reference features require compatible physical scales and the same checkpoint/preprocessing. Image-feature similarity does not establish prediction accuracy on a new tissue. See the [FOV tool guide](tools/fov_similarity/README.md).

## Reproduction, integrity and licensing

Reproducing a result requires its actual data, registrations, selected genes, method settings and evaluation domains. The batch runner does not infer these. Clinical data and third-party models retain their respective access conditions.

Check the extracted files on Linux with:

```bash
sha256sum -c MANIFEST.sha256
```

Original SRBench contributions use [MIT](LICENSE). Adapted iSTAR/scstGCN files retain GPL-3.0; XFuse files retain upstream MIT notices. See [third-party notices](THIRD_PARTY_NOTICES.md). Cite the manuscript, archived software DOI and original methods as applicable.
