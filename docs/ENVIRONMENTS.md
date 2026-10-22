# Software environments

## Analysis environment

Use Python 3.12 and `requirements-analysis.txt` for CSV/Source Data analyses,
statistical summaries and plots. `requirements-analysis-tested.txt` records the
package versions for numerical analysis. Prediction models use separate environments.

```bash
python -m venv .venv
# Activate the environment using your operating system's normal command.
python -m pip install -r requirements-analysis.txt
python -m pipeline.analysis.cli --help
```

Use this environment for numerical analysis. Prediction methods require their own environments.

## Platform preprocessing

Install `requirements-preprocessing.txt` for the noninteractive Python
preprocessing and format-export steps. PyTorch supports CPU aggregation.
Visium HD retains sparse counts and has a CPU path without CuPy/Numba;
GPU acceleration additionally requires Numba and a CuPy distribution matched
to the CUDA toolkit. Install GPU-specific dependencies separately.

Gene characterization and the standalone tools have additional requirements,
including Scanpy/Squidpy/AnnData where specified by their individual workflows.

Notebook execution uses `requirements-notebooks.txt`. `1.selectSVG.ipynb`
requires an R/IRkernel environment with HEARTSVG and its dependencies, not a
Python kernel. CosMx alignment requires the napari/napari-cosmx desktop
environment and exported image-registration settings. See the platform notes
in `DATA_REQUIREMENTS.md` before running raw-data notebooks.

## Prediction environments

For prediction-metric evaluation, install `requirements-metrics.txt`. It adds
Scanpy, PyArrow and H5 support to the analysis dependencies. Run evaluator files
directly, for example `python pipeline/metrics/Xenium_eval.py --help`; see the
root README for full commands.

Install each method according to its original documentation in a separate
environment. Set its repository as the working directory or `--method-dir`,
and include the appropriate environment in its batch command. The five methods
have different PyTorch/CUDA requirements.

| Method | Repository | Retained study wrapper settings |
| --- | --- | --- |
| iSTAR | https://github.com/daviddaiweizhang/istar | pixel_size=0.5; 400 epochs; existing gene list retained |
| scstGCN | https://github.com/wenwenmin/scstGCN | 500 epochs; existing gene list retained; local UNI_CHECKPOINT |
| TESLA | https://github.com/jianhuupenn/TESLA | supplied patch size; s=1, k=2, num_nbs=10 |
| XFuse | https://github.com/ludvb/xfuse | supplied TOML and conversion scale; no automatic tuning |
| SuperST | https://github.com/portrai-io/SuperST | 320 iterations; is_cut_bg=True; explicit target_s; both outputs retained |

Install the supplied study modifications as described in `methods/patches/README.md`.
For scstGCN, set `UNI_CHECKPOINT` to the separately obtained local UNI checkpoint.
Authentication tokens and model checkpoints must not be committed to the repository.

## Standalone tools

Each tool has its own `requirements.txt`. Histology feature extraction requires
the selected UNI checkpoint, compatible torch/torchvision/timm versions, and a
reference feature vector generated with exactly the same checkpoint/preprocessing.
Dataset SVG characterization also requires HEARTSVG through `Rscript`.
