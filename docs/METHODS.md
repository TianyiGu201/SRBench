# Running the methods

Run each command in the method's installed environment. Replace `SRBENCH`,
`UPSTREAM`, `INPUT` and `OUTPUT` with actual absolute paths. iSTAR and XFuse
launchers require Bash. Method installations, checkpoints and prepared inputs
are supplied separately. Install the supplied [study modifications](../methods/patches/README.md)
before running iSTAR/scstGCN or converting XFuse inputs.

## iSTAR

Run from the installed iSTAR repository:

```bash
cd UPSTREAM/istar
bash SRBENCH/methods/iSTAR_run.sh INPUT/iSTAR/rep1
```

Inputs are `he-raw.tiff`, `cnts.tsv`, `locs-raw.tsv`, `radius-raw.txt`,
`pixel-size-raw.txt` and the actual `gene-names.txt`. Predictions are written
under `INPUT/iSTAR/rep1/cnts-super/`. The script retains 400 training epochs
and a target image pixel size of 0.5 micrometres. Set `SRBENCH_PYTHON` if the
environment's interpreter has a different executable name. Optional `PID_FILE`
records process IDs for an external resource monitor.

## scstGCN

```bash
python SRBENCH/methods/scstGCN_run.py INPUT/scstGCN/rep1 --method-dir UPSTREAM/scstGCN
```

Supply the matrix-format input files and study gene list. Set `UNI_CHECKPOINT`
to the separately obtained local UNI checkpoint, or pass `--uni-checkpoint PATH`. The script retains
500 training epochs and exports predictions under `cnts-super/`.

## TESLA

```bash
python SRBENCH/methods/TESLA_run.py INPUT/TESLA/rep1 INPUT/TESLA/rep1/he-raw.tiff INPUT/TESLA/rep1 PATCH_SIZE INPUT/TESLA/rep1/gene-names.txt
```

Arguments are the Visium-format data directory, H&E image, working directory,
actual patch size and gene list. Use the same prepared directory as the working
directory to retain the evaluation layout. Predictions are exported under
`tesla_enhanced/`. Install the upstream TESLA package in this environment.

## XFuse

Prepare `data.h5` using the original XFuse conversion workflow and set its actual
location in the study TOML. Run from a working directory consistent with any
relative data paths in that TOML:

```bash
cd INPUT/XFuse/rep1
bash SRBENCH/methods/XFuse_run.sh INPUT/XFuse/rep1/config.toml OUTPUT/my-run
python SRBENCH/methods/xfuse_export.py --input-dir OUTPUT/my-run/analyses/final/analysis-gene_maps/section0
```

The second command exports posterior-sample mean maps using the supplied
conversion. Apply the finalized XFuse edge domain explicitly before metric
evaluation. The TOML supplied in `methods/` is a template; preserve the actual
study scale, parameters and input paths for reproduction.

## SuperST

```bash
python SRBENCH/methods/SuperST_run.py INPUT/SuperST/rep1 TARGET_S --method-dir UPSTREAM/SuperST
```

Supply `filtered_feature_bc_matrix.h5`, the matching `spatial/` directory and
the finalized image-size setup. `TARGET_S` must be the actual study setting.
The original 320-iteration call produces `SuperST_enhanced_t/` and
`SuperST_enhanced_t_g/`; both outputs are retained. Image-size preparation is
in `methods/prepare_superst_image_size.ipynb`.

## Batch command list

`run_methods.py` sequentially invokes the commands in a JSON file. A job has an
ID, method name, working directory, command, declared inputs and expected outputs:

```json
{
  "id": "scstGCN_dataset_rep1",
  "method": "scstGCN",
  "cwd": "DEPENDENCIES/scstGCN",
  "command": ["conda", "run", "--no-capture-output", "-n", "scstgcn", "python", "{package}/methods/scstGCN_run.py", "{root}/DATA/scstGCN/rep1", "--method-dir", "{root}/DEPENDENCIES/scstGCN"],
  "inputs": ["DATA/scstGCN/rep1/gene-names.txt"],
  "outputs": ["DATA/scstGCN/rep1/cnts-super"]
}
```

Use `steps` instead of `command` for a job with multiple sequential commands.
The full example includes all five methods and the required model input paths.
An optional `env` object adds local environment variables to that job. Avoid
storing authentication credentials in configuration files.

`--dry-run` prints commands. `--check` checks declared paths; it does not load
model libraries. During execution, failures and exit codes are saved to a unique
log directory. A failed job stops its remaining steps; `--keep-going` allows the
other jobs to continue. A successful command must produce its declared outputs.
Existing outputs require a new run directory. Computational repeats require
separate jobs and, for dropout, the corresponding pre-generated sampled inputs.
