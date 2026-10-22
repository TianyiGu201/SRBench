# Study modifications to the upstream methods

The files here contain the implementation modifications described in Supplementary Note 1. Obtain and install each original method separately, then install the study files before running its launcher. Use the matching method environment.

```bash
python SRBENCH/methods/install_patches.py --method iSTAR --method-dir UPSTREAM/istar --apply
python SRBENCH/methods/install_patches.py --method scstGCN --method-dir UPSTREAM/scstGCN --apply
python SRBENCH/methods/install_patches.py --method XFuse --method-dir UPSTREAM/xfuse --apply
```

The installer preserves each original file as `*.srbench-original`. Use `--verify` instead of `--apply` to check installed files without changing them. iSTAR and scstGCN launchers verify their installed study files before starting.

| Supplied file | Destination relative to the upstream source checkout |
| --- | --- |
| `iSTAR/impute.py` | `impute.py` |
| `scstGCN/impute.py` | `impute.py` |
| `scstGCN/histological_feature_map.py` | `histological_feature_map.py` |
| `XFuse/visium.py` | `xfuse/convert/visium.py` |
| `XFuse/utility.py` | `xfuse/convert/utility.py` |

For XFuse, use an editable installation of this checkout (`python -m pip install -e UPSTREAM/xfuse`) so that its converter uses the supplied files. Apply these files before generating `data.h5`. Prediction from an already prepared `data.h5` does not invoke the converter again.

iSTAR and scstGCN retain failed patch positions as NaNs until the features, expression targets and coordinates are filtered together. iSTAR records patch errors under the corresponding state directory; scstGCN uses the run directory. The supplied scstGCN square-patch extraction is retained.

Set `UNI_CHECKPOINT` to the local UNI checkpoint used for feature extraction, or pass `--uni-checkpoint` to `scstGCN_run.py`. Checkpoints are obtained separately.

The iSTAR/scstGCN files retain GPL-3.0 notices. The XFuse files retain the upstream MIT notice in `licenses/XFuse-MIT.txt`. Original SRBench additions use MIT.
