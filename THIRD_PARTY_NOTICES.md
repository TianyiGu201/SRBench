# License scope and upstream software

The root MIT license applies to original SRBench contributions. It does not
replace an upstream license or authorize redistribution of model weights or data.

## Adapted files included in this archive

| File | Origin | Applicable notice |
| --- | --- | --- |
| `methods/iSTAR_run.sh` | Adapted from the iSTAR execution workflow | GPL-3.0; `licenses/iSTAR-GPL-3.0.txt` |
| `methods/scstGCN_run.py` | Adapted from the scstGCN execution workflow | GPL-3.0; `licenses/scstGCN-GPL-3.0.txt` |
| `methods/patches/iSTAR/impute.py` | Supplied iSTAR implementation with study patch handling | GPL-3.0; `licenses/iSTAR-GPL-3.0.txt` |
| `methods/patches/scstGCN/*.py` | Supplied scstGCN implementation and local feature-model loading | GPL-3.0; `licenses/scstGCN-GPL-3.0.txt` |
| `methods/patches/XFuse/*.py` | Supplied XFuse conversion modifications | MIT; `licenses/XFuse-MIT.txt` |
| `methods/TESLA_run.py` | TESLA-based execution and SRBench output export | Original TESLA MIT notice retained in `licenses/TESLA-MIT.txt`; original SRBench changes use MIT |
| `methods/xfuse_config.example.toml` | XFuse configuration template | Original XFuse MIT notice in `licenses/XFuse-MIT.txt`; SRBench settings do not change the upstream license |

The adapted launchers and iSTAR/scstGCN implementation files retain their GPL
notices. They are not relicensed under MIT.
The original SRBench orchestration, simulations, analysis code and standalone tools
use the root MIT license unless a file explicitly carries a different notice.

## External implementations and weights

- iSTAR: https://github.com/daviddaiweizhang/istar (GPL-3.0).
- scstGCN: https://github.com/wenwenmin/scstGCN (GPL-3.0).
- TESLA: https://github.com/jianhuupenn/TESLA (MIT).
- XFuse: https://github.com/ludvb/xfuse (MIT).
- SuperST: https://github.com/portrai-io/SuperST. No explicit repository license
  was located in the license audit on 2026-10-22. Its algorithm implementation is
  not bundled; the included SRBench launcher invokes a separately obtained copy.
- UNI and other third-party checkpoints are external resources governed by their
  providers' access and license terms. No checkpoint or authentication token is included.

License notices were retrieved from upstream repositories; `licenses/SOURCES.json`
records their URLs and SHA-256 digests. These retrievals identify license sources,
not the historical commits used for the manuscript experiments. Cite the original
methods and respect the terms attached to the particular versions installed.
