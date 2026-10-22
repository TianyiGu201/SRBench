#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-only
set -e
# Adapted from iSTAR run.sh; GPL-3.0. See THIRD_PARTY_NOTICES.md.
if [ "${1:-}" = "--help" ]; then
    printf '%s\n' 'Usage: bash iSTAR_run.sh WORK_DIR (run from the upstream iSTAR repository)'
    exit 0
fi
SRBENCH_PYTHON="${SRBENCH_PYTHON:-python}"

: "${1:?Usage: bash iSTAR_run.sh work_dir}"
prefix="${1%/}/"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
"$SRBENCH_PYTHON" "$script_dir/install_patches.py" --method iSTAR --method-dir "$PWD" --verify

# Optional PID output for an external runtime/memory monitor.
record_pid() {
    if [ -n "${PID_FILE:-}" ]; then
        printf '%s\n' "$1" >> "$PID_FILE"
    fi
}

device="cuda"  # "cuda" or "cpu"
pixel_size=0.5  # desired pixel size for the whole analysis
n_genes=1000  # number of most variable genes to impute


# preprocess histology image
printf '%s\n' "$pixel_size" > "${prefix}pixel-size.txt"
"$SRBENCH_PYTHON" rescale.py "${prefix}" --image &
record_pid "$!"
wait "$!"

"$SRBENCH_PYTHON" preprocess.py "${prefix}" --image &
record_pid "$!"
wait "$!"

# extract histology features
"$SRBENCH_PYTHON" extract_features.py "${prefix}" --device=${device} &
record_pid "$!"
wait "$!"

# auto detect tissue mask
"$SRBENCH_PYTHON" get_mask.py "${prefix}embeddings-hist.pickle" "${prefix}mask-small.png" &
record_pid "$!"
wait "$!"


# Preserve a supplied gene list; select genes only when none is present.
if [ ! -f "${prefix}gene-names.txt" ]; then
    "$SRBENCH_PYTHON" select_genes.py --n-top=${n_genes} "${prefix}cnts.tsv" "${prefix}gene-names.txt" &
    record_pid "$!"
    wait "$!"
fi

# predict super-resolution gene expression
# rescale coordinates and spot radius
"$SRBENCH_PYTHON" rescale.py "${prefix}" --locs --radius &
record_pid "$!"
wait "$!"

# train gene expression prediction model and predict at super-resolution
"$SRBENCH_PYTHON" impute.py "${prefix}" --epochs=400 --device=${device} &
record_pid "$!"  # train model from scratch
wait "$!"

# visualize imputed gene expression
