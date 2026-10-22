#!/usr/bin/env bash
set -euo pipefail
if [ "${1:-}" = "--help" ]; then
    printf '%s\n' 'Usage: bash XFuse_run.sh CONFIG.toml OUTPUT_DIR'
    exit 0
fi
: "${1:?Supply a study XFuse config}"
: "${2:?Supply an output directory}"
xfuse run "$1" --save-path "$2"
