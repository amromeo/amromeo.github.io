#!/usr/bin/env bash
set -euo pipefail
hold="$(dirname "${BASH_SOURCE[0]}")/reliability-v2/ACCURACY-HOLD.txt"
if [[ -f "$hold" ]]; then
  cat "$hold" >&2
  exit 1
fi
# Run from anchoredRI/. The Rscript executable may be supplied through RSCRIPT.
config=${1:?Supply the study configuration RDS}
rscript=${RSCRIPT:-Rscript}
"$rscript" tools/boxcox/test-models.R
"$rscript" tools/boxcox/run-study.R "$config" validate
"$rscript" tools/boxcox/run-study.R "$config" fit
"$rscript" tools/boxcox/score-study.R "$config"
