#!/usr/bin/env bash
set -euo pipefail
# Run from anchoredRI/. The Rscript executable may be supplied through RSCRIPT.
config=${1:?Supply the study configuration RDS}
rscript=${RSCRIPT:-Rscript}
"$rscript" tools/boxcox/test-models.R
"$rscript" tools/boxcox/run-study.R "$config" validate
"$rscript" tools/boxcox/run-study.R "$config" fit
"$rscript" tools/boxcox/score-study.R "$config"
