#!/usr/bin/env bash
set -eu
# Run from anchoredRI. The fit controller owns run.lock for the full fitting stage.
rscript=$1
config=$2
out=$3
flock -s "$out/run.lock" true
test -f "$out/FITTED.txt"
before=$(date +%s)
"$rscript" tools/boxcox/accuracy/extract.R "$config" "$out"
"$rscript" tools/boxcox/accuracy/score.R "$config" "$out"
after=$(date +%s)
printf '%s\n' "$((after-before))" > "$out/scoring-seconds.txt"
small=$("$rscript" -e 'cat(readRDS(commandArgs(TRUE)[1])$small_archive)' "$config")
large=$("$rscript" -e 'cat(readRDS(commandArgs(TRUE)[1])$large_archive)' "$config")
python tools/boxcox/accuracy/finish.py "$out" "$small" "$large"
