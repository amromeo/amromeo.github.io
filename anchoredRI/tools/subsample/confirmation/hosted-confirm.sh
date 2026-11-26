#!/usr/bin/env bash
# Use the compiled, SHA-256-verified Task H runtime inside its original container.
set -euo pipefail
study=anchoredRI/tools/subsample
confirm="$study/confirmation"
case "$1" in
    setup)
        runtime="$study/results/pilot-45854822638/runtime"
        (cd "$runtime" && sha256sum --quiet -c runtime.sha256)
        tar xzf "$runtime/runtime.tgz"
        docker pull --platform linux/amd64 rocker/r-ver@sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3
        ;;
    fixture)
        host="$2"
        mkdir -p "fixture-output"
        python3 - <<'PY'
import json
import os
import platform
from pathlib import Path
data = dict(boot_id=Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
            hostname=platform.node(), runner=os.environ["RUNNER_NAME"], image=os.getenv("ImageVersion"),
            run=os.environ["GITHUB_RUN_ID"], attempt=os.environ["GITHUB_RUN_ATTEMPT"])
Path("fixture-output/host.json").write_text(json.dumps(data, indent=2) + "\n")
PY
        bash "$study/full-hosted.sh" r "$confirm/fixtures.R" fixture-output/first
        if [[ $host == primary ]]; then
            bash "$study/full-hosted.sh" r "$confirm/fixtures.R" fixture-output/repeat
            bash "$study/full-hosted.sh" r "$confirm/compare-fixtures.R" same fixture-output/first fixture-output/repeat | tee fixture-output/same.txt
        fi
        ;;
    case)
        j="$2"
        root="blockdata/$j"
        bash "$study/full-hosted.sh" r "$confirm/worker.R" inputs "$root/records" "$j"
        bash "$study/full-hosted.sh" r "$confirm/verify-records.R" inputs "$root/records" "$root/verification" "$j"
        python3 "$confirm/schedule.py" seal "$j"
        ;;
    *) exit 2 ;;
esac
