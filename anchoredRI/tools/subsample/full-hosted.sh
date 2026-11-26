#!/usr/bin/env bash
# Restore the verified pilot runtime once per block. Each R process uses the
# same container, compiler build and single-threaded settings as the pilot.
set -euo pipefail
study=anchoredRI/tools/subsample
image_digest=sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3
image=rocker/r-ver@$image_digest
export TZ=UTC LC_ALL=C.UTF-8
if [[ $1 == setup || $1 == smoke ]]; then
    bundle="$study/results/pilot-45854822638/runtime"
    (cd "$bundle" && sha256sum -c runtime.sha256)
    tar xzf "$bundle/runtime.tgz"
    docker pull --platform linux/amd64 "$image"
    if [[ $1 == setup ]]; then
        python3 "$study/full-inputs.py" inputs --block "$2"
    else
        python3 "$study/prepare-inputs.py" inputs --ids 2065
        mkdir -p preflight/records/2065
        tar xzf "$study/results/pilot-45854822638/pilot/data/2065.tar.gz" -C preflight/records/2065
        python3 "$study/full-worker-check.py" seal preflight/records/2065 2065
        python3 "$study/full-worker-check.py" verify preflight/records/2065 2065
        bash "$study/full-hosted.sh" r "$study/full-worker.R" --verify preflight/records/2065 2065
        bash "$study/full-hosted.sh" r "$study/full-verify.R" preflight/records inputs preflight/verified pilot 2065
    fi
elif [[ $1 == r ]]; then
    shift
    exec docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" \
      -v "$PWD:/work" -w /work -e HOME=/tmp -e TZ -e LC_ALL -e GITHUB_JOB -e RUNNER_NAME \
      -e GITHUB_RUN_ID -e GITHUB_RUN_ATTEMPT -e TASK_H_CONTAINER_DIGEST="$image_digest" \
      -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -e OPENBLAS_CORETYPE=CORE2 \
      -e GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX512F,-AVX2,-AVX,-FMA,-FMA4 \
      -e R_LIBS=/work/runtime/package-lib:/work/runtime/deps "$image" Rscript --vanilla "$@"
else
    exit 2
fi
