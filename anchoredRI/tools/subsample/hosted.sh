#!/usr/bin/env bash
# Run from the checkout root on a GitHub-hosted VM. All R CPU timings are
# taken inside the container, where R and its child processes execute.
set -euo pipefail
export LC_ALL=C.UTF-8 TZ=UTC
study=anchoredRI/tools/subsample
image=rocker/r-ver@sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3
stage0=https://raw.githubusercontent.com/amromeo/anchorRefine/638e2d3500faa01d4e2816895a95d7a6685567a0/stage0/frozen
phase=$1
shift
mkdir -p output/account tmp-work
job_key="${TASK_H_ROLE:-${GITHUB_JOB:-local}}-${TASK_H_KEY:-single}"
timelog="output/account/$job_key.tsv"
measure() {
    local label=$1
    shift
    local TIMEFORMAT=$'%R\t%U\t%S'
    local receipt
    receipt=$(mktemp)
    { time "$@" 2>&3; } 3>&2 2>"$receipt"
    printf '%s\t%s\n' "$label" "$(cat "$receipt")" >> "$timelog"
    rm "$receipt"
}
container() {
    local label=$1
    shift
    docker run --rm --platform linux/amd64 --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work \
        -e HOME=/tmp -e TZ -e LC_ALL -e GITHUB_JOB -e RUNNER_NAME -e GITHUB_RUN_ID -e GITHUB_RUN_ATTEMPT \
        -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 -e OPENBLAS_CORETYPE=CORE2 \
        -e GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX512F,-AVX2,-AVX,-FMA,-FMA4 \
        -e R_LIBS=/work/runtime/package-lib:/work/runtime/deps "$image" \
        bash -ec 'TIMEFORMAT=$(printf "%%R\t%%U\t%%S"); log=$1; label=$2; shift 2
            exec 3>&2
            { time "$@" 2>&3; } 2>tmp-work/timing
            printf "%s\t%s\n" "$label" "$(cat tmp-work/timing)" >> "$log"' bash "$timelog" "$label" "$@"
}
prepare() {
    measure prepare python3 "$study/prepare-inputs.py" inputs "$@"
}
restore() {
    sha256sum -c runtime.sha256
    measure unpack tar xzf runtime.tgz
    docker pull --platform linux/amd64 "$image" >/dev/null
}
case "$phase" in
    build)
        python3 "$study/freeze.py" verify
        mkdir -p runtime/package-lib tmp-work/anchoredRI output/runtime
        curl -fsSL --retry 4 "$stage0/bundle.tar.gz" -o tmp-work/bundle.tgz
        curl -fsSL --retry 4 "$stage0/installed-files.sha256" -o tmp-work/upstream.sha256
        echo 'ab482e6a48109c476c8d4e6d14231e64f28470259708cdba728b7b8ec0e91027  tmp-work/bundle.tgz' | sha256sum -c -
        echo '00595b3a1a9f8ccba0ce89931d055020bd3d876c392323d7bbd459ad8a734e77  tmp-work/upstream.sha256' | sha256sum -c -
        tar xzf tmp-work/bundle.tgz -C runtime deps
        (cd runtime; awk '$2 ~ /^deps\//' ../tmp-work/upstream.sha256 | sha256sum -c --quiet -)
        cp -R anchoredRI/{DESCRIPTION,NAMESPACE,R,src,inst,man,LICENSE} tmp-work/anchoredRI/
        docker pull --platform linux/amd64 "$image"
        docker image inspect "$image" > output/runtime/image.json
        container build R CMD INSTALL --no-multiarch -l runtime/package-lib tmp-work/anchoredRI
        container versions Rscript --vanilla "$study/runtime.R" "$study/library-lock.csv" output/runtime
        measure pack tar czf runtime.tgz runtime
        sha256sum runtime.tgz > runtime.sha256
        cp runtime.sha256 output/runtime/
        ;;
    fixtures)
        restore
        for repetition in "$@"; do
            container "fixture-$repetition" Rscript --vanilla "$study/fixtures.R" "output/$repetition"
        done
        ;;
    audit)
        restore
        prepare --shards "$1"
        mkdir -p output/audit
        container archive Rscript --vanilla "$study/archive-audit.R" inputs output/audit "$1"
        ;;
    pilot)
        restore
        prepare --ids "$1"
        container pilot Rscript --vanilla "$study/pilot.R" inputs output/records "$1"
        ;;
    fixture-report)
        measure fixture-comparison python3 "$study/compare-fixtures.py" incoming/fixtures-primary/first incoming/fixtures-primary/second incoming/fixtures-other/other output/verification
        ;;
    audit-report)
        measure audit-report python3 "$study/report-audit.py" incoming output/verification
        ;;
    report-R)
        container score Rscript --vanilla "$@"
        ;;
    report)
        restore
        mkdir -p tmp-work/bin
        cat > tmp-work/bin/Rscript <<'SHIM'
#!/usr/bin/env bash
set -eu
args=()
for value in "$@"; do args+=("${value#"$PWD/"}"); done
exec bash anchoredRI/tools/subsample/hosted.sh report-R "${args[@]}"
SHIM
        chmod +x tmp-work/bin/Rscript
        export PATH="$PWD/tmp-work/bin:$PATH"
        prepare --metadata-only
        python3 "$study/job-account.py" jobs.json incoming output/jobs.csv
        measure pilot-report python3 "$study/pilot-report.py" all-records inputs output/pilot-report --roster "$study/pilot-roster.csv" --jobs output/jobs.csv --concurrency 20 --indices-per-job 1 --report-wall-seconds 120 --report-cpu-seconds 30
        ;;
    *) exit 2 ;;
esac
hostname > "output/account/$job_key.host"
cat /proc/sys/kernel/random/boot_id >> "output/account/$job_key.host"
printf '%s\n' "${RUNNER_NAME:-local}" "${ImageVersion:-local}" >> "output/account/$job_key.host"
