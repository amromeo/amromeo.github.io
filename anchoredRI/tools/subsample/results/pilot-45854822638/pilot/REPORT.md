# Task H timing pilot

The 20 prespecified datasets span ten analytes and both pathology strata. Each used both banks of 50 draws. These paired scores use the pilot cohort for RIbench’s analyte-wide handling of zero limits. The approved full-cohort scorer requires all 1,440 datasets.

| Full first stage | Maximum components | Eligible | Majority releases | W releases | Common triples | Mean paired error change, median | Mean paired error change, W | Extreme intervals, reference / median / W |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Released | 5 | 14 | 14 | 14 | 14 | -0.023 | -0.007 | 0 / 0 / 0 |
| Declined | 5 | 6 | 5 | 4 | 4 | 0.026 | 0.063 | 0 / 0 / 0 |
| Released | 3 | 14 | 14 | 14 | 14 | 0.058 | 0.057 | 0 / 0 / 0 |
| Declined | 3 | 6 | 3 | 4 | 3 | -0.237 | -0.239 | 0 / 0 / 0 |

Seed set A uses every usable draw. Panel A uses the standard full first-stage limits. Panel B requires both a released median and an actual released full-package fallback. 0 full first stages were technically withheld. W’s own check determines the W release count. Both paired differences use the same scorable triples. An extreme interval has at least one requested endpoint with |z| > 5.

The paired tables include set B, passing-draw medians, all-eligible scores with nonreturns scored 5, W’s standalone releases, endpoint-specific extreme counts and changes in analyte units. [Per-Index pairs](paired-index.csv) · [Paired summaries](paired-summary.csv) · [Changes in analyte units](paired-limits-by-analyte.csv)

## Measured phase times

| Phase | Indices timed | Runner seconds, total | Runner seconds per Index, range | CPU seconds, total | CPU seconds per Index, range |
|---|---:|---:|---:|---:|---:|
| All timed per-Index work | 20 | 5630.6 | 77.3–857.0 | 5630.2 | 77.3–856.9 |
| W | 20 | 83.3 | 1.2–13.3 | 83.3 | 1.2–13.3 |
| aggregation | 20 | 2.8 | 0.1–0.2 | 2.8 | 0.1–0.2 |
| bank_A | 20 | 2709.8 | 37.1–414.4 | 2709.6 | 37.0–414.4 |
| bank_B | 20 | 2705.9 | 37.1–414.2 | 2705.7 | 37.1–414.2 |
| full_cap3 | 20 | 17.0 | 0.3–1.9 | 17.0 | 0.3–1.9 |
| full_public | 20 | 105.8 | 1.2–12.7 | 105.8 | 1.2–12.7 |
| output | 20 | 4.1 | 0.2–0.3 | 4.1 | 0.2–0.3 |
| preflight | 20 | 1.3 | 0.1–0.1 | 1.3 | 0.0–0.1 |
| fallback (included above) | 6 | 22.5 | 0.0–5.1 | 22.5 | 0.0–5.1 |
| sample (included above) | 20 | 1.6 | 0.1–0.1 | 1.6 | 0.1–0.1 |
| subsample_fit (included above) | 20 | 5405.0 | 73.7–827.9 | 5404.6 | 73.7–827.8 |

The standard full-data package called the fallback on 6 of the 20 datasets. Each stratum contributes a single observed time, so infrequent fallback calls, slow optimizations and differences in runner load contribute uncertainty.

## Compute projection

The total CPU column is a lower bound. Seven interrupted jobs used additional unmeasured CPU. The remaining projection uses the 20 completed pilot jobs. [Interrupted job accounting](ACCOUNTING.md)

| Work | Jobs | Occupied runner-hours | CPU-hours |
|---|---:|---:|---:|
| Measured archive | 48 | 0.662 | 0.344 |
| Measured archivecheck | 1 | 0.002 | 0.000 |
| Measured build | 1 | 0.008 | 0.001 |
| Measured failed build (CPU lower bound) | 1 | 0.007 | 0.001 |
| Measured failed pilot (CPU lower bound) | 6 | 0.053 | 0.000 |
| Measured fixturecheck | 1 | 0.002 | 0.000 |
| Measured fixtures | 2 | 0.018 | 0.005 |
| Measured pilot | 20 | 1.701 | 1.567 |
| Measured recovery preflight | 1 | 0.009 | 0.000 |
| Measured report | 1 | 0.024 | 0.000 |
| Projected remaining per-Index work | 1420 | 111.049 | 111.041 |
| Projected per-job overhead | 1420 | 9.691 | 0.228 |
| Projected final reporting | 1 | 0.007 | 0.000 |
| **Remaining** | | **120.746** | **111.269** |
| **Total** | | **123.233** | **113.187** |

The finished pilot-report job occupied its runner for 86.0 seconds and used 1.30 measured CPU seconds. It is included once in the spent balance. The future report allowance uses 25.0 runner seconds and 1.30 measured CPU seconds. The 61.0 seconds spent waiting for the original pilot run appear only in spent time. This replaces the preliminary allowance of 120 runner seconds and 30 CPU seconds. The full cohort will require larger transfers and scoring tables. The remaining report allowance assumes equivalent runtime and is an estimate.
Estimated remaining elapsed time is approximately **6.12 hours** with 20 concurrent runners, one dataset per job and 1,420 jobs in Index order, plus one final reporting job. This assumes immediate dispatch whenever a runner frees up. Queue delays and future failed-job retries would increase elapsed time.

Each remaining dataset receives the time of its pilot stratum (71 per stratum). Measured mean setup, export and other job overhead adds 24.6 runner seconds and 0.6 measured CPU seconds per job. The 20 completed pilot records are reused. The projection retains the timed standalone cap-3 fit in future per-Index work and counts prior archive verification and archived-seed reruns only in the measured total. Shared candidate fits for the two caps are timed once in each bank. CPU covers timed study processes. Runner time includes setup and artifact transfer. [Job overhead](job-overhead.csv) · [Projection inputs](compute-projection.csv)

## Pilot output checks

Readback verified draw counts, input hashes, the 50 votes per bank and cap, selected component bounds, and the coordinatewise raw medians against the exported draws.

Six passing-only medians had zero contributors. Both seed sets had zero passing fits for CRP 4968 under cap 3, CRP 5017 under cap 5, and LDH 5510 under cap 3. All 80 primary medians used all 50 draws.

Index, cap, bank and reason are recorded for each occurrence. [Pilot diagnostics](pilot-diagnostics.csv)
