# Task E accuracy at small N and on the simulation grid

The improved three-component baseline and five-component Arm B were evaluated on 2,880 RIbench datasets and 900 grid datasets. Each fit used 40 screened anchors before any RIbench filtering. Arm A stopped during reliability testing.

The 50-dataset replay passed for 49/50 improved baseline fits and 49/50 Arm B fits. Two distinct datasets missed. The prespecified rule therefore required Actions fits for every dataset. Index 4298 changed from three to two shifted-log components for the baseline, and Index 3033 changed from four to five for Arm B. [Replay results](runtime-replay.csv)

All 120 fitting batches succeeded. The workflow took 59 minutes 50 seconds including setup, queues and the replay decision. The two Ubuntu 24.04 runner image revisions fitted 2,062 and 1,718 datasets, respectively (20261116.133.1 and 20261116.376.1). Both used R 4.5.0, RIbench 1.0.2 and the pinned fitting source. The frozen local baseline was fitted in the earlier local runtime. Numerical runtime effects can contribute to differences from this comparator. [Per-dataset provenance](runtime-provenance.csv) · [GitHub Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45257620744)

## Small-N RIbench

| Method | Before | After | Reported / attempted | Declined | Technical |
| --- | --- | --- | --- | --- | --- |
| Improved baseline | 0.271 | 0.187 | 2406/2880 | 474 (16.5%) | 0 |
| Arm B | 0.222 | 0.208 | 2781/2880 | 99 (3.4%) | 0 |
| Original frozen local fit | 0.277 | 0.192 | 2347/2880 | 533 (18.5%) | 0 |
| Archived C2 | 0.275 | 0.191 | 2347/2880 | 533 (18.5%) | 0 |
| Anchors only | 0.290 | 0.290 | 2880/2880 | 0 (0.0%) | 0 |
| refineR | 0.370 | 0.370 | 2880/2880 | 0 (0.0%) | 0 |

Values are native RIbench mean absolute z deviations after its cutoff. The native cutoff omits dataset errors above 5. The uncapped mean, number included by the cutoff, extreme endpoints and counts of implausible fits are reported separately. After-safeguard means describe the intervals retained by each method. The archived anchors-only and refineR methods use their original reporting behavior.

| Small-N criterion | Improved score | Frozen baseline threshold | Decision |
| --- | --- | --- | --- |
| Improved baseline | 0.271 | 0.277 | Passed |
| Arm B | 0.222 | 0.277 | Passed |

| Method | Filter | Scored | Uncapped mean absolute z | Mean absolute z > 5 |
| --- | --- | --- | --- | --- |
| Improved baseline | before | 2880 | 0.438 | 22 |
| Improved baseline | after | 2406 | 0.187 | 0 |
| Arm B | before | 2880 | 0.241 | 2 |
| Arm B | after | 2781 | 0.208 | 0 |
| Original frozen local fit | before | 2880 | 0.441 | 23 |
| Original frozen local fit | after | 2347 | 0.192 | 0 |
| Archived C2 | before | 2880 | 0.439 | 23 |
| Archived C2 | after | 2347 | 0.191 | 0 |
| Anchors only | before | 2880 | 0.290 | 0 |
| Anchors only | after | 2880 | 0.290 | 0 |
| refineR | before | 2880 | 0.955 | 32 |
| refineR | after | 2880 | 0.955 | 32 |

### Skewed sets with lambda = 0.432

| Method | Before | After | Reported / attempted | Declined | Technical |
| --- | --- | --- | --- | --- | --- |
| Improved baseline | 0.236 | 0.180 | 267/288 | 21 (7.3%) | 0 |
| Arm B | 0.185 | 0.185 | 288/288 | 0 (0.0%) | 0 |
| Original frozen local fit | 0.240 | 0.181 | 268/288 | 20 (6.9%) | 0 |
| Archived C2 | 0.240 | 0.181 | 268/288 | 20 (6.9%) | 0 |
| Anchors only | 0.228 | 0.228 | 288/288 | 0 (0.0%) | 0 |
| refineR | 0.295 | 0.295 | 288/288 | 0 (0.0%) | 0 |

## Simulation grid

| Method | Before | After | Reported / attempted | Declined | Technical |
| --- | --- | --- | --- | --- | --- |
| Improved baseline | 14.039 | 9.446 | 809/900 | 91 (10.1%) | 0 |
| Arm B | 14.039 | 9.446 | 809/900 | 91 (10.1%) | 0 |
| Original frozen local fit | 14.060 | 9.469 | 809/900 | 91 (10.1%) | 0 |
| Archived grid fit | 14.060 | 9.469 | 809/900 | 91 (10.1%) | 0 |
| Anchors only | 12.625 | 12.625 | 900/900 | 0 (0.0%) | 0 |
| refineR | 95.615 | 95.615 | 900/900 | 0 (0.0%) | 0 |

Grid errors are the larger absolute endpoint error as a percentage of the true reference interval width. Replicates 1 to 20 were used in development. The paired comparisons use the same routine values and anchors in every method.

## Change from the original fits

| Dataset group | Comparator | Filter | Common pairs | Mean change | Wins / ties / losses |
| --- | --- | --- | --- | --- | --- |
| Small N, uncapped mean absolute z | Original frozen local fit | before | 2880 | -0.003 | 1509 / 3 / 1368 |
| Small N, uncapped mean absolute z | Original frozen local fit | after | 2309 | -0.005 | 1221 / 3 / 1085 |
| Small N, uncapped mean absolute z | Archived C2 | before | 2880 | -0.001 | 1523 / 1 / 1356 |
| Small N, uncapped mean absolute z | Archived C2 | after | 2308 | -0.005 | 1225 / 1 / 1082 |
| Grid, width-error percentage points | Original frozen local fit | before | 900 | -0.021 | 455 / 2 / 443 |
| Grid, width-error percentage points | Original frozen local fit | after | 809 | -0.024 | 410 / 0 / 399 |
| Grid, width-error percentage points | Archived grid fit | before | 900 | -0.021 | 447 / 0 / 453 |
| Grid, width-error percentage points | Archived grid fit | after | 809 | -0.023 | 400 / 0 / 409 |

A negative paired change favors the improved baseline. Paired after-safeguard results include only datasets on which both methods reported an interval.

Improved baseline had 3 of 45 grid cells with mean error higher than the frozen local baseline by more than two paired Monte Carlo SEs.

Arm B had 3 of 45 grid cells with mean error higher than the frozen local baseline by more than two paired Monte Carlo SEs.

The distribution and pathological-fraction table reports paired differences, wins, ties, losses and the common reporting denominator before and after safeguarding. [RIbench paired comparisons](ribench-pairs.csv) · [Grid comparisons and MCSE](grid-pairs.csv)

## BIC near ties

| Stage | Method | Runner-up within 2 BIC | Median largest endpoint difference in near ties (% of interval width) |
| --- | --- | --- | --- |
| small | Improved baseline | 251/2880 (8.7%) | 4.28% |
| small | Arm B | 517/2880 (18.0%) | 3.85% |
| grid | Improved baseline | 44/900 (4.9%) | 18.31% |
| grid | Arm B | 44/900 (4.9%) | 18.31% |

One-sided CRP fits use the zero-to-upper reporting span. [Near-tie summary](runner-up-summary.csv)

## Winning optimizer exits

| Stage | Method | Iteration-limit wins | Line-search-error wins |
| --- | --- | --- | --- |
| small | Improved baseline | 9 | 11 |
| small | Arm B | 245 | 11 |
| grid | Improved baseline | 0 | 23 |
| grid | Arm B | 0 | 23 |

The exit log includes the source start of the winning candidate. A logged exit leaves the prespecified estimator and safeguard decision intact.

## Runtime and N = 50,000 estimate

| Stage | Completed | Fresh | Cached | Wall time (min) | Fresh fitting CPU (hours) |
| --- | --- | --- | --- | --- | --- |
| small | 2880 | 2880 | 0 | 56.6 | 2.44 |
| grid | 900 | 900 | 0 | 57.9 | 16.23 |
| combined | 3780 | 3780 | 0 | 58.2 | 18.67 |

The two stages ran concurrently on GitHub Actions with four workers per batch. Wall times span the first to last accepted fit and include queue gaps. Extraction and scoring took 4.6 additional minutes. The projected fitting time for the remaining 1,296 large datasets is 0.08 hours at observed throughput (0.08 hours using median archived timing ratios). Repeating all 1,440 from scratch would take approximately 0.09 hours. This uses fresh small-N fitting times and the earlier C2 50,000-to-5,000 timing ratio for each analyte. A cross-analyte correction supplies the 5,000 estimate where all corresponding small fits were cached. This is a planning estimate for fitting. Scoring and publication add time.

The comparator CSV now includes 2,880 small-N datasets, 900 grid datasets, and the 144 pre-existing large-N Hb datasets. The large-N accuracy comparison awaits approval.

[Comparator CSV](../../anchorRefine-comparator/anchoredRI-comparator.csv) · [Summary and paired tables](.)
