# RIbench N = 50000

RIbench evaluateAlgorithmResults stopped at a sample-size subset plot. Native per-dataset metrics and category scores were recovered using its unchanged computePerfMeasAll and getBenchmarkResults functions. The archived estimates were used for numerical scoring.

All 1,440 datasets were evaluated with each of three methods. This is a subset evaluation, not the published full RIbench benchmark score. AnchoredRI and anchorsOnly receive simulated screened anchors; refineR receives only routine data.

Runtime below is the sum of per-fit elapsed seconds, excluding generation, installation and workflow queueing. Parallel execution wall time should be reported separately.

| Subset | Method | N | Failed | Mean absolute Z deviation | Capped mean | Fit runtime (minutes) |
|---|---|---:|---:|---:|---:|---:|
| all | anchoredRI | 1440 | 0 | 0.9572 | 0.3656 | 40.25 |
| lambda=0.432 | anchoredRI | 144 | 0 | 0.5154 | 0.2575 | 3.86 |
| all | anchorsOnly | 1440 | 0 | 0.2961 | 0.2961 | 0.47 |
| lambda=0.432 | anchorsOnly | 144 | 0 | 0.2456 | 0.2456 | 0.04 |
| all | refineR | 1440 | 0 | 0.3099 | 0.2282 | 93.61 |
| lambda=0.432 | refineR | 144 | 0 | 0.1792 | 0.1792 | 8.62 |

The uncapped metric summarizes available results; inspect failure counts with it. Capped deviations and RIbench subset scores retain RIbench conventions. FT4 is the lambda=0.432 subset. CRP uses its one-sided 95th percentile.

Stage sample size is specified above. Runtime is the sum of fit times. Source versions and native outputs are retained.

## Generation and elapsed runtime

The run completed all 60 fitting batches and scoring in approximately 25 minutes, from 9:40 p.m. to 10:05 p.m. Eastern on November 22, 2026.

Generation took 155.21 minutes summed across 60 jobs. Median generation time per job was 151.39 seconds, with a range of 107.57 to 197.32 seconds. Each job ran RIbench's generator for its 24-dataset subset, including the generator's random-sequence advancement through the full design. These overlapping job times should not be interpreted as serial generation time for the whole stage.

The fit runtime table above sums per-dataset elapsed time for each method. It excludes generation, dependency installation and workflow queueing.

## Diagnostics

The diagnostic CSV preserves the small-N column names and adds LRL, URL, routine_min and lower_flag. The lower flag is true when the saved lower limit is strictly below the smallest routine value. CRP uses an upper-only 95th percentile, so its lower flag is false. Routine vectors were regenerated through RIbench 1.0.2 and all 1,440 hashes checked against its definitions to recover the minima. The saved estimates supply the limits and z deviations.

LRT p-values and statistics apply to anchoredRI. Missing LRT entries for the other methods are NA. RefineR uses a Box-Cox transform specified by Est_Lambda and shift, with one modeled reference component. RefineR bound-hit information is unavailable.

| Method | Lower-limit flags |
|---|---:|
| anchoredRI | 21/1440 |
| anchorsOnly | 0/1440 |
| refineR | 0/1440 |

Results use independent anchors seeded with 17000000 + Index. Source run and full artifact are available at https://github.com/amromeo/amromeo.github.io/actions/runs/45101564181 .
