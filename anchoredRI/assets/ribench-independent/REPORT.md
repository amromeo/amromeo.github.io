# RIbench: N <= 5,000

RIbench evaluateAlgorithmResults stopped at a sample-size subset plot. Native per-dataset metrics and category scores were recovered using its unchanged computePerfMeasAll and getBenchmarkResults functions. No fits were rerun.

All 2,880 datasets were evaluated with each of three methods. This is a subset evaluation, not the published full RIbench benchmark score. AnchoredRI and anchorsOnly receive simulated screened anchors; refineR receives only routine data.

Runtime below is the sum of per-fit elapsed seconds, excluding generation, installation and workflow queueing. Parallel execution wall time should be reported separately.

| Subset | Method | N | Failed | Mean absolute Z deviation | Capped mean | Fit runtime (minutes) |
|---|---|---:|---:|---:|---:|---:|
| all | anchoredRI | 2880 | 0 | 0.4393 | 0.2754 | 64.58 |
| lambda=0.432 | anchoredRI | 288 | 0 | 0.2618 | 0.2402 | 6.15 |
| all | anchorsOnly | 2880 | 0 | 0.2899 | 0.2899 | 0.90 |
| lambda=0.432 | anchorsOnly | 288 | 0 | 0.2277 | 0.2277 | 0.09 |
| all | refineR | 2880 | 0 | 0.9553 | 0.3696 | 252.10 |
| lambda=0.432 | refineR | 288 | 0 | 0.2954 | 0.2954 | 25.26 |

The uncapped metric summarizes available results; inspect failure counts with it. Capped deviations and RIbench subset scores retain RIbench conventions. FT4 is the lambda=0.432 subset. CRP uses its one-sided 95th percentile.

No N=50,000 or N=500,000 fit was run. Report these runtime results before proceeding to those stages. All source versions, native Rdata outputs and per-analyte CSVs are retained.
