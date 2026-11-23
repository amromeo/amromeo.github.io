# Task E reliability review after round 3

**Arm B passed all three reproducibility requirements on the 156 hemoglobin datasets.** The original baseline discrepancy at Hb 142 arose because the GitHub fit missed a better solution for the three-component log model. Arm A stopped under its prespecified weak-identification rule (12/156, 7.7%). The improved baseline and Arm B are eligible for the accuracy comparison after review. [Protocol](../../PROTOCOL.md) · [Arm A stopping result](../../reliability-v3/results/REPORT.md)

## Why the baseline differed at Hb 142

Both runtimes used the same 1,000 routine results and 40 screened anchors. Each selected its lowest-BIC model. The identity model reached the same log-likelihood to nine decimal places. For the log model, the local fit reached a log-likelihood 3.905 units higher than the original GitHub fit. That difference changed which model won.

| Three-component model | Local log-likelihood | GitHub log-likelihood | Local BIC | GitHub BIC |
| --- | ---: | ---: | ---: | ---: |
| Identity | −2331.101 | −2331.101 | 4717.778 | **4717.778** |
| Log | **−2329.540** | −2333.445 | **4714.655** | 4722.465 |
| Shifted log | −2329.429 | −2333.141 | 4721.380 | 4728.805 |

The local runtime selected log and estimated an interval of **11.920 to 17.389**. GitHub selected identity and estimated **11.038 to 17.327**. The shifted-log model incurred an additional BIC penalty for its shift parameter. [All nine original candidates per runtime](hb142-original-candidates.csv)

Each runtime then restarted the identity and log candidates from both sets of saved parameters. The eight restarts converged with code 0. Every transferred starting likelihood reproduced exactly. For a given start, both runtimes returned identical limits and likelihoods at stored precision.

| Model and source of starting parameters | Final log-likelihood on either runtime | Lower limit | Upper limit |
| --- | ---: | ---: | ---: |
| Identity, local | −2331.101 | 11.038148 | 17.326609 |
| Identity, GitHub | −2331.101 | 11.038150 | 17.326610 |
| Log, local | −2329.540 | 11.919906 | 17.389050 |
| Log, GitHub | −2333.445 | 11.829069 | 17.473382 |

With access to both starts, both runtimes selected **three-component log, 11.919906 to 17.389050**. The better observed solution beat identity by 3.123 BIC points. The cross-start was diagnostic. The comparator file uses the original local baseline fit. [All cross-starts](hb142-cross-starts.csv) · [Optimizer codes](hb142-optimizer-blocks.csv) · [Selections after cross-starts](hb142-best-witnessed.csv)

## Arm B reproducibility

Both runtimes independently refit Index 1 through 156 from the routine and anchor data. Each evaluated 15 candidates, covering identity, log and shifted log with one through five components. The fixed round 3 fitter enforced the same original-scale SD floor as the improved baseline.

| Prespecified requirement | Result | Decision |
| --- | --- | --- |
| Same transform and component count in at least 152/156 | 155/156 (99.4%) | Passed |
| Both limits within 0.5% of the GitHub interval width in at least 149/156 | 155/156 (99.4%) | Passed |
| Zero nesting violations in each runtime | Local 0/2,322 eligible checks. GitHub 0/2,322. | Passed |

**Hb 150 was the sole exception to both agreement criteria.** The local fit selected five-component identity (12.0056 to 15.9738). GitHub selected four-component shifted log (12.1508 to 16.0904). The larger endpoint difference was 0.1452, or 3.69% of the GitHub interval width. The other 155 datasets all agreed within 0.012% of width. [All 156 pairs](arm-b-pairs.csv) · [Exception](arm-b-exceptions.csv)

All 156 pairs had finite limits. Each runtime had 330 structurally excluded nesting comparisons and zero unavailable eligible comparisons. All 780 zero-shift tests per runtime gave identical log-likelihoods and SD floors. Each runtime exactly reproduced its earlier likelihoods and endpoints for all nine one- through three-component baseline candidates on every dataset. [Gate counts](arm-b-gate.csv) · [Nesting checks by dataset](arm-b-checks-by-dataset.csv) · [Baseline repeats](baseline-repeatability.csv)

The routine and anchor hashes matched for every pair. All 17 source and harness hashes matched. GitHub reproduced the archived C2 image, R 4.5.0 and numerical libraries. The replay used the frozen fitter at commit ec3dd35. [Source hashes](source-hashes.csv) · [Local environment](local-runtime.txt) · [GitHub environment](Actions-runtime.txt) · [Completed Actions run](https://github.com/amromeo/amromeo.github.io/actions/runs/45185210923)

## Baseline estimates available to anchorRefine

The comparator CSV contains **540 unique local improved-baseline results** already fitted during reliability work. There are 396 small-N RIbench rows and 144 rows at N = 50,000. Each row includes the dataset key, routine and anchor hashes, selected model, fitted limits, safeguard decision and reported limits. The 156 overlapping rows match the dedicated baseline checkpoints exactly. The full 5,220-dataset export, including the simulation grid, awaits the accuracy study. [Comparator CSV](../../anchorRefine-comparator/anchoredRI-comparator.csv) · [Coverage and column definitions](../../anchorRefine-comparator/README.md)
