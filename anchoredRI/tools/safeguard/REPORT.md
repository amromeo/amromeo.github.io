# Verification of Arm B with the Task F safeguard

The proposed package selects up to five components with the Task E fitter. After the point fit and every bootstrap fit, Task F checks for disagreement with the screened sample (likelihood-ratio p < 0.01) and for a two-sided lower limit below the smallest routine result.

## Replay of Task E decisions

All 5,220 archived Arm B fits gave the same reporting decision in the installed package as in Task E. Each dataset was matched by ID. The replay verified the hashes of the routine and anchor vectors and used the archived reference-component mean, SD, transform and offset to calculate the LRT and the requested limits. It checked both flags, the one-sided skip, technical withholding, the reasons, and whether the package released or withheld the limits. The LRT was recalculated for all 432 one-sided CRP fits.

| Task E datasets | Fits | LRT flags | Lower-limit flags | Both flags | Declined | Reported |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RIbench, N = 1,000 or 5,000 | 2,880 | 99 | 0 | 0 | 99 | 2,781 |
| Simulation grid, replicates 1–20 | 900 | 91 | 0 | 0 | 91 | 809 |
| RIbench, N = 50,000 | 1,440 | 305 | 1 | 1 | 305 | 1,135 |

The comparison uses the published Task E results at commit `8bcd49c`. It replays the selected fits without rerunning the optimizer. [Small-N and grid archive](https://github.com/amromeo/amromeo.github.io/tree/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/results) · [Large-N archive](https://github.com/amromeo/amromeo.github.io/tree/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results)

## Refitting against the Task E source

On seven existing datasets, the public `fit_anchored_ri()` and the original Task E fitter agreed on the same R 4.5.0 runtime. Both chose the same transform and component count and returned the same safeguard flags. All candidate BICs, the selected weights and offset, diagnostic endpoints, and LRT statistic agreed to an absolute tolerance of 1e-10. The test used Task E's original multistart seed for each input and let the package infer the measurement increment.

Six of the seven local decisions also matched the archived GitHub Actions decisions. The exception was calcium index 870. The archived fit selected five normal components and triggered both flags. Both local implementations selected four shifted-log components and reported the interval. Applying the package safeguard to the archived calcium fit reproduced its decline. A fresh numerical fit on a different runtime can therefore change which interval reaches the safeguard and whether it is reported. The seven datasets are a software check, not an estimate of the frequency of this difference.

## Original Task F regressions

The original Task F branch passed its package tests and its three comparisons against frozen main on this runtime. After the rebase, the package also reproduced the historical Task F archive decisions below. Those decisions came from the earlier fitter. Task F used the Step B LRT, while B2's separate rule using a separation statistic declined 122/900.

| Historical fits | Fits | LRT flags | Lower-limit flags | Declined |
| --- | ---: | ---: | ---: | ---: |
| B2 confirmation, replicates 21–40 | 900 | 91 | 0 | 91 |
| Independent-anchor RIbench, N ≤ 5,000 | 2,880 | 533 | 9 | 533 |
| Independent-anchor RIbench, N = 50,000 | 1,440 | 402 | 21 | 402 |

Among B2 confirmation datasets, the LRT flagged 86/312 with error above 10% and 5/453 with error below 5%. On historical small-N RIbench it flagged all 23 fits with mean absolute z deviation above 5, 224/1,952 below 0.25 and 32/33 with either endpoint more than 5 z units from truth. For CRP in these older archives, only the upper limit and LRT statistic were available. The historical replay used the saved statistic for CRP and recalculated it from central limits and original anchors for the remaining fits.

## Package tests and reproduction

`R CMD check --no-manual` runs the model, offset and safeguard tests. The rebased tests cover the five-component candidates, gradients, one- and two-setting fits, one-sided limits, both decline reasons, technical withholding, the research override, and the large-sample recommendation. Bootstrap tests cover withheld point estimates, per-draw decisions and diagnostics, conditional confidence intervals, accepted-draw counts, and reproducibility with sequential and parallel workers.

Install the proposed package and run the following from its `anchoredRI/` directory. Set `task_e` to the archived Task E package source (the `anchoredRI/` directory), `small_data` to the RIbench directory containing its `Data/` tree, and `large_inputs` to the 24 archived Task E input shards. The archive tests use RIbench 1.0.2, data.table and digest. The historical replay uses the same configuration fields as the original Task F script: `b2`, `small`, `large`, `small_diagnostics`, `large_diagnostics`, and `output`.

```sh
Rscript tools/safeguard/replay-task-e.R "$task_e" "$small_data" "$large_inputs"
(cd "$task_e" && Rscript tools/boxcox/reliability-v3/build.R "$study_so")
Rscript tools/safeguard/test-parity.R "$task_e" "$small_data" "$study_so" "$large_inputs"
Rscript tools/safeguard/replay-archives.R "$task_f_config"
```

Use an absolute output filename for `study_so` and an existing Task F configuration RDS for `task_f_config`.
