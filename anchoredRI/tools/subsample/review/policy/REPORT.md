# Package-first fallback policy on RIbench

The revised policy returned 1,436 of 1,440 intervals under the primary subsample seeds. It kept all 1,423 package returns and added 13. On the common returns, mean endpoint error fell from 0.1476 to 0.1308. The sensitivity seeds gave 1,435 returns and a common-return mean of 0.1311.

Each result below uses the original 40 anchors and cap-5 fits from Task H at nominal N = 50,000. The subsample estimate is the coordinatewise median of raw endpoints from all usable draws among the 50. The full-data first-stage decision takes precedence. Step a returns the median when the full-data first stage passes. After a first-stage decline, step b returns the full package result when available. If the package withholds, step c returns the median only when the subsample majority approves. Step d otherwise withholds. The two seed sets change only the subsamples. These are retrospective comparisons on the datasets that informed the rule.

| Method | Returns | Nonreturns | Mean error on returns | All-assignment score | Returned intervals with an extreme endpoint |
|---|---:|---:|---:|---:|---:|
| Full package | 1,423 | 17 | 0.1476 | 0.2049 | 0 |
| Revised policy, A | 1,436 | 4 | 0.1316 | 0.1451 | 0 |
| Revised policy, B | 1,435 | 5 | 0.1317 | 0.1486 | 0 |

Each dataset contributes one error, the mean absolute RIbench z score of its requested endpoints. CRP contributes its upper endpoint alone. An extreme means at least one requested endpoint had absolute z above 5. The all-assignment score caps each returned dataset error at 5 and assigns 5 to a nonreturn. [Cohort summaries](summary.csv)

## Paired accuracy on the 1,423 common returns

| Seed set | Policy mean | Package mean | Mean paired difference, policy minus package | Policy better | Exact ties | Policy worse |
|---|---:|---:|---:|---:|---:|---:|
| A | 0.13077 | 0.14761 | −0.01684 | 572 | 297 | 554 |
| B | 0.13106 | 0.14761 | −0.01655 | 558 | 297 | 568 |

The 297 ties in each bank are the package fallbacks copied at step b. All paired changes therefore came from step a. Among its 1,126 datasets, the package mean was 0.13112, and the median estimates had mean errors of 0.10984 with A and 0.11020 with B. [Paired results](common-returns.csv)

## Routes and additional returns

| Step | Source | Datasets A / B | Mean error on returns A / B |
|---|---|---:|---:|
| a | Subsample median after full fit passed | 1,126 / 1,126 | 0.1098 / 0.1102 |
| b | Returned package fallback | 297 / 297 | 0.2101 / 0.2101 |
| c | Majority-approved median after package withheld | 13 / 12 | 0.2221 / 0.2102 |
| d | No interval | 4 / 5 | — |

Step c added 13 of 17 package nonreturns under A and 12 of 17 under B. Its maximum dataset error was 0.5642 under A and 0.5546 under B. Neither set had a step-c endpoint beyond absolute z of 5. Index 3871 was released under A and withheld under B. All other step assignments were identical across seeds. [Route summaries](routes.csv) · [Per-Index gates, routes, limits, endpoint errors and paired differences](by-index.csv)

The analysis reads the published Task H paired records and verifies their SHA-256 hash before routing. The input includes the full first-stage decision panel, package return status and returned reference endpoints, all-draw median endpoints and majority decision under each seed set. The command `python3 tabulate.py` in this directory regenerates all four CSV files using Python 3.11 or later. [Tabulation code](tabulate.py) · [Source records at the results commit](https://github.com/amromeo/amromeo.github.io/blob/0b525d24f71e972032649be538dd434288b49d7b/anchoredRI/tools/subsample/results/full-45860591866/paired-index.csv)

As an independent scoring check, we assembled each policy’s actual endpoints over all 1,440 assignments and ran the pinned RIbench 1.0.2 scorer again for both seed sets and the package. The upper estimates were positive for every returned interval in every analyte. All 4,320 dataset-method assignments reproduced the saved requested endpoint z scores, dataset errors and missing values. The command `python3 verify-scores.py /path/to/testsets.csv` repeats this check in the frozen Task H R environment. [Scoring readback by analyte and method](score-readback.csv) · [Verification code](verify-scores.py)
