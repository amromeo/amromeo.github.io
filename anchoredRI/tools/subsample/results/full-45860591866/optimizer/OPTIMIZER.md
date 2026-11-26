# Selected optimizer exits and endpoint error

The analysis covers 1,440 datasets and 293,760 first-stage fit records. Each subsample contributes a record under each cap. The full and W fits enter once per cap. The denominator for an exit rate includes every attempted fit in that group. A failed selection appears as unavailable.

The frozen fitter labels code 0 as converged. It labels code 1 or an iteration-limit message as iteration_limit, a line-search or LNSRCH message as line_search_failure, and other codes as other_exit. The label belongs to the selected mixture candidate. Either check can decline a fit after any optimizer exit. [Status mapping in anchoredRI 0.4.0](https://github.com/amromeo/amromeo.github.io/blob/fe8105ac6961f1b45bbfbe034648938f227c4726/anchoredRI/R/mixture_fit.R)

## Before the safeguards

| Fit | Cap | Attempts | Converged | Iteration limit | Line-search failure | Other exit | Unavailable |
|---|---:|---:|---:|---:|---:|---:|---:|
| Full | 3 | 1,440 | 1,428 | 10 | 2 | 0 | 0 |
| Full | 5 | 1,440 | 957 | 481 | 2 | 0 | 0 |
| W | 3 | 1,440 | 1,430 | 8 | 2 | 0 | 0 |
| W | 5 | 1,440 | 1,220 | 218 | 2 | 0 | 0 |
| Subsamples A | 3 | 72,000 | 71,534 | 262 | 204 | 0 | 0 |
| Subsamples B | 3 | 72,000 | 71,582 | 223 | 195 | 0 | 0 |
| Subsamples A | 5 | 72,000 | 62,813 | 8,988 | 199 | 0 | 0 |
| Subsamples B | 5 | 72,000 | 62,666 | 9,142 | 192 | 0 | 0 |

Each cell below gives iteration-limit exits / line-search exits / fits selecting that K. The values use all attempted fits before the safeguards.

| Fit | Cap | K = 1 | K = 2 | K = 3 | K = 4 | K = 5 |
|---|---:|---:|---:|---:|---:|---:|
| Full | 3 | 0 / 0 / 0 | 0 / 2 / 80 | 10 / 0 / 1,360 | — | — |
| Full | 5 | 0 / 0 / 0 | 0 / 2 / 80 | 2 / 0 / 60 | 104 / 0 / 366 | 375 / 0 / 934 |
| W | 3 | 0 / 1 / 21 | 0 / 1 / 146 | 8 / 0 / 1,273 | — | — |
| W | 5 | 0 / 1 / 19 | 0 / 1 / 135 | 8 / 0 / 253 | 95 / 0 / 660 | 115 / 0 / 373 |
| Subsamples A | 3 | 0 / 64 / 1,019 | 0 / 139 / 7,311 | 262 / 1 / 63,670 | — | — |
| Subsamples B | 3 | 0 / 63 / 1,008 | 0 / 130 / 7,251 | 223 / 2 / 63,741 | — | — |
| Subsamples A | 5 | 0 / 62 / 955 | 0 / 137 / 6,860 | 150 / 0 / 12,203 | 3,303 / 0 / 32,353 | 5,535 / 0 / 19,629 |
| Subsamples B | 5 | 0 / 62 / 966 | 0 / 128 / 6,784 | 128 / 2 / 12,239 | 3,391 / 0 / 32,242 | 5,623 / 0 / 19,769 |

Signed lower and upper z scores measure direction. Absolute z measures endpoint error. The dataset score d is the uncapped mean absolute z across its requested endpoints. CRP uses only its upper endpoint. An extreme fit has at least one requested endpoint with absolute z above 5. Each metric has its own denominator.

| Fit | Cap | Exit | Scored fits | Mean d | Converged mean d | Difference | Fits with an extreme endpoint / assessed |
|---|---:|---|---:|---:|---:|---:|---:|
| Full | 3 | iteration_limit | 10 | 0.5868 | 0.8540 | -0.2671 | 0 / 10 |
| Full | 3 | line_search_failure | 2 | 0.0163 | 0.8540 | -0.8376 | 0 / 2 |
| Full | 5 | iteration_limit | 481 | 0.2941 | 0.4699 | -0.1758 | 5 / 481 |
| Full | 5 | line_search_failure | 2 | 0.0163 | 0.4699 | -0.4535 | 0 / 2 |
| W | 3 | iteration_limit | 8 | 0.1681 | 0.5156 | -0.3475 | 0 / 8 |
| W | 3 | line_search_failure | 2 | 0.0113 | 0.5156 | -0.5042 | 0 / 2 |
| W | 5 | iteration_limit | 218 | 0.3297 | 0.1750 | 0.1546 | 0 / 218 |
| W | 5 | line_search_failure | 2 | 0.0113 | 0.1750 | -0.1637 | 0 / 2 |
| Subsamples A | 3 | iteration_limit | 262 | 0.5962 | 0.4786 | 0.1176 | 4 / 262 |
| Subsamples A | 3 | line_search_failure | 204 | 0.0330 | 0.4786 | -0.4456 | 0 / 204 |
| Subsamples B | 3 | iteration_limit | 223 | 0.3102 | 0.5003 | -0.1901 | 1 / 223 |
| Subsamples B | 3 | line_search_failure | 195 | 0.0293 | 0.5003 | -0.4710 | 0 / 195 |
| Subsamples A | 5 | iteration_limit | 8,988 | 0.3092 | 0.2001 | 0.1091 | 3 / 8,988 |
| Subsamples A | 5 | line_search_failure | 199 | 0.0316 | 0.2001 | -0.1686 | 0 / 199 |
| Subsamples B | 5 | iteration_limit | 9,142 | 0.2978 | 0.1985 | 0.0992 | 2 / 9,142 |
| Subsamples B | 5 | line_search_failure | 192 | 0.0291 | 0.1985 | -0.1695 | 0 / 192 |

The detailed table gives counts, denominators, exit rates, signed and absolute endpoint errors, and differences from converged fits for each cap and selected component count. The released analysis uses the same first-stage scores on fits passing both applicable safeguards. [All exit and error summaries](optimizer-summary.csv)

## Subsamples compared within the same dataset

For each dataset with both exit types, we average each metric over its draws separately by exit and subtract the converged mean from the other exit mean. Each dataset with both means contributes once. The selected-K rows also hold the selected component count fixed. Full and W have one fit per cap per dataset.

| Seed set | Cap | Exit | Datasets with paired d | Mean within-dataset difference in d |
|---|---:|---|---:|---:|
| A | 3 | iteration_limit | 139 | -0.1362 |
| A | 3 | line_search_failure | 110 | -0.0197 |
| B | 3 | iteration_limit | 128 | -0.3183 |
| B | 3 | line_search_failure | 105 | -0.0522 |
| A | 5 | iteration_limit | 1,151 | 0.0063 |
| A | 5 | line_search_failure | 106 | -0.0101 |
| B | 5 | iteration_limit | 1,144 | -0.0023 |
| B | 5 | line_search_failure | 103 | -0.0190 |

These comparisons describe the saved fits. Optimizer exit and model selection both depend on the draw. [Within-dataset contrasts before and after the checks](optimizer-within-index.csv)
