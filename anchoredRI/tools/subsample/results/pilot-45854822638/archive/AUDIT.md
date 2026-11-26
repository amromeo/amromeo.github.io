# Full-data cap-3 agreement with the archived Task E baseline

The study-seed fit passed every required check on **1,354 of 1,440** matched datasets. The archived seed resolved **37 of 86** original exceptions.

The study seed was 1. Each diagnostic rerun used the original RIbench Index as its fitting seed. Both runs used the same archived routine values, anchors, rounding and requested percentiles. [Task E baseline](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/baseline-fits.csv)

## Study seed

| Check | Comparable | Pass | Fail | Unavailable | Inapplicable | Max abs. difference | Max fraction of width |
|---|---:|---:|---:|---:|---:|---:|---:|
| transform | 1440 | 1437 | 3 | 0 | 0 | — | — |
| components | 1440 | 1440 | 0 | 0 | 0 | — | — |
| lrt_flag | 1440 | 1429 | 11 | 0 | 0 | — | — |
| lower_flag | 1296 | 1293 | 3 | 0 | 144 | — | — |
| status | 1440 | 1429 | 11 | 0 | 0 | — | — |
| raw_lower_missingness | 1440 | 1440 | 0 | 0 | 0 | — | — |
| raw_upper_missingness | 1440 | 1440 | 0 | 0 | 0 | — | — |
| lower_missingness | 1440 | 1429 | 11 | 0 | 0 | — | — |
| upper_missingness | 1440 | 1429 | 11 | 0 | 0 | — | — |
| raw_lower | 1296 | 1261 | 35 | 0 | 144 | 17.44429 | 0.5169185 |
| raw_upper | 1440 | 1387 | 53 | 0 | 0 | 536.8502 | 227.7947 |
| lower | 1051 | 1032 | 19 | 11 | 378 | 17.44429 | 0.2236868 |
| upper | 1062 | 1036 | 26 | 11 | 367 | 22.00223 | 0.1953546 |
| p | 1440 | 1369 | 71 | 0 | 0 | 0.9749348 | — |

Endpoint differences are in analyte units. For two-sided fits the denominator is the archived raw interval width. For CRP it is the archived upper limit minus zero. The endpoint tolerance is 0.001 of that span. The absolute LRT p-value tolerance is 0.0001. Lower-endpoint missingness is also checked for CRP. Two matching declines have no released endpoints to compare.

Component 1's mean study-minus-archive share difference was 0.03952756 percentage points across 1,440 matched shares. The maximum absolute difference was 96.74578 percentage points. Share differences are descriptive.

## Archived-seed diagnostic

| Check | Comparable | Pass | Fail | Unavailable | Inapplicable | Max abs. difference | Max fraction of width |
|---|---:|---:|---:|---:|---:|---:|---:|
| transform | 86 | 85 | 1 | 0 | 0 | — | — |
| components | 86 | 86 | 0 | 0 | 0 | — | — |
| lrt_flag | 86 | 81 | 5 | 0 | 0 | — | — |
| lower_flag | 82 | 80 | 2 | 0 | 4 | — | — |
| status | 86 | 81 | 5 | 0 | 0 | — | — |
| raw_lower_missingness | 86 | 86 | 0 | 0 | 0 | — | — |
| raw_upper_missingness | 86 | 86 | 0 | 0 | 0 | — | — |
| lower_missingness | 86 | 81 | 5 | 0 | 0 | — | — |
| upper_missingness | 86 | 81 | 5 | 0 | 0 | — | — |
| raw_lower | 82 | 73 | 9 | 0 | 4 | 17.44504 | 0.2236964 |
| raw_upper | 86 | 66 | 20 | 0 | 0 | 536.8488 | 227.7947 |
| lower | 61 | 58 | 3 | 5 | 20 | 17.44504 | 0.2236964 |
| upper | 61 | 56 | 5 | 5 | 20 | 14.66443 | 0.129288 |
| p | 86 | 47 | 39 | 0 | 0 | 0.5687476 | — |

Endpoint differences are in analyte units. For two-sided fits the denominator is the archived raw interval width. For CRP it is the archived upper limit minus zero. The endpoint tolerance is 0.001 of that span. The absolute LRT p-value tolerance is 0.0001. Lower-endpoint missingness is also checked for CRP. Two matching declines have no released endpoints to compare.

Component 1's mean study-minus-archive share difference was 1.040547 percentage points across 86 matched shares. The maximum absolute difference was 96.74533 percentage points. Share differences are descriptive.

## Exceptions by original Index

| Index | Analyte | Study-seed differences | Archived-seed differences |
|---:|---|---|---|
| 427 | Hb | raw_lower, raw_upper | resolved |
| 865 | Ca | p | p |
| 866 | Ca | p | resolved |
| 881 | Ca | p | p |
| 882 | Ca | p | p |
| 891 | Ca | raw_lower, raw_upper | resolved |
| 898 | Ca | p | p |
| 914 | Ca | p | p |
| 929 | Ca | p | p |
| 930 | Ca | p | p |
| 946 | Ca | p | p |
| 961 | Ca | p | p |
| 977 | Ca | p | p |
| 993 | Ca | p | p |
| 994 | Ca | p | resolved |
| 1444 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1458 | FT4 | p | p |
| 1460 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1476 | FT4 | p | resolved |
| 1480 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1493 | FT4 | transform, raw_upper, upper, p | resolved |
| 1494 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1505 | FT4 | p | p |
| 1506 | FT4 | p | p |
| 1507 | FT4 | p | p |
| 1526 | FT4 | transform, raw_lower, raw_upper, lower, upper, p | resolved |
| 1537 | FT4 | p | p |
| 1542 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1554 | FT4 | p | p |
| 1556 | FT4 | raw_upper, upper, p | p |
| 1558 | FT4 | raw_lower, raw_upper, lower, upper, p | resolved |
| 1576 | FT4 | raw_lower, raw_upper, p | resolved |
| 2147 | AST | raw_lower, raw_upper | raw_lower, raw_upper |
| 2611 | LACT | raw_lower, raw_upper, lower, upper, p | p |
| 2621 | LACT | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p |
| 2623 | LACT | raw_lower, raw_upper | resolved |
| 2626 | LACT | raw_upper, upper, p | resolved |
| 2645 | LACT | lrt_flag, lower_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 2677 | LACT | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 2687 | LACT | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 2709 | LACT | raw_lower, raw_upper, lower, upper, p | resolved |
| 3219 | GGT | raw_lower, raw_upper, lower, upper, p | resolved |
| 3223 | GGT | raw_lower, raw_upper, lower, upper, p | raw_lower, raw_upper, lower, upper, p |
| 3237 | GGT | raw_lower, raw_upper, lower, upper, p | raw_lower, raw_upper, lower, upper, p |
| 3267 | GGT | raw_lower, raw_upper, lower, upper, p | resolved |
| 3268 | GGT | p | resolved |
| 3270 | GGT | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 3271 | GGT | raw_lower, raw_upper, lower, upper, p | resolved |
| 3284 | GGT | raw_lower, raw_upper, lower, upper, p | resolved |
| 3766 | TSH | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p |
| 3768 | TSH | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p |
| 3864 | TSH | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 4322 | IgE | p | p |
| 4323 | IgE | raw_upper, upper, p | resolved |
| 4325 | IgE | raw_upper | resolved |
| 4338 | IgE | raw_upper, upper, p | p |
| 4339 | IgE | lrt_flag, lower_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | lrt_flag, lower_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p |
| 4340 | IgE | raw_upper | raw_upper |
| 4353 | IgE | p | p |
| 4354 | IgE | p | p |
| 4355 | IgE | lrt_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | resolved |
| 4360 | IgE | raw_lower, raw_upper, lower, upper, p | resolved |
| 4362 | IgE | raw_upper | raw_upper |
| 4366 | IgE | raw_lower, raw_upper, lower, upper, p | resolved |
| 4369 | IgE | p | p |
| 4370 | IgE | p | p |
| 4371 | IgE | lrt_flag, lower_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p | lrt_flag, lower_flag, status, raw_lower, raw_upper, lower_missingness, upper_missingness, p |
| 4373 | IgE | raw_upper | raw_upper |
| 4378 | IgE | raw_lower, raw_upper, lower, upper, p | resolved |
| 4385 | IgE | raw_upper | raw_upper |
| 4386 | IgE | p | resolved |
| 4394 | IgE | raw_upper | raw_upper |
| 4404 | IgE | raw_upper, upper, p | raw_upper, upper, p |
| 4433 | IgE | p | p |
| 4434 | IgE | p | p |
| 4436 | IgE | raw_upper | resolved |
| 4438 | IgE | raw_upper, upper, p | raw_upper, upper, p |
| 4449 | IgE | p | p |
| 4450 | IgE | p | resolved |
| 4943 | CRP | raw_upper | raw_upper |
| 4989 | CRP | raw_upper | raw_upper |
| 5026 | CRP | raw_upper | raw_upper |
| 5031 | CRP | raw_upper | raw_upper |
| 5509 | LDH | p | p |
| 5539 | LDH | p | resolved |
| 5561 | LDH | transform, raw_lower, raw_upper, lower, upper, p | transform, raw_lower, raw_upper, lower, upper, p |

The files retain the Task E and study values on the same row, the two input hashes, fit seed, first-stage flags and decisions, component-1 shares and fit times. [All original comparisons](archive-audit.csv.gz) · [Original exceptions](exceptions.csv) · [Archived-seed reruns](seed-repeats.csv) · [Check summaries](audit-checks.csv)
