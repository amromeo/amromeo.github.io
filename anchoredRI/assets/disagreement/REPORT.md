# Step B: disagreement check

All 900 existing datasets were evaluated. Each of 60 distinct anchor samples was bootstrapped 1,000 times, with family and offset selection repeated. Shared bootstrap calculations were checked against identical anchors in all 15 corresponding scenarios. No thresholds were changed; replicates 21–40 were therefore not used.

Primary rule: report the original anchoredRI interval only if both limits fall within their respective marginal 95% bootstrap intervals. Otherwise decline and recommend at least 120 screened patients or verification of a suitable published interval with 20. Anchors-only estimates are diagnostic, never a fallback. This is a study evaluation; production package functions have not been changed.

## Primary check

Flagged 214/900 (23.8%); unavailable 0; declined 214/900 (23.8%).

Flagged 129/289 datasets with original error >10% (44.6%), and 64/485 with error <5% (13.2%).

Retained 686 intervals: mean error 9.08%, median 3.52%, 90th percentile 21.24%, maximum 96.74%; 23.3% still had error >10%.

### all: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 9/60 (15%) | 9/60 (15%) | 8/60 (13%) | 8/60 (13%) | 7/60 (12%) |
| moderate | 9/60 (15%) | 9/60 (15%) | 8/60 (13%) | 10/60 (17%) | 5/60 (8%) |
| heavy | 12/60 (20%) | 18/60 (30%) | 25/60 (42%) | 36/60 (60%) | 41/60 (68%) |

### normal: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 2/20 (10%) | 2/20 (10%) | 1/20 (5%) | 1/20 (5%) | 1/20 (5%) |
| moderate | 2/20 (10%) | 1/20 (5%) | 1/20 (5%) | 3/20 (15%) | 1/20 (5%) |
| heavy | 3/20 (15%) | 6/20 (30%) | 6/20 (30%) | 10/20 (50%) | 15/20 (75%) |

### lognormal: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) |
| moderate | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) | 4/20 (20%) |
| heavy | 4/20 (20%) | 7/20 (35%) | 13/20 (65%) | 13/20 (65%) | 12/20 (60%) |

### shifted: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 3/20 (15%) | 3/20 (15%) | 3/20 (15%) | 3/20 (15%) | 2/20 (10%) |
| moderate | 3/20 (15%) | 4/20 (20%) | 3/20 (15%) | 3/20 (15%) | 0/20 (0%) |
| heavy | 5/20 (25%) | 5/20 (25%) | 6/20 (30%) | 13/20 (65%) | 14/20 (70%) |

## Secondary check

Flagged 91/900 (10.1%); unavailable 0; declined 91/900 (10.1%).

Flagged 84/289 datasets with original error >10% (29.1%), and 7/485 with error <5% (1.4%).

Retained 809 intervals: mean error 9.47%, median 3.89%, 90th percentile 25.23%, maximum 86.36%; 25.3% still had error >10%.

Secondary declines and retained errors are hypothetical applications of the LRT alone. They do not change the primary output. The 2-df chi-squared reference is exploratory because the transformation was selected and anchors contributed to the original fit.

### all: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 1/60 (2%) | 1/60 (2%) | 1/60 (2%) | 1/60 (2%) | 1/60 (2%) |
| moderate | 1/60 (2%) | 1/60 (2%) | 0/60 (0%) | 0/60 (0%) | 0/60 (0%) |
| heavy | 1/60 (2%) | 1/60 (2%) | 9/60 (15%) | 26/60 (43%) | 47/60 (78%) |

### normal: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) |
| moderate | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) |
| heavy | 0/20 (0%) | 0/20 (0%) | 1/20 (5%) | 8/20 (40%) | 15/20 (75%) |

### lognormal: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 1/20 (5%) | 1/20 (5%) | 1/20 (5%) | 1/20 (5%) | 1/20 (5%) |
| moderate | 1/20 (5%) | 1/20 (5%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) |
| heavy | 1/20 (5%) | 1/20 (5%) | 4/20 (20%) | 8/20 (40%) | 13/20 (65%) |

### shifted: flag counts and rates

| Overlap | 10% affected | 30% | 50% | 70% | 90% |
|---|---:|---:|---:|---:|---:|
| separated | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) |
| moderate | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) | 0/20 (0%) |
| heavy | 0/20 (0%) | 0/20 (0%) | 4/20 (20%) | 10/20 (50%) | 19/20 (95%) |

## Interpretation limits

These percentages describe this simulation, with shared anchors across scenarios. The two marginal CIs do not imply a 5% overall flag rate. Performance among reported intervals is conditional on passing; it must be read together with the decline rate and the errors that remain undetected. No thresholds were tuned after inspection. The original anchoredRI error is measured against simulated truth; that truth is never used by the disagreement rule.

Full per-scenario counts, denominators and retained-error summaries are in summary.csv. datasets.csv separates diagnostic estimates from reported limits and includes recommendations. bootstrap-estimates.csv contains all 60,000 bootstrap fits.
