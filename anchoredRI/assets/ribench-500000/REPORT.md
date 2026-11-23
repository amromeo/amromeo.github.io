# RIbench results for 500,000 routine results

All 60 fitting batches and final scoring finished. Each method was evaluated on the same 1,440 independent RIbench test sets. AnchoredRI and anchorsOnly each received the same 40 screened draws per test set, generated with seed 17000000 + Index. Invalid draws were removed before rounding. All 1,440 anchor sample hashes were distinct. RefineR used the routine data alone.

RefineR had the lowest mean absolute z across all 1,440 datasets and in the FT4 subset. AnchoredRI had the smaller error on 480/1,440 pairs against refineR and 772/1,440 against anchorsOnly. The extreme-error counts below show the contribution of large misses.

## Performance

Lower standardized error is better. RIbench excludes datasets with mean absolute z above 5 from its native cutoff mean. The uncapped mean and the extreme-error counts retain those datasets. FT4 accounts for the 144 test sets with lambda = 0.432. CRP is scored on the upper limit only.

| Dataset group | Method | Scored / attempted | Mean absolute z | Native cutoff mean (n) | Total fit time (minutes) |
|---|---|---:|---:|---:|---:|
| All test sets | anchoredRI | 1440/1440 | 1.6267 | 0.4434 (1385) | 60.22 |
| All test sets | anchorsOnly | 1440/1440 | 0.2903 | 0.2903 (1440) | 0.48 |
| All test sets | refineR | 1440/1440 | 0.2027 | 0.2027 (1440) | 80.31 |
| FT4, lambda = 0.432 | anchoredRI | 144/144 | 1.5456 | 0.3684 (135) | 5.65 |
| FT4, lambda = 0.432 | anchorsOnly | 144/144 | 0.2440 | 0.2440 (144) | 0.04 |
| FT4, lambda = 0.432 | refineR | 144/144 | 0.1362 | 0.1362 (144) | 7.56 |

## Runtime

GitHub elapsed time from the original launch through the completed export was 28 min 49 s, from 11:09:25 PM to 11:38:14 PM Eastern. Generation took 154.26 minutes summed over parallel batches (median per batch 150.3 seconds, range 107.3 to 187.9). The method-specific times above sum the elapsed time of individual fits. All 60 fitting jobs completed 23 min 47 s after launch. Native scoring finished 25 min 46 s after launch. The remaining elapsed time includes fixing and rerunning the report export.

## Paired performance

Each difference is anchoredRI mean absolute z minus the comparator value for the same test set. Negative values favor anchoredRI. Values in the two wins columns count strictly smaller errors.

| Distribution | Pathological fraction | Comparator | Complete pairs | Mean difference | AnchoredRI wins | Ties | Comparator wins |
|---|---|---|---:|---:|---:|---:|---:|
| normal | ≤20% | anchorsOnly | 216/216 | 4.5848 | 110 | 0 | 106 |
| normal | >20% | anchorsOnly | 216/216 | 0.5824 | 57 | 0 | 159 |
| skewed | ≤20% | anchorsOnly | 216/216 | -0.1424 | 183 | 0 | 33 |
| skewed | >20% | anchorsOnly | 216/216 | 0.3987 | 120 | 0 | 96 |
| heavilySkewed | ≤20% | anchorsOnly | 216/216 | 1.3609 | 140 | 0 | 76 |
| heavilySkewed | >20% | anchorsOnly | 216/216 | 0.0094 | 119 | 0 | 97 |
| shifted | ≤20% | anchorsOnly | 72/72 | 1.5210 | 16 | 0 | 56 |
| shifted | >20% | anchorsOnly | 72/72 | 4.8260 | 27 | 0 | 45 |
| normal | ≤20% | refineR | 216/216 | 4.7922 | 20 | 0 | 196 |
| normal | >20% | refineR | 216/216 | 0.5993 | 28 | 0 | 188 |
| skewed | ≤20% | refineR | 216/216 | 0.0709 | 82 | 0 | 134 |
| skewed | >20% | refineR | 216/216 | 0.2841 | 102 | 0 | 114 |
| heavilySkewed | ≤20% | refineR | 216/216 | 1.6746 | 48 | 0 | 168 |
| heavilySkewed | >20% | refineR | 216/216 | 0.1114 | 112 | 0 | 104 |
| shifted | ≤20% | refineR | 72/72 | 1.3576 | 38 | 0 | 34 |
| shifted | >20% | refineR | 72/72 | 4.5256 | 50 | 0 | 22 |

## Extreme errors and diagnostic flags

| Method | Mean absolute z > 5 | Either eligible endpoint with absolute z > 5 |
|---|---:|---:|
| anchoredRI | 55 | 64 |
| anchorsOnly | 0 | 0 |
| refineR | 0 | 0 |

The LRT and lower-limit checks were recorded as diagnostics. All fitted intervals entered scoring. The lower-limit flag compares the lower estimate with the smallest routine value and is skipped for CRP.

| AnchoredRI group | Datasets | Lower flag | LRT p < 0.01 | Either flag | Both flags |
|---|---:|---:|---:|---:|---:|
| All anchoredRI fits | 1440 | 12 | 473 | 473 | 12 |
| Mean absolute z > 5 | 55 | 12 | 55 | 55 | 12 |
| Mean absolute z < 0.25 | 823 | 0 | 126 | 126 | 0 |
| Either eligible endpoint absolute z > 5 | 64 | 12 | 64 | 64 | 12 |

## Scoring and files

All 4,320 diagnostic rows were checked against the native per-dataset tables. The summary scores, paired wins and losses, error lists, independent seeds and lower-limit flags were independently recalculated from those rows. The tables preserve selected transforms, component counts, available bound indicators, LRT statistics, p-values and both standardized endpoint errors.

The RIbench sample-size subset plot stopped after computing its inputs. The scoring job used the package’s native per-dataset and category scoring functions to recover both tables from the original fitted intervals. Categories absent at this sample size have missing scores.

The 5,760 rounded anchor samples were distinct across all three RIbench stages.

The original job completed native scoring and saved all fitted objects. Its endpoint export stopped when refineR had no extreme rows. The score-only recovery handled the empty list and exported the diagnostics from the saved results. [Completed scoring archive](https://github.com/amromeo/amromeo.github.io/actions/runs/45112864524)

[Summary](summary.csv) · [RIbench category scores](stage-category-score.csv) · [Paired comparisons](paired-strata.csv) · [All endpoint extremes](extreme-errors.csv) · [Per-dataset diagnostics](per-dataset-diagnostics.csv) · [Diagnostic flag counts](diagnostic-flags.csv)

[Original fitting and native scoring run](https://github.com/amromeo/amromeo.github.io/actions/runs/45110114996)
