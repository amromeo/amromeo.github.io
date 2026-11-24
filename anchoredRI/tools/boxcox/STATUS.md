# Task E implementation status

Task E closed on 23 November 2026. At N = 50,000, neither the improved three-component baseline (0.3618) nor Arm B (0.3744) met the prespecified before-safeguard RIbench score of 0.2282 set by refineR. Lower scores indicate better accuracy. The baseline had 37 implausible fits among 1,440 attempts. Arm B had 1 and passed the count criterion of at most 14. Arm B selected the five-component maximum in 941 of 1,440 fits (65.3%).

At N ≤ 5,000, Arm B had the lowest before-safeguard score among the tested methods (0.2224 across 2,880 datasets, compared with 0.3696 for refineR). At N = 50,000, use refineR or anchorRefine. Task E directly compared refineR. Further anchoredRI method development is on hold pending the project lead's decision.

[Small-N and grid results](accuracy/results/REPORT.md) · [Large-N results](accuracy/large/results/REPORT.md)

Arm B and the improved baseline tied on all 900 simulation-grid datasets. Arm A stopped during reliability testing. [Protocol](PROTOCOL.md) · [Comparator for anchorRefine](anchorRefine-comparator/README.md)

## Task G (planned, not started)

Recorded 24 November 2026. The proposed method returns the anchoredRI interval when its safeguard passes. If anchoredRI declines and there are at least 5,000 routine results, try anchorRefine (or refineR plus the anchor check). Report the fallback interval only if it passes the anchor check. Below 5,000 routine results, a decline stands.

On the RIbench discovery data, anchoredRI’s reported intervals had lower error than refineR on the same datasets at each of the three sizes. A separate discovery comparison used anchoredRI when it reported an interval and refineR when it declined. Lower scores indicate better accuracy.

| Routine results per dataset | anchoredRI then refineR | refineR alone |
|---:|---:|---:|
| 1,000 | 0.24 | 0.51 |
| 5,000 | 0.20 | 0.33 |
| 50,000 | 0.16 | 0.24 |

The 1,000-result row is part of the discovery comparison. The proposed Task G fallback starts at 5,000. Testing will wait for anchorRefine Stage 3 and will use its confirmation datasets.

## Earlier validation runs

Implementation started after user approval on 22 November 2026. No package files under `R/` were changed.

The analytic gradient checks, multiplicity checks, inverse-domain checks and initial frozen identity/log/automatic-fit checks passed. A synthetic fitted-shape smoke test also completed.

Full input and equivalence validation began and stopped in its first batch. All 4,320 independent anchor hashes matched the archived hashes. Routine-data hashes for the first six datasets matched. The study-local fitter agreed with the frozen source on those datasets, but the archived C2 endpoint check failed for indices 2, 3 and 5. Indices 1, 4 and 6 passed all checks.

For Hb Index 2, the archived limits are 11.68842156845287 and 16.22519510342150. The local reproduction gives 11.68819692220466 and 16.22547400917307. Both use three shifted-log components. The discrepancies exceed the predeclared tolerance and cannot be ignored.

The frozen C2 source hashes match the artifact's source manifest exactly. Local R and archived R are both version 4.5.0, but the local system is Debian 13 with reference BLAS and the archived run used Ubuntu with OpenBLAS. This is a possible numerical explanation, not yet a demonstrated cause.

No experimental benchmark fits or performance comparisons have run. The validation gate remains closed. No threshold, optimizer setting, estimator or scoring rule has been changed to accommodate this discrepancy. The next step is to reproduce the baseline in the archived runtime and resolve the discrepancy before continuing.

## Validation amendment approved

The user approved the cross-machine validation amendment on 22 November 2026. It is recorded in PROTOCOL.md. A fresh validation run preserves the original failed run and uses the unchanged estimators. Every same-machine comparison must pass absolute tolerance 1e-10. Cross-machine differences are audited against 0.1% of archived interval width, with a review gate for any model change or more than 1% outside tolerance. The baseline for scoring is now the frozen C2 source refitted in the experimental runtime. Experimental fits remain gated on completion of this validation.

## Amended review gate triggered

Validation paused on 22 November 2026 after 156 of 5,220 datasets. All completed same-machine checks passed at absolute tolerance 1e-10. Among these 156 RIbench datasets, two selected a different transform than the archived C2 result. Both also exceed the 0.1%-of-archived-width limit tolerance.

| Analyte | Index | N | Local transform | Archived transform | Components, local / archived | Largest limit difference / archived width |
|---|---:|---:|---|---|---|---:|
| Hb | 100 | 1,000 | log | identity | 3 / 3 | 4.170293% |
| Hb | 148 | 5,000 | identity | shifted_log | 3 / 3 | 4.207564% |

The largest absolute limit difference is 0.1707087 for Index 100. The largest width-normalized difference is 4.207564% for Index 148. These are partial validation findings, not estimates for the full benchmark. They trigger the user-approved pause because any selected model change requires review.

The controller has stopped after completing its active batch. `HALT.txt` remains in the run directory. No experimental benchmark fits or scoring have started. No estimator, fitting setting, outcome or success criterion was changed. The exception and summary tables are saved under `validation/`. Continuing requires the user's review of this gate; it must not be bypassed automatically.


## Review decision and local BICs

The full audit now runs through all 4,320 archived comparisons. It pauses if more than 216 datasets change model. The exact same-machine checks cover all 5,220 datasets and stop on any failure. The original 156 checkpoints did not retain the other families’ candidate limits. A fresh run records all candidate limits so each near tie can be evaluated. The earlier checkpoints are preserved.

| Hb index | Local winner | Next local candidate | BIC gap | Upper limit, winner → runner-up |
|---|---|---|---:|---:|
| 100 | log, 3 components | identity, 3 | 5.7511 | 16.1338 → 15.9278 |
| 148 | identity, 3 components | shifted-log, 2 | 37.2374 | 15.9583 → 16.2142 |

Neither local result is within 2 BIC points. For Index 148, the archived winner (shifted-log, 3 components) ranks third locally, 41.3002 points behind identity. The two CSVs in `validation/` hold all 18 BICs, candidate limits and the two margins. The original C2 percentile objects contain only percentiles, estimates and runtime. The original custom CSVs contain the winning transform and component count. Candidate BIC tables from the GitHub Actions reproduction will be identified as new calculations.
