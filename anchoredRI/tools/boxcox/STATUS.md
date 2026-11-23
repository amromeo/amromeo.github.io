# Task E implementation status

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
