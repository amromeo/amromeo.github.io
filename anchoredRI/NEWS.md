# anchoredRI 0.3.0

- The point fit and bootstrap now default to Task E Arm B. Original-scale BIC selects among one to five components across eligible normal, lognormal and shifted-lognormal families. Set `rounding_unit` to the known assay increment and `fit_seed` to reproduce the multistart search.
- On RIbench before the safeguard, Arm B had the lowest score tested at N ≤ 5,000 (0.2224 vs 0.3696 for refineR). At N = 50,000, refineR scored 0.2282 and Arm B scored 0.3744. Use refineR or anchorRefine for datasets this large. Arm B had one implausible fit in 1,440 and selected five components in 941 (65.3%).

- Added default-on Step B LRT and observed-minimum checks for point and bootstrap fits, with `safeguard = FALSE` for research use.
- Withhold flagged intervals, retain anchoredRI and anchors-only estimates as diagnostics, record both reasons when applicable, and give the screened-patient recommendations. Add the indirect-method suggestion for at least 50,000 routine results.
- Distinguish bootstrap declines from numerical failures and document the accepted-draw confidence-interval policy in the README.
- Added the one-group and one-sided reporting support needed to reproduce the archived RIbench safeguard results, with reporting checks against archived fits.
- Added exact B2 confirmation and RIbench decline-count regression checks, and document the blind spot in heavily mixed populations with many mildly abnormal patients.

The checks are enabled by default. Declined estimates remain diagnostic only.

# anchoredRI 0.2.0

- Added `transform = "shifted_log"`, which jointly estimates a shared offset and reference/disease mixture parameters.
- Added `transform = "auto"`, which compares identity, log and shifted-log families with one to three components using BIC on the original measurement scale. The existing point-fit default remains identity.
- Added `bootstrap_anchored_ri()`, including resampling of both routine services and screened observations, full model reselection, offset refitting, reproducible parallel draws, failure records and boundary flags.
- Preserved the original four worked examples and comparator estimates. Added the updated model, 200 bootstrap draws per example and a separate paired old/new evaluation on 220 fresh synthetic datasets.
- Added tests for transformation likelihoods, analytic gradients, unit conversion and bootstrap reproducibility. Clinical performance and confidence-interval coverage remain unvalidated.
