# anchoredRI 0.4.0

- `fit_anchored_ri()` tries the checked anchorRefine Stage 3 post-fit blend when safeguarded anchoredRI declines or withholds an interval. The blend uses the Stage 2 peak-choice fit with weight 0.25. Set `fallback = FALSE` for the anchoredRI-only point fit.
- The result identifies the reporting method and retains the first decision and the diagnostics from both attempted methods. If neither reports an interval, the result gives the screened-patient advice.
- The fallback includes the refineR 2.0.0 point-fit functions and uses Stage 1's rounding-aware randomized PIT check. The external refineR package is optional for comparisons.
- Fallback accepts at least 40 screened patients and evaluates the exact Stage 1 null at the supplied count. Point and bootstrap fits warn at counts other than 40, the only count evaluated for both methods in Task G.
- Bootstrap uses the combined point fit. When the fallback reports, confidence limits are missing and `ci_reason` explains why. Confidence-interval coverage has yet to be validated for either method.
- The size message records that accuracy between 5,000 and 50,000 routine results has not been tested directly.

# anchoredRI 0.3.1

- Large-dataset messages now recommend refineR; anchorRefine is labeled as in development. No change to fitting or safeguard decisions.

# anchoredRI 0.3.0

- The point fit and bootstrap now default to Task E Arm B. Original-scale BIC selects among one to five components across eligible normal, lognormal and shifted-lognormal families.
- Both functions return the assay rounding unit. If you omit it, a message gives the coarsest shared step of the form {1, 2, 2.5, 5} × 10^k for integer k from −8 through 3, or zero if no step fits.
- Results include the winning optimizer’s exit status, code and message, and the runner-up model and its BIC gap. A reported fit that reaches the iteration limit gives a warning and retains its interval. The fit seed and computing environment can occasionally affect model selection. The README describes how the same seed selected different models for calcium dataset 870 in two environments.
- On RIbench before the safeguard, Arm B had the lowest score tested at N ≤ 5,000 (0.2224 vs 0.3696 for refineR). At N = 50,000, refineR scored 0.2282 and Arm B scored 0.3744. Arm B had one implausible fit in 1,440 and selected five components in 941 (65.3%). Every call with more than 5,000 routine results now identifies the untested range between 5,000 and 50,000 and suggests refineR or anchorRefine.

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
