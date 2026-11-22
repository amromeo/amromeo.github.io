# anchoredRI 0.2.0

- Added `transform = "shifted_log"`, which jointly estimates a shared offset and reference/disease mixture parameters.
- Added `transform = "auto"`, which compares identity, log and shifted-log families with one to three components using BIC on the original measurement scale. The existing point-fit default remains identity.
- Added `bootstrap_anchored_ri()`, including resampling of both routine services and screened observations, full model reselection, offset refitting, reproducible parallel draws, failure records and boundary flags.
- Preserved the original four worked examples and comparator estimates. Added the updated model, 200 bootstrap draws per example and a separate paired old/new evaluation on 220 fresh synthetic datasets.
- Added tests for transformation likelihoods, analytic gradients, unit conversion and bootstrap reproducibility. Clinical performance and confidence-interval coverage remain unvalidated.
