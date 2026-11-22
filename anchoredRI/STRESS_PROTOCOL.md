# Outside favorable conditions: first pass

Frozen on 22 November 2026, before inspecting study results. The anchoredRI 0.2 estimator is unchanged.

There are 45 scenarios: three reference shapes (normal, lognormal and shifted lognormal), three degrees of disease overlap and five affected fractions (10%, 30%, 50%, 70%, 90%). Each dataset contains 5,000 routine patients from two equally sized services, with affected fractions five percentage points below and above the overall fraction. There is one observation per person and a single stable age/sex and assay partition.

On the reference distribution's standardized normal scale, disease has mean 4, 2 or 0.75 and SD 1.25. The reference has mean 0 and SD 1. The same monotone transformation applies to both populations. This keeps the mixture within anchoredRI's model family; it does not test distributional misspecification. `assets/stress-grid/design.csv` records the fraction of disease values inside the true reference interval.

Twenty independently generated routine datasets per scenario test anchoredRI automatic selection and the empirical screened-only estimates from 40, 120 and 400 people. The screened samples are independent of routine patients, perfectly labeled and nested. For a given family and replicate, the same screened sample is reused across scenarios to hold anchor quality fixed. Different replicates have independent screened samples. No true parameters or affected labels enter a fitter.

Replicate 1 in every scenario also receives the full bootstrap comparison: released refineR 2.0.0 with 200 resamples and medianBS, anchoredRI with 200 full-pipeline resamples, and screened-only percentile estimates with 200 resamples. refineR uses BoxCox for normal/lognormal and modBoxCox for shifted lognormal. These 45 worked comparisons cannot establish repeated superiority or confidence-interval coverage rates. Repeated refineR comparisons remain a later extension.

Primary error is the larger absolute endpoint error, divided by the true interval width and multiplied by 100. A large error exceeds 10% of that width; this is an exploratory reporting threshold, not a clinical acceptance criterion. Mean percentage endpoint error, signed endpoint errors, both reference tail probabilities, numerical failures and boundary flags are also retained. Failed estimates must remain in attempted denominators. Mean errors among successful fits must be accompanied by failure counts.

The first pass holds routine sample size at 5,000. The proposed 1,000- and 500-patient extensions, imperfect anchors and incompatible distribution shapes are not included in this run.

Run from `anchoredRI/`:

```sh
Rscript tools/test-stress-design.R
Rscript tools/run-stress-grid.R both 0 4 20
```

The GitHub Actions workflow splits the 45 cases across standard public runners. It uploads outputs for review; it does not automatically publish results to the website. Source hashes reject changed estimators or study runners when using cached results.
