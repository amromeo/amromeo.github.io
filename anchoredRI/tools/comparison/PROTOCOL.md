# Repeated comparator study

Specified on 22 November 2026 after reviewing the first stress study and before examining the new comparator results. This is a follow-up using the existing simulated datasets, not an independent validation set.

The production estimator is frozen at commit f463a2cb84439ef9a9bb71e44f76661ceb0df00f. No files in R/ or the original stress design are changed. All 45 scenarios and 20 replicates per scenario are regenerated from their original seeds: 900 datasets, each with 5,000 routine results and an independent nested screened sample. No disease labels, true parameters or error scores enter mixture fitting or component selection.

## Methods

1. Released refineR 2.0.0, BoxCox for normal/lognormal scenarios and modBoxCox for shifted-lognormal scenarios, as in the initial study. NBootstrap=0, getRI(pointEst='fullDataEst', UCMargins=FALSE). This estimates accuracy of the full-data fit; it is distinct from the previous medianBS estimate. Bootstrap coverage and the effect of bootstrap aggregation require a subsequent analysis. refineR retains its original family-informed configuration, rather than forcing an automatic transformation procedure it does not share with the other methods.
2. Conventional finite Gaussian mixture on a selected transformation, pooling routine results without anchors. The component with greatest fitted weight is designated healthy. A second reported interval uses the component with the lowest interval midpoint; this sensitivity rule favors this simulation's upward disease shifts and is not a general clinical identification rule. Both rules are fixed before the full run and neither chooses the component by proximity to truth.
3. The same pooled mixture likelihood plus 40 known healthy observations, assigned to component 1. The contribution of these observations is an ordinary likelihood term, without extra weight.
4. Frozen anchoredRI with those same 40 anchors and the actual two service labels.
5. Empirical 2.5th and 97.5th percentiles from the nested 40, 120 and 400 screened observations.

The three mixture methods consider identity, log where valid, and shifted-log transformations and one, two or three components, selected by BIC with measurement-scale Jacobians. Bounds, multistart patterns, optimizer tolerances and iteration budgets are inherited from the frozen estimator. Without anchors, the first component's initial mean and SD come from the routine sample instead of anchors. All methods are deterministic given their inputs/seeds. These comparators are implemented as study-only ablations of the same likelihood, not as an independent third-party mixture package benchmark. Pooled weights use fewer parameters, with the BIC penalty adjusted accordingly.

Tests verify agreement between the generalized two-group likelihood and frozen anchoredRI, finite-difference agreement with the analytic shifted-log gradient for one/two groups and zero/known anchors, and successful pooled estimation. Source hashes are stored with every cached fit and reject reuse after source changes.

## Outcomes and interpretation

Report mean, median and 90th percentile of the larger absolute endpoint error normalized by true interval width; signed lower/upper errors, reference tail probabilities, fit failures, boundary hits and errors above the exploratory 10%-of-width threshold. Retain failures in attempted denominators. Report paired mean error differences and paired wins for each scenario, separating failed pairs from finite comparisons. Negative differences favor the first named method. Monte Carlo SE of a paired mean difference is the SD of the 20 paired differences divided by sqrt(20), not a confidence interval for a clinical population.

Primary contrasts: anchoredRI versus refineR; pooled anchored versus unanchored largest-component mixture; anchoredRI versus pooled anchored mixture. The lowest-component sensitivity and screened-only benchmarks are secondary. Aggregate scenario summaries give every scenario equal weight. Replicates share anchors across scenarios within a family, so 900 rows must not be treated as independent for a pooled inferential confidence interval. The same screened values across conditions make screened-only curves flat within each family.

The study does not establish novelty, clinical validity, behavior with misspecified distributions, imperfect anchors or realistic service heterogeneity. The simulated service proportions differ by only ten percentage points. Results will be reviewed before changing the estimator or updating claims on the public results page.
