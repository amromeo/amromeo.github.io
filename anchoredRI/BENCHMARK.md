# Comparing anchoredRI with refineR as intended

The primary experiment uses the authors' practical workflow: a prepared physiological partition, enough routine observations, bootstrap estimates and diagnostic plots. The earlier clinic/ward experiment is retained on `stress.html`.

## What changed

| Item | Earlier experiment | Primary comparison |
| --- | --- | --- |
| Routine cohort | 1,000 or 4,000 results; default mixture 61.25% affected | 5,000 results in each fit; 10% or 25% expected affected |
| Filtering | Outpatient arm had 250 or 1,000 observations | Both methods receive exactly the same eligible routine patients |
| Partitions | Single age and assay assumed | Single age/sex partition, assay and stable period assumed; one result per person |
| refineR | Default BoxCox; no bootstrap | 200 bootstrap resamples; medianBS; 95% confidence intervals |
| Model family | Normal reference distribution only | Normal, lognormal and shifted lognormal; modBoxCox for the prespecified shifted profile |
| Screened-only 120/400 | Mean ± 1.96 SD | Empirical 2.5th/97.5th percentiles, R quantile type 7 |
| Evaluation | Average larger endpoint error | Endpoint errors and bias, mean percentage error, reference-population coverage, CI coverage, failures |

## Sources and their implications

1. [Ammer et al. 2021](https://doi.org/10.1038/s41598-021-95301-2), [PubMed](https://pubmed.ncbi.nlm.nih.gov/34362961/). The original method assumes a majority of non-pathological results and a concentration region with negligible pathological contamination. Its simulations cover several distribution shapes and include deliberately adverse cases. The direct benchmarks use 120 and 400 draws from the reference distribution and estimate its percentiles. The published evaluation uses 100 seeds for indirect methods and 10,000 direct-method draws.
2. [Ammer et al. 2023 practical guide](https://doi.org/10.1093/jalm/jfac101). Select an appropriate analyte and compatible assay data; clean invalid results; use one result per patient; partition by age and sex where needed; exclude or stratify high-risk departments as appropriate. The guide reports better performance with at least 5,000 observations and pathological fractions at or below 30%. These are favorable conditions, not guarantees or an absolute 5,000-result minimum. Start with BoxCox; consider modBoxCox for a skewed distribution shifted from zero. Inspect plots and clinical plausibility.
3. [refineR 2.0.0 vignette](https://cloud.r-project.org/web/packages/refineR/vignettes/refineR_package.html) and [manual](https://cran.r-project.org/web/packages/refineR/refman/refineR.html). Preparation and partitioning occur before findRI. The vignette recommends at least 200 bootstrap resamples and medianBS point estimates. Bootstrap confidence limits are CILow/CIHigh; the runner disables the separate margins with the implementation’s `UCMargins = FALSE` argument. The separate uncertainty margins added in 2.0.0 are not substituted for bootstrap CIs here.

The prior pooled default has 25% clinic patients with 20% affected results and 75% ward patients with 75% affected results: 0.25 × 0.20 + 0.75 × 0.75 = 61.25%. This violates the majority-non-pathological assumption. The original authors themselves examined assumption violations; showing them is useful, but their results need a separate interpretation.

## Prespecified design

Each replicate contains 2,500 patients from each of two eligible services. Their reference distributions are identical; their disease fractions differ. Both methods receive those same 5,000 routine measurements. anchoredRI additionally receives 40 independent screened reference measurements and the two service labels. No simulated disease labels or true parameters enter either fitter. Routine patients and screened patients do not overlap. Screening does not use the analyte value. The 40-, 120- and 400-person groups are nested for paired comparisons.

| Profile | Reference distribution | Disease distribution | Fractions in services A/B | refineR model | anchoredRI transform |
| --- | --- | --- | --- | --- | --- |
| Normal, separated | Normal(30, SD 3) | Normal(42, SD 4) | 5% / 15% | BoxCox | Identity |
| Normal, overlap | Normal(30, SD 3) | Normal(36, SD 4.5) | 20% / 30% | BoxCox | Identity |
| Right-skewed | Lognormal(log 20, SDlog 0.35) | Lognormal(log 60, SDlog 0.30) | 5% / 15% | BoxCox | Log |
| Right-skewed, shifted | 10 + Lognormal(log 15, SDlog 0.50) | 10 + Lognormal(log 65, SDlog 0.25) | 5% / 15% | modBoxCox | Log |

All values use synthetic units and are rounded to 0.01. Truth is defined before rounding. The shifted profile deliberately falls outside anchoredRI's current identity/log Gaussian family. Neither method's model is selected by looking at error against the known reference limits. The shifted profile uses modBoxCox by design; its initial default BoxCox plot is also supplied. In clinical data the choice would follow inspection rather than a known simulation profile.

refineR uses its released findRI and getRI functions, including its median-model bootstrap point estimate. anchoredRI retains its original full-data point estimate and model selection; its percentile CIs resample within each service and resample the 40 screened patients independently, repeating model selection each time. Direct-method CIs use ordinary percentile bootstrap. Those CIs can be inaccurate in small tail samples; the exports record whether each CI contains its true endpoint. Bootstrap failures remain counted. A numerical boundary flag is not automatically a failed fit.

The published first run contains one worked dataset per profile, each with 200 bootstrap resamples. It shows fitted limits and uncertainty; it cannot estimate across-dataset bias, failure rates or CI coverage. The script also writes summary columns for larger runs, but single-dataset coverage indicators must not be read as performance rates. This is a demonstration of the practical workflow, not a replication of the original paper or proof of superiority. No method is tuned after examining results.

## Reading the metrics

- Endpoint bias is the mean signed difference from the known lower or upper limit. Endpoint MAE is the mean absolute difference.
- Mean percentage error averages the absolute relative errors at the two endpoints, multiplied by 100. Its Monte Carlo standard error describes variation across simulated datasets.
- Average larger limit error takes the larger absolute endpoint error within each dataset and then averages it.
- Reference-population coverage is the probability that a new reference individual's result falls inside the estimated interval, computed from the known generating distribution. The target is 95%. This differs from confidence-interval coverage, which records how often a 95% CI contains its true endpoint across datasets.
- Fits with errors are retained as failures. Performance summaries use successful fits and report the denominator. There is no favorable-error threshold, because these synthetic profiles have no prespecified analyte-specific total-error specification.

## Reproduce

The estimator implementation remains anchoredRI 0.1.0. This change replaces the evaluation workflow, not its fitting algorithm. `compare_ri_methods()` and `tools/build-demo.R` retain the earlier stress-test comparisons. Use the new script below for the prepared-cohort workflow.

Install refineR 2.0.0 and run from the package directory:

```sh
Rscript tools/build-guided-benchmark.R 8 1 200
# Replicated evaluation (longer run):
Rscript tools/build-guided-benchmark.R 8 100 200
```

The first argument selects parallel workers, the second independent datasets per profile, and the third bootstrap resamples (at least 200). Outputs include CSVs, software versions and the original package's fit plots. Cached cases are stored in `tools/guided-cache`; remove that cache before changing the generator or estimation code. `tools/build-demo.R` regenerates only the earlier stress experiment.

Clinical data would still need actual checks of assay stability, analyte-specific eligibility, physiological partitions and independent validation. Their absence cannot be repaired by a simulation.
