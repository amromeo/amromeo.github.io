# anchoredRI

anchoredRI estimates a laboratory reference interval from routine results and results from about 40 chart-reviewed patients who meet the reference criteria. It first fits the anchoredRI model, which uses the screened patients to pick out the reference group. If that fit fails its checks, it tries anchorRefine, a refineR-based fallback adjusted toward the screened patients and checked against them. If neither passes, it returns no interval and suggests what to do next.

**Research software.** anchoredRI has been tested only on simulated data. Clinical performance on patient data has yet to be established.

**Fully AI-developed.** A pediatric laboratory physician posed the question and directed the exploration. OpenAI ChatGPT designed and wrote the method implementation, simulation experiments, tests, package documentation and website.

[Home](https://amromeo.github.io/anchoredRI/) · [Evidence](https://amromeo.github.io/anchoredRI/evidence.html) · [Package help](https://amromeo.github.io/anchoredRI/help.html) · [Earlier studies](https://amromeo.github.io/anchoredRI/earlier-studies.html) · [Source](https://github.com/amromeo/amromeo.github.io/tree/main/anchoredRI)

## What you need

- Routine results from one group that can share a reference interval, with the same analyte, assay, age group and sex.
- One result per patient.
- A care-setting label for each routine result, with one or two settings.
- Results from 40 patients who meet your reference criteria after chart review.
- The assay’s reporting increment, if you know it. Otherwise the package infers it.

The simulations drew the routine and screened results separately, representing different patients in each group. The tested routine dataset sizes were 1,000, 1,500, 5,000 and 50,000. [Study report](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md)

[Project roadmap](ROADMAP.md)

## Fitting and reporting a reference interval

`fit_anchored_ri()` compares normal, lognormal and shifted-lognormal mixtures with one to five components per family, using BIC on the measurement scale. It reports the selected fit’s limits when the fit passes its checks. The default is the central 95% interval. Use `probs = .95` for a one-sided upper limit.

| Outcome | `fit$method` | `fit$status` |
| --- | --- | --- |
| First fit supplies the interval | `"anchoredRI"` | `"reported"` |
| Fallback supplies the interval | `"anchorRefine"` | `"reported"` |
| No interval | `NA` | `"declined"` or `"withheld"` |

Read the limits in `fit$limits`. The first decision and its reason stay in `fit$anchoredRI_status` and `fit$anchoredRI_reason` even if the fallback reports. When no interval is available, `fit$recommendation` suggests obtaining at least 120 screened patients for a direct interval or verifying a suitable published interval with 20 screened patients.

The fallback pools the routine results for the Box-Cox fit. It uses the screened patients to help choose the reference peak and uses the original Box-Cox fit if peak selection fails. It adjusts the estimated reference mean and spread toward the screened results with a fixed weight of 0.25, then checks how well the adjusted distribution fits their results, accounting for rounding. The package includes the refineR 2.0.0 code it needs. The external refineR package is optional for the separate comparator. [Fallback research record](https://github.com/amromeo/anchorRefine)

The fallback requires at least 40 nonnegative screened results, at least 11 nonnegative routine results and a positive reporting increment. Its check accounts for the actual number of screened patients. The package warns at counts other than 40 because the combined-method study evaluated exactly 40. Set `fallback = FALSE` to use the first fit alone. The `fallback_seed` controls the fallback search and `pit_seed` controls the random draws used to check rounded results.

Supply the assay’s reporting increment with `rounding_unit` when known. It sets a minimum component spread in measurement units. If you omit it, the fitter finds the largest shared step of the form {1, 2, 2.5, 5} × 10^k, for integer k from −8 through 3, that fits all routine and screened values. The fit reports the unit in `rounding_unit` and sends a message when it infers the value. Zero means there is no shared step or that you supplied zero for continuous data. All bootstrap draws use the same unit.

`runner_up` identifies the alternative model with the next-lowest BIC. Its `bic_gap` is the runner-up BIC minus the winning BIC. The fit also returns the winning optimizer’s exit information, including `optimizer_status`. An interval that passes its checks can be reported with an iteration-limit warning. The optimizer’s starting points and the computing environment can occasionally change the selected model. RIbench calcium dataset 870 is an example. [Validation and reproduction](tools/safeguard/REPORT.md)

## Safeguard

The likelihood-ratio test checks whether the screened results agree with the fitted reference group. It compares the screened results under that fit with a model fitted to the screened results alone on the same transformed scale. A p-value below 0.01 causes a decline. A two-sided lower limit below every routine result also causes a decline. The lower-limit check is skipped for one-sided requests. If the test statistic cannot be calculated, the interval is withheld.

The first fit keeps its joint and screened-only diagnostics in `diagnostics$anchored` and `diagnostics$anchors_only`. The fallback’s diagnostics are in `diagnostics$anchorRefine`. For research use, `safeguard = FALSE` skips the checks and the fallback. When many patients have disease and their results overlap heavily with the reference group, an inaccurate interval can pass both checks.

`bootstrap_anchored_ri()` resamples routine results within settings and resamples the screened patients. It uses the first model for the bootstrap draws. If the fallback supplies the point estimate, its limits are available and the confidence limits are missing. Read `ci_reason` for the explanation. For a reported first fit, confidence intervals use the finite draws that passed the checks and require at least 80% of the requested draws. Coverage on either method has yet to be validated. [Package help](https://amromeo.github.io/anchoredRI/help.html#fallback-ci)

## Evidence

In the independent confirmation simulations, the combined method met all 16 prespecified conditions across 7,488 datasets. It returned an interval for 98–100% of datasets at every size tested. On descriptive RIbench, when the share of patients with disease was 20% or less, refineR had lower mean error at 5,000 and 50,000 results among the intervals it returned. Accuracy between 5,000 and 50,000 has yet to be tested directly. The package prints a message above 5,000 routine results. [Tables and limitations](https://amromeo.github.io/anchoredRI/evidence.html) · [Study record (Task G)](https://github.com/amromeo/amromeo.github.io/blob/655a007f6d5ec67d41262c45f21293b3d0b945b5/anchoredRI/tools/combined/results/REPORT.md)

## Install from GitHub

Use R 4.1.0 or later and a C++17 compiler. On Windows, install the Rtools version for your R release. On a Mac, use the Xcode Command Line Tools or Xcode. [Package requirements](https://github.com/amromeo/amromeo.github.io/blob/53d7b947bb40849b9b0b7556f97c676a0845928f/anchoredRI/DESCRIPTION) · [Rtools](https://cran.r-project.org/bin/windows/Rtools/) · [Mac instructions](https://developer.apple.com/documentation/xcode/installing-the-command-line-tools/)

```r
install.packages(c("remotes", "ash", "future.apply"))
remotes::install_github(
  "amromeo/amromeo.github.io",
  ref = "anchoredRI-v0.4.0", subdir = "anchoredRI"
)
```

Or download the source tarball and run the following from the folder where you saved it. [Download anchoredRI 0.4.0](https://amromeo.github.io/anchoredRI/downloads/anchoredRI_0.4.0.tar.gz) · [Previous versions](https://amromeo.github.io/anchoredRI/help.html#install)

```r
install.packages(c("ash", "future.apply"))
install.packages("anchoredRI_0.4.0.tar.gz", repos = NULL, type = "source")
```

Put `value` and `setting` columns in `routine.csv` and a `value` column in `screened.csv`. Replace `0.1` with your assay’s reporting increment or omit that argument.

```r
library(anchoredRI)
routine <- read.csv("routine.csv")
screened <- read.csv("screened.csv")$value
fit <- fit_anchored_ri(routine$value, routine$setting, screened, rounding_unit = 0.1)
fit[c("method", "status", "limits", "recommendation")]
```

## How to cite

Obstfeld A. anchoredRI: reference intervals from routine results and screened patients. R package version 0.4.0, 2026. [https://amromeo.github.io/anchoredRI/](https://amromeo.github.io/anchoredRI/) · [Commit `53d7b947bb40849b9b0b7556f97c676a0845928f`](https://github.com/amromeo/amromeo.github.io/commit/53d7b947bb40849b9b0b7556f97c676a0845928f). Fully AI-developed.

[Change log](NEWS.md) · [Questions or issues](https://github.com/amromeo/amromeo.github.io/issues)

## Earlier results

<a id="ribench-accuracy-by-routine-sample-size"></a>

### RIbench accuracy for the version 0.3.0 fitter

The earlier RIbench comparison used 40 independently screened patients. It compared the updated anchoredRI fitter before the safeguard with refineR.

| Routine results per dataset | Datasets | anchoredRI | refineR |
| ---: | ---: | ---: | ---: |
| 1,000 | 1,440 | 0.2334 | 0.4318 |
| 5,000 | 1,440 | 0.2114 | 0.3082 |
| 50,000 | 1,440 | 0.3744 | 0.2282 |

The scores are mean absolute endpoint z deviations before the safeguard. Lower scores indicate better accuracy. RIbench excludes dataset errors above 5 from this mean. Across both smaller sizes the scores were 0.2224 for anchoredRI and 0.3696 for refineR. At N = 50,000, anchoredRI had one implausible fit and selected five components in 941/1,440 fits (65.3%). The safeguard declined 305 fits and released 1,135. The score of 0.1307 among released fits applies to those 1,135 only. [Small-N results](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/results/REPORT.md) · [Large-N results](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/REPORT.md)

The checks reproduced the 5,220 archived reporting decisions. A same-environment repeat on seven datasets matched the original models and decisions. [Validation and reproduction](tools/safeguard/REPORT.md)

### Earlier results with version 0.2

Version 0.2 compared normal, lognormal and shifted-lognormal mixtures with one to three components when `transform = "auto"`. The shifted model estimates its offset jointly with the component shapes and service proportions. BIC includes the transformation Jacobian and the extra offset parameter, so scores are comparable on the measurement scale. Its point-fit default was `"identity"`.

`bootstrap_anchored_ri()` now resamples both routine services and screened observations, repeating the complete model selection and offset estimation. It records failed fits, boundary flags, offsets and selected models for every draw. Its default transform is `"auto"`.

On the original shifted example, the new model estimates an offset of 10.30 and limits of 15.70–49.51, against known limits of 15.63–49.97. Mean percentage endpoint error falls from 5.62% to 0.69%; the previously published refineR estimate has 2.35% error. The other three worked point estimates are unchanged. This is one example per distribution. The separate [fresh-data evaluation](assets/offset/validation-summary.csv) compares the two anchoredRI versions across 220 datasets, without rerunning refineR or estimating bootstrap-CI coverage in those datasets.

Across the 20 fresh datasets matching that shifted scenario, mean endpoint percentage error decreased from 7.12% to 1.04%, with improvement in all 20 pairs. The three unshifted scenarios were unchanged in all 60 pairs. The seven shifted scenarios without strong disease overlap improved in all 140 pairs. The shifted overlap scenario worsened on average, from 8.51% to 11.83%, although 9 of its 20 pairs improved. All 440 fits returned estimates without selected numerical-boundary flags. The new transformation fixes a distributional limitation; it does not resolve ambiguous separation of reference and disease components.

## Primary comparison: refineR as intended

The worked examples and the earlier clinic/ward results below used earlier package versions. To reproduce the published figures, run their scripts from the archived source. [Source for the earlier comparisons](https://github.com/amromeo/amromeo.github.io/tree/55eee86b3f4111197b0043cc062ef0b781dc2057/anchoredRI)

The [prepared-cohort comparison](https://amromeo.github.io/anchoredRI/prepared-cohort.html) uses the authors' practical workflow: the same 5,000 eligible routine patients for both methods, one age and sex group, one result per person, 200 bootstrap resamples and `medianBS` refineR estimates with 95% CIs. anchoredRI adds 40 independent screened patients and service labels. Direct benchmarks use empirical percentiles from 120 and 400 screened patients. The shifted profile uses `modBoxCox` for refineR. Both the original identity/log anchoredRI estimates and the new automatic estimates are shown.

The four worked examples show estimates, confidence intervals and original refineR diagnostic plots. The updated runner reuses the published comparator rows without changing their estimates or resampling. It adds 200 bootstrap draws for each new anchoredRI fit.

```sh
Rscript tools/build-guided-benchmark.R 8 1 200
# Add the updated model to those exact datasets and comparator results:
Rscript tools/build-offset-benchmark.R 6 200
# Paired old/new anchoredRI point fits on 20 fresh datasets per scenario:
Rscript tools/validate-offset.R 2 20
```

Read [BENCHMARK.md](BENCHMARK.md) for the prespecified distributions, metrics and sources. The [2023 practical guide](https://pubmed.ncbi.nlm.nih.gov/36610416/) and [package vignette](https://cloud.r-project.org/web/packages/refineR/vignettes/refineR_package.html) define the workflow. The original clinic/ward experiment below remains a stress test: its default pooled disease fraction is 61.25%, and the outpatient subset has only 250 observations.

### Apply the workflow to a prepared partition

`eligible` below is a data frame of compatible routine results with one row per patient and two eligible service labels. `screened` contains independent reference-eligible patients. Decide eligibility and the age and sex group before fitting. Select screened patients using clinical criteria.

```r
ref <- refineR::findRI(eligible$value, model = "BoxCox",
                      NBootstrap = 200, seed = 123, NCores = 2)
refineR::getRI(ref, pointEst = "medianBS", CIprop = 0.95)
plot(ref, pointEst = "medianBS")

anch <- anchoredRI::bootstrap_anchored_ri(
  eligible$value, eligible$service, screened$value, transform = "auto",
  n_boot = 200, seed = 123, workers = 2)
anch$limits
anch$ci
anch$fit$transform
anch$fit$shift
table(anch$bootstrap$transform)
```

Inspect the distribution and fitted curve before accepting the result. Consider refineR's `modBoxCox` for a skewed distribution shifted from zero. For anchoredRI, inspect the selected family, boundary flags and bootstrap model choices. The shifted-lognormal model estimates an offset but no Box-Cox power. `compare_ri_methods()` retains the older stress-test workflow.

## Earlier stress example: 40 screened and 1,000 routine patients

```r
library(anchoredRI)

d <- simulate_ri_study(scenario = 3, replicate = 1,
                       n_clinic = 250, n_ward = 750, n_anchor = 40)
fit <- fit_anchored_ri(d$routine$value, d$routine$setting, d$anchor)
fit$limits                 # reference limits on the original scale
fit$weights                # estimated component fractions in clinic and ward
fit$components             # model size selected by BIC
fit$at_bound               # optimizer reached a numeric parameter bound
d$truth                    # available for simulated data

# Legacy stress-test comparison; use the prepared-cohort script above for the new workflow.
compare_ri_methods(d$routine$value, d$routine$setting, d$anchor)
```

A screened patient contributes once. The simulations drew the groups separately. Their results describe a study design with different patients in the two sets. The simulator counts 40 **additional** screened patients. Eligibility should be decided from clinical information relevant to the analyte, while blinded to the measurement. Taking only values that fall inside an existing interval can narrow the apparent reference distribution.

## How the method works

1. Both settings share the component means and standard deviations. Each setting gets its own component proportions.
2. A routine observation contributes the probability density of the mixture for its setting. A screened observation contributes the probability density of the reference component. Each observation contributes once, with unit weight.
3. The default fitter compares one through five components across the eligible normal, lognormal and shifted-lognormal families by original-scale BIC. It reports the requested limits in measurement units when the selected fit passes both safeguards.

The offset is shared across components and constrained below every observation. Default bounds run from the minimum observation minus five routine MADs to that minimum minus 0.001 MAD; SD substitutes when MAD is zero. These numerical constraints avoid a fit collapsing onto the smallest result, but they can influence estimates. Optional `shift_bounds` must use measurement units. Selected boundary fits are retained and flagged. Mixture BIC and percentile bootstrap intervals are practical choices with limitations, especially when models fit poorly or screened samples are small. They do not guarantee correct model selection or confidence-interval coverage.

The initial prototype labeled the component enriched in the outpatient cohort. It struggled in simulations where both settings had the same share of patients with disease or disease shifted values in both directions. The next version introduced clinically screened observations and allowed two disease components. The original routine-only model and a fixed two-component comparison using screened patients remain available for investigating what each change contributes.

## Earlier clinic/ward stress experiment

| Method | Data used |
| --- | --- |
| refineR, all results | Released refineR on the entire routine sample |
| refineR, outpatients | The same release restricted to outpatient results |
| anchoredRI | All routine results plus independent screened results, selecting one to three components |
| Screened alone, selected size | Screened sample mean ± 1.96 sample SD |
| Screened alone, 120 / 400 | The same Gaussian calculation on fixed 120- and 400-patient screened cohorts |
| Settings alone | Two-component Gaussian mixture with setting-specific fractions and a clinic enrichment assumption |
| Pooled mixture | Same two-component family with pooled fractions, using setting labels afterward to name the putative reference component |
| With screened patients, two components | Joint screened/routine likelihood with the component count fixed at two |

The page displays R fits calculated in advance on 96 synthetic datasets. The table includes screened-only estimates from 120 and 400 patients alongside the model using the selected number of screened patients, using the same screening rule. Smaller screened cohorts are nested within the generated 400-patient cohort. Choose one of four case mixes, 1,000 or 4,000 routine results, 0/20/40/80/120/400 screened individuals, clinical screening or an intentionally restrictive 26–34 s cutoff, and an individual draw or the average of 12. The page uses the real `refineR::findRI(..., NBootstrap = 0L)` and `refineR::getRI()` point estimates. The average worst-limit error is the mean of each draw's larger absolute endpoint error. Averaging signed errors first would hide errors in opposite directions.

```r
# Reproduce a benchmark for 40 screened and 1,000 routine patients.
out <- benchmark_setting_ri(reps = 12, n_clinic = 250, n_ward = 750, n_anchor = 40)
summarize_ri_benchmark(out)
```

To regenerate the earlier stress-test results, run `Rscript tools/build-demo.R` from the package source directory with refineR installed. On Unix, an optional positive integer argument sets the worker count. The CSV records every fitted pair of endpoints. The website renders the CSV using plain JavaScript and can be served by any static HTTP server. [Data generator](R/simulate.R) · [Benchmark script](tools/build-demo.R) · [Recorded package versions](assets/versions.txt)

## Interpretation for pediatric laboratories

Both experiments assume one age and sex group, one assay, a common reference distribution across settings, independent patients and perfect clinical eligibility labels. The earlier stress generator uses Gaussian components, which favor anchoredRI's model. The prepared-cohort comparison also includes lognormal and shifted lognormal distributions. Neither models developmental trajectories, censored reporting, preanalytical differences or mislabeled clinical histories. The earlier 1,000-result stress design gives the outpatient-only comparator 250 observations, below refineR's 1,000-result caution threshold.

The appropriate number of screened patients still needs to be established for each intended application. The exercise tests errors in estimated reference limits. It does not evaluate classification of individual patients. A clinical study would need an independent group of screened patients covering both tails and an evaluation of uncertainty before reporting an interval. No CHOP data were used.

## Credit and license

We acknowledge Tatjana Ammer and colleagues for the published **refineR** method. The R package by Matthias Beck, Tatjana Ammer, Christopher M. Rank and André Schützenmeister supplies the baseline estimates and the search code for the checked fallback. The port retains the GPL-3 notices. [Ported code and license](inst/NOTICE). anchoredRI is an independent research project. [Ammer et al., *Scientific Reports* (2021)](https://pubmed.ncbi.nlm.nih.gov/34362961/) · [refineR on CRAN](https://cloud.r-project.org/web/packages/refineR/index.html)

Package code, simulations and website materials are licensed under the [GNU General Public License, version 3](LICENSE). The Hedvig Letters Serif font is distributed under its [SIL Open Font License](assets/FONT-LICENSE.txt).
