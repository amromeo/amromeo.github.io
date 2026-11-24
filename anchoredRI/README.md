# anchoredRI

An R research package for asking whether a small group of independently screened children can improve a reference interval estimated from mixed outpatient and ward results. The main function estimates the central 95% interval for one analyte, assay method and comparable age group. All published experiments use synthetic data.

**Fully AI-developed.** A pediatric laboratory physician posed the question and directed the exploration. OpenAI ChatGPT designed and wrote the method implementation, simulation experiments, tests, package documentation and website. The package uses maximum likelihood with an R interface and a compiled C++ fitter. Clinical performance has yet to be established.

[Prepared-cohort comparison](https://amromeo.github.io/anchoredRI/) · [Earlier stress tests](https://amromeo.github.io/anchoredRI/stress.html) · [Original prepared-cohort results (CSV)](https://amromeo.github.io/anchoredRI/assets/offset/fits.csv) · [Download anchoredRI 0.3.0](https://amromeo.github.io/anchoredRI/downloads/anchoredRI_0.3.0.tar.gz) · [Previous versions and installation](https://amromeo.github.io/anchoredRI/stress.html#install)

## Fitting and reporting a reference interval

The point estimate and bootstrap use Task E Arm B. The fitter considers one through five components for each eligible family (normal, lognormal and shifted lognormal) and selects the lowest BIC on the measurement scale. Every fit runs the two safeguards below. The default interval covers the central 95% of the screened component. Use `probs = .95` for a one-sided upper limit.

Supply the assay’s reporting increment with `rounding_unit` when known. It sets a minimum component SD in measurement units. If you omit it, the fitter finds the largest step of the form {1, 2, 2.5, 5} × 10^k, for integer k from −8 through 3, that divides every routine and screened value. For example, results reported to the nearest 0.5 can give an inferred unit of 0.5. The fitter sends a message with the inferred unit and always returns the unit used in `rounding_unit`. A value of zero means that only the other SD constraints apply. You can supply zero for continuous data, and the fitter uses zero if no shared step fits. All bootstrap draws use the point fit’s unit.

`runner_up` identifies the alternative model with the next-lowest original-scale BIC. Its `bic_gap` is the runner-up BIC minus the winning BIC. Zero indicates a tie. The fit also returns the winning optimizer’s `optimizer_status`, numeric `optimizer_code`, and original `optimizer_message`. A reported fit that reached the iteration limit produces a warning and retains its interval. Line-search failures have their own status. Bootstrap output includes these diagnostics for the point fit and records the optimizer exit and runner-up for each fitted draw.

`fit_seed` (default 1) controls the optimizer starts. The bootstrap `seed` controls patient resampling. A change in fit seed or computing environment can occasionally change the selected model and the reporting decision. On RIbench calcium dataset 870, the same fit seed selected five normal components in the archived run and four shifted-lognormal components in the local run. Both safeguards flagged the archived fit. The local fit passed both. [Validation and reproduction](tools/safeguard/REPORT.md)

### RIbench accuracy by routine sample size

With 40 independently screened anchors, anchoredRI had the lowest error of the methods Task E tested at N = 1,000 and 5,000. RIbench accuracy between 5,000 and 50,000 routine results has not been tested. At N = 50,000, refineR was more accurate. For datasets that large, use refineR or anchorRefine. For every fit with more than 5,000 routine results, the package sends a message about this limitation whether it reports or withholds the interval. Task E measured refineR directly. [anchorRefine](https://github.com/amromeo/anchorRefine)

| Routine results per dataset | Datasets | anchoredRI (Arm B) | refineR |
| ---: | ---: | ---: | ---: |
| 1,000 | 1,440 | 0.2334 | 0.4318 |
| 5,000 | 1,440 | 0.2114 | 0.3082 |
| 50,000 | 1,440 | 0.3744 | 0.2282 |

The scores are mean absolute endpoint z deviations before the safeguard. Lower scores indicate better accuracy. RIbench excludes dataset errors above 5 from this mean. Across both smaller sizes the scores were 0.2224 for anchoredRI and 0.3696 for refineR. At N = 50,000, Arm B had one implausible fit and selected five components in 941/1,440 fits (65.3%). The safeguard declined 305 fits and released 1,135. The score of 0.1307 among released fits applies to those 1,135 only. [Small-N results](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/results/REPORT.md) · [Large-N results](https://github.com/amromeo/amromeo.github.io/blob/8bcd49c3c7977c50982524b9dbc48d16bdfa5914/anchoredRI/tools/boxcox/accuracy/large/results/REPORT.md)

## Safeguard

Version 0.3 adds `safeguard = TRUE` to `fit_anchored_ri()` and `bootstrap_anchored_ri()`. After the selected point fit and every bootstrap fit, two checks decide whether an interval can be reported. The Step B likelihood-ratio test compares the anchors under component 1 with their own normal MLE on the selected transform and fixed offset. A nominal chi-squared p-value with 2 df below 0.01 declines the interval. A lower limit strictly below the smallest routine value also declines it. The lower-limit check is skipped for one-sided limits.

A declined result has missing reported limits and a recorded reason. The reason is “routine data disagree with the screened sample” for the LRT flag or “estimated lower limit is below every observed result” for the lower-limit flag. Both reasons are retained if both checks flag. The fitted anchoredRI interval and the anchors-only interval remain available under `diagnostics$anchored` and `diagnostics$anchors_only`. Read `status`, `reason`, `recommendation`, and `safeguard` for the reporting decision. The anchors-only interval is never substituted into `limits`. A declined result recommends at least 120 screened patients or verifying a suitable published interval with 20 screened patients. For at least 50,000 routine results, it recommends refineR or anchorRefine.

Bootstrap output distinguishes declined fits from numerical failures and retains both checks for every draw. A declined point fit withholds the reported interval and confidence intervals. For an accepted point fit, confidence intervals use accepted, finite bootstrap draws and require at least 80% of requested draws. Report the accepted and declined counts and identify these as intervals conditional on passing the safeguard.

Researchers can set `safeguard = FALSE` to return the fitted limits. The package skips both checks, sets their flags and the LRT p-value to `NA`, and records `safeguard$enabled = FALSE`. With checks enabled, an unavailable statistic gives a `withheld` status.

Heavily mixed routine data with many mildly abnormal patients, including ICU-only populations, can produce wrong reference intervals that pass both checks. Passing the safeguard does not establish clinical validity.

The safeguard reproduced all 5,220 archived Arm B reporting decisions. On seven datasets refitted on the same runtime, the package and the original Task E fitter chose the same models and made the same decisions. [Validation and reproduction](tools/safeguard/REPORT.md)

## Earlier results with version 0.2

Version 0.2 compared normal, lognormal and shifted-lognormal mixtures with one to three components when `transform = "auto"`. The shifted model estimates its offset jointly with the component shapes and service proportions. BIC includes the transformation Jacobian and the extra offset parameter, so scores are comparable on the measurement scale. Its point-fit default was `"identity"`.

`bootstrap_anchored_ri()` now resamples both routine services and screened observations, repeating the complete model selection and offset estimation. It records failed fits, boundary flags, offsets and selected models for every draw. Its default transform is `"auto"`.

On the original shifted example, the new model estimates an offset of 10.30 and limits of 15.70–49.51, against known limits of 15.63–49.97. Mean percentage endpoint error falls from 5.62% to 0.69%; the previously published refineR estimate has 2.35% error. The other three worked point estimates are unchanged. This is one example per distribution. The separate [fresh-data evaluation](assets/offset/validation-summary.csv) compares the two anchoredRI versions across 220 datasets, without rerunning refineR or estimating bootstrap-CI coverage in those datasets.

Across the 20 fresh datasets matching that shifted scenario, mean endpoint percentage error decreased from 7.12% to 1.04%, with improvement in all 20 pairs. The three unshifted scenarios were unchanged in all 60 pairs. The seven shifted scenarios without strong disease overlap improved in all 140 pairs. The shifted overlap scenario worsened on average, from 8.51% to 11.83%, although 9 of its 20 pairs improved. All 440 fits returned estimates without selected numerical-boundary flags. The new transformation fixes a distributional limitation; it does not resolve ambiguous separation of reference and disease components.

## Install from GitHub

To install the current version from the main branch in R 4.1 or later, run

```r
install.packages("remotes")  # once
remotes::install_github("amromeo/amromeo.github.io", subdir = "anchoredRI",
                        ref = "main")
```

Building the package from source requires a C++17 compiler and R’s package build tools. Install the checkout before sourcing its R files directly. To run the published comparator, install `refineR` with `install.packages("refineR")`.

## Primary comparison: refineR as intended

The worked examples and the earlier clinic/ward results below used earlier package versions. To reproduce the published figures, run their scripts from the archived source. [Source for the earlier comparisons](https://github.com/amromeo/amromeo.github.io/tree/55eee86b3f4111197b0043cc062ef0b781dc2057/anchoredRI)

The main page uses the authors' practical workflow: the same 5,000 eligible routine patients for both methods, a single physiological partition, one result per person, 200 bootstrap resamples and `medianBS` refineR estimates with 95% CIs. anchoredRI adds 40 independent screened patients and service labels. Direct benchmarks use empirical percentiles from 120 and 400 screened patients. The shifted profile uses `modBoxCox` for refineR. Both the original identity/log anchoredRI estimates and the new automatic estimates are shown.

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

`eligible` below is a data frame of compatible routine results with one row per patient and two eligible service labels. `screened` contains independent reference-eligible patients. Decide eligibility and physiological partitions before fitting, without filtering by an existing reference interval.

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

A screened patient contributes once. If you review charts from the same set of routine patients, remove those screened individuals from the unlabeled input before fitting. The simulator counts 40 **additional** screened patients. Eligibility should be decided from clinical information relevant to the analyte, while blinded to the measurement. Taking only values that fall inside an existing interval can narrow the apparent reference distribution.

## How the method works

1. Both settings share the component means and standard deviations. Each setting gets its own component proportions.
2. A routine observation contributes the probability density of the mixture for its setting. A screened observation contributes the probability density of the reference component. Each observation contributes once, with unit weight.
3. The default fitter compares one through five components across the eligible normal, lognormal and shifted-lognormal families by original-scale BIC. It reports the requested limits in measurement units when the selected fit passes both safeguards.

The offset is shared across components and constrained below every observation. Default bounds run from the minimum observation minus five routine MADs to that minimum minus 0.001 MAD; SD substitutes when MAD is zero. These numerical constraints avoid a fit collapsing onto the smallest result, but they can influence estimates. Optional `shift_bounds` must use measurement units. Selected boundary fits are retained and flagged. Mixture BIC and percentile bootstrap intervals are practical choices with limitations, especially when models fit poorly or screened samples are small. They do not guarantee correct model selection or confidence-interval coverage.

The initial prototype labeled the component enriched in the outpatient cohort. It struggled in simulations where both settings had the same fraction of affected patients or disease shifted values in both directions. The next version introduced clinically screened observations and allowed two disease components. The original routine-only model and a fixed two-component anchored comparison remain available for investigating what each change contributes.

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
| Anchored, two components | Joint screened/routine likelihood with the component count fixed at two |

The page displays R fits calculated in advance on 96 synthetic datasets. The table includes screened-only estimates from 120 and 400 patients alongside the selected-size anchored model, using the same screening rule. Smaller screened cohorts are nested within the generated 400-patient cohort. Choose one of four case mixes, 1,000 or 4,000 routine results, 0/20/40/80/120/400 screened individuals, clinical screening or an intentionally restrictive 26–34 s cutoff, and an individual draw or the average of 12. The page uses the real `refineR::findRI(..., NBootstrap = 0L)` and `refineR::getRI()` point estimates. The average worst-limit error is the mean of each draw's larger absolute endpoint error. Averaging signed errors first would hide errors in opposite directions.

```r
# Reproduce a benchmark for 40 screened and 1,000 routine patients.
out <- benchmark_setting_ri(reps = 12, n_clinic = 250, n_ward = 750, n_anchor = 40)
summarize_ri_benchmark(out)
```

To regenerate the earlier stress-test results, run `Rscript tools/build-demo.R` from the package source directory with refineR installed. On Unix, an optional positive integer argument sets the worker count. The CSV records every fitted pair of endpoints. The website renders the CSV using plain JavaScript and can be served by any static HTTP server. [Data generator](R/simulate.R) · [Benchmark script](tools/build-demo.R) · [Recorded package versions](assets/versions.txt)

## Interpretation for pediatric laboratories

Both experiments assume a single physiological partition, one assay, a common reference distribution across settings, independent patients and perfect clinical eligibility labels. The earlier stress generator uses Gaussian components, which favor anchoredRI's model. The prepared-cohort comparison also includes lognormal and shifted lognormal distributions. Neither models developmental trajectories, censored reporting, preanalytical differences or mislabeled clinical histories. The earlier 1,000-result stress design gives the outpatient-only comparator 250 observations, below refineR's 1,000-result caution threshold.

The appropriate number of screened observations still needs to be established for each intended application. The exercise tests errors in estimated reference limits. It does not evaluate classification of individual patients. A clinical study would need an independent reference cohort covering both tails and an evaluation of uncertainty before reporting an interval. No CHOP data were used.

## Credit and license

We acknowledge Tatjana Ammer and colleagues for the published **refineR** method. The released R package by Matthias Beck, Tatjana Ammer, Christopher M. Rank and André Schützenmeister supplies the two baseline estimates. anchoredRI is an independent research project. [Ammer et al., *Scientific Reports* (2021)](https://pubmed.ncbi.nlm.nih.gov/34362961/) · [refineR on CRAN](https://cran.r-project.org/package=refineR)

Package code, simulations and website materials are licensed under the [GNU General Public License, version 3](LICENSE). The Hedvig Letters Serif font is distributed under its [SIL Open Font License](assets/FONT-LICENSE.txt).
