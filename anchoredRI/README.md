# anchoredRI

An R research package for asking whether a small group of independently screened children can improve a reference interval estimated from mixed outpatient and ward results. The main function estimates the central 95% interval for one analyte, assay method and comparable age group. All published experiments use synthetic data.

**Fully AI-developed.** A pediatric laboratory physician posed the question and directed the exploration. OpenAI ChatGPT designed and wrote the method implementation, simulation experiments, tests, package documentation and website. The estimator uses conventional maximum likelihood in R. Clinical performance has yet to be established.

[Prepared-cohort comparison](https://amromeo.github.io/anchoredRI/) · [Earlier stress tests](https://amromeo.github.io/anchoredRI/stress.html) · [Current results (CSV)](https://amromeo.github.io/anchoredRI/assets/offset/fits.csv) · [R source archive](https://amromeo.github.io/anchoredRI/downloads/anchoredRI_0.2.0.tar.gz)

## Proposed safeguard (Task F, awaiting review)

The proposed release adds `safeguard = TRUE` to `fit_anchored_ri()` and `bootstrap_anchored_ri()`. After the selected point fit and every bootstrap fit, two checks decide whether an interval can be reported. The Step B likelihood-ratio test compares the anchors under component 1 with their own normal MLE on the selected transform and fixed offset. A nominal chi-squared p-value with 2 df below 0.01 declines the interval. A lower limit strictly below the smallest routine value also declines it. The lower-limit check is skipped for one-sided limits.

A declined result has missing reported limits and a recorded reason. The reason is “routine data disagree with the screened sample” for the LRT flag or “estimated lower limit is below every observed result” for the lower-limit flag. Both reasons are retained if both checks flag. The fitted anchoredRI interval and the anchors-only interval remain available as diagnostics. A declined result recommends at least 120 screened patients or verifying a suitable published interval with 20 screened patients. For at least 50,000 routine results, it also suggests an established indirect method such as refineR.

Bootstrap output distinguishes declined fits from numerical failures and retains both checks for every draw. A declined point fit withholds the reported interval and confidence intervals. For an accepted point fit, confidence intervals use accepted, finite bootstrap draws and require at least 80% of requested draws. Report the accepted and declined counts and identify these as intervals conditional on passing the safeguard. This bootstrap policy is proposed for review before implementation.

Researchers can request `safeguard = FALSE` to return the fitted intervals regardless of these checks. The output records that the checks were disabled. With checks enabled, unavailable check statistics cause a separately identified technical withholding.

Heavily mixed routine data with many mildly abnormal patients, including ICU-only populations, can produce wrong reference intervals that pass both checks. Passing the safeguard does not establish clinical validity.

Validation will reproduce the committed B2 confirmation decisions and the archived RIbench decline counts using the same raw estimates and inputs. The tests will distinguish native mean-absolute-z implausibility from an extreme error at either endpoint. Current main requires two routine settings, while the archived RIbench studies used a one-group adapter. Supporting that one-group input and one-sided reporting is a prerequisite for reproducing the RIbench counts through the package entry points. Those additions will preserve the existing two-setting fits and central-95% defaults and will be checked against the frozen fits.

## Version 0.2: estimate the offset

`transform = "auto"` compares normal, lognormal and shifted-lognormal mixtures, each with one to three components. The shifted model estimates its offset jointly with the component shapes and service proportions. BIC includes the transformation Jacobian and the extra offset parameter, so scores are comparable on the measurement scale. The default point-fit transform remains `"identity"` for compatibility.

`bootstrap_anchored_ri()` now resamples both routine services and screened observations, repeating the complete model selection and offset estimation. It records failed fits, boundary flags, offsets and selected models for every draw. Its default transform is `"auto"`.

On the original shifted example, the new model estimates an offset of 10.30 and limits of 15.70–49.51, against known limits of 15.63–49.97. Mean percentage endpoint error falls from 5.62% to 0.69%; the previously published refineR estimate has 2.35% error. The other three worked point estimates are unchanged. This is one example per distribution. The separate [fresh-data evaluation](assets/offset/validation-summary.csv) compares the two anchoredRI versions across 220 datasets, without rerunning refineR or estimating bootstrap-CI coverage in those datasets.

Across the 20 fresh datasets matching that shifted scenario, mean endpoint percentage error decreased from 7.12% to 1.04%, with improvement in all 20 pairs. The three unshifted scenarios were unchanged in all 60 pairs. The seven shifted scenarios without strong disease overlap improved in all 140 pairs. The shifted overlap scenario worsened on average, from 8.51% to 11.83%, although 9 of its 20 pairs improved. All 440 fits returned estimates without selected numerical-boundary flags. The new transformation fixes a distributional limitation; it does not resolve ambiguous separation of reference and disease components.

## Install from GitHub

The package occupies the `anchoredRI` directory of the author's public GitHub Pages repository. In R 4.1 or later, run

```r
install.packages("remotes")  # once
remotes::install_github("amromeo/amromeo.github.io", subdir = "anchoredRI")
```

Model fitting uses R's standard packages. To run the published comparator, install `refineR` with `install.packages("refineR")`.

## Primary comparison: refineR as intended

The main page uses the authors' practical workflow: the same 5,000 eligible routine patients for both methods, a single physiological partition, one result per person, 200 bootstrap resamples and `medianBS` refineR estimates with 95% CIs. anchoredRI adds 40 independent screened patients and service labels. Direct benchmarks use empirical percentiles from 120 and 400 screened patients. The shifted profile uses `modBoxCox` for refineR. Both the original identity/log anchoredRI estimates and the new automatic estimates are shown.

The four worked examples show estimates, confidence intervals and original refineR diagnostic plots. The updated runner reuses the published comparator rows without changing their estimates or resampling. It adds 200 bootstrap draws for each new anchoredRI fit.

```sh
Rscript tools/build-guided-benchmark.R 8 1 200
# Add the updated model to those exact datasets and comparator results:
Rscript tools/build-offset-benchmark.R 6 200
# Paired old/new anchoredRI point fits on 20 fresh datasets per scenario:
Rscript tools/validate-offset.R 2 20
```

Read [BENCHMARK.md](BENCHMARK.md) for the prespecified distributions, metrics and sources. The [2023 practical guide](https://doi.org/10.1093/jalm/jfac101) and [package vignette](https://cloud.r-project.org/web/packages/refineR/vignettes/refineR_package.html) define the workflow. The original clinic/ward experiment below remains a stress test: its default pooled disease fraction is 61.25%, and the outpatient subset has only 250 observations.

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
3. The fitter compares one, two and three components using BIC. Automatic selection also compares identity, log and shifted-log families. It reports the central 95% of the selected reference component in the original units.

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

We acknowledge Tatjana Ammer and colleagues for the published **refineR** method. The released R package by Matthias Beck, Tatjana Ammer, Christopher M. Rank and André Schützenmeister supplies the two baseline estimates. anchoredRI is an independent research project. [Ammer et al., *Scientific Reports* (2021)](https://doi.org/10.1038/s41598-021-95301-2) · [refineR on CRAN](https://cran.r-project.org/package=refineR)

Package code, simulations and website materials are licensed under the [GNU General Public License, version 3](LICENSE). The Hedvig Letters Serif font is distributed under its [SIL Open Font License](assets/FONT-LICENSE.txt).
