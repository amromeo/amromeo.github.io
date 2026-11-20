# anchoredRI

An R research package for asking whether a small group of independently screened children can improve a reference interval estimated from mixed outpatient and ward results. The main function estimates the central 95% interval for one analyte, assay method and comparable age group. All published experiments use synthetic data.

**Fully AI-developed.** A pediatric laboratory physician posed the question and directed the exploration. OpenAI ChatGPT designed and wrote the method implementation, simulation experiments, tests, package documentation and website. The estimator uses conventional maximum likelihood in R. Clinical performance has yet to be established.

[Interactive experiment](https://amromeo.github.io/anchoredRI/) · [Simulation results (CSV)](https://amromeo.github.io/anchoredRI/assets/fits.csv) · [R source archive](https://amromeo.github.io/anchoredRI/downloads/anchoredRI_0.1.0.tar.gz)

## Install from GitHub

The package occupies the `anchoredRI` directory of the author's public GitHub Pages repository. In R 4.1 or later, run

```r
install.packages("remotes")  # once
remotes::install_github("amromeo/amromeo.github.io", subdir = "anchoredRI")
```

Model fitting uses R's standard packages. To run the published comparator, install `refineR` with `install.packages("refineR")`.

## Try 40 screened patients with 1,000 routine results

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

compare_ri_methods(d$routine$value, d$routine$setting, d$anchor)
```

A screened patient contributes once. If you review charts from the same set of routine patients, remove those screened individuals from the unlabeled input before fitting. The simulator counts 40 **additional** screened patients. Eligibility should be decided from clinical information relevant to the analyte, while blinded to the measurement. Taking only values that fall inside an existing interval can narrow the apparent reference distribution.

## How the method works

1. Both settings share the component means and standard deviations. Each setting gets its own component proportions.
2. A routine observation contributes the probability density of the mixture for its setting. A screened observation contributes the probability density of the reference component. Each observation contributes once, with unit weight.
3. The fitter compares one, two and three Gaussian components using BIC. It reports the central 95% of the selected reference component. A log-scale fit is available for positive results. The reference limits always use the original units.

The initial prototype labeled the component enriched in the outpatient cohort. It struggled in simulations where both settings had the same fraction of affected patients or disease shifted values in both directions. The next version introduced clinically screened observations and allowed two disease components. The original routine-only model and a fixed two-component anchored comparison remain available for investigating what each change contributes.

## What the interactive experiment compares

| Method | Data used |
| --- | --- |
| refineR, all results | Released refineR on the entire routine sample |
| refineR, outpatients | The same release restricted to outpatient results |
| anchoredRI | All routine results plus independent screened results, selecting one to three components |
| Screened alone | Screened sample mean ± 1.96 sample SD |
| Settings alone | Two-component Gaussian mixture with setting-specific fractions and a clinic enrichment assumption |
| Pooled mixture | Same two-component family with pooled fractions, using setting labels afterward to name the putative reference component |
| Anchored, two components | Joint screened/routine likelihood with the component count fixed at two |

The page displays R fits calculated in advance on 96 synthetic datasets. Choose one of four case mixes, 1,000 or 4,000 routine results, 0/20/40/80 screened individuals, clinical screening or an intentionally restrictive 26–34 s cutoff, and an individual draw or the average of 12. The page uses the real `refineR::findRI(..., NBootstrap = 0L)` and `refineR::getRI()` point estimates. The average worst-limit error is the mean of each draw's larger absolute endpoint error. Averaging signed errors first would hide errors in opposite directions.

```r
# Reproduce a benchmark for 40 screened and 1,000 routine patients.
out <- benchmark_setting_ri(reps = 12, n_clinic = 250, n_ward = 750, n_anchor = 40)
summarize_ri_benchmark(out)
```

To regenerate every website result, run `Rscript tools/build-demo.R` from the package source directory with refineR installed. On Unix, an optional positive integer argument sets the worker count. The CSV records every fitted pair of endpoints. The website renders the CSV using plain JavaScript and can be served by any static HTTP server. [Data generator](R/simulate.R) · [Benchmark script](tools/build-demo.R) · [Recorded package versions](assets/versions.txt)

## Interpretation for pediatric laboratories

The generator assumes a single age group, one assay, common healthy physiology across settings, one independent result per individual and perfect clinical eligibility labels. It uses Gaussian mixture components, which favor the assumptions of the anchored model. It does not model developmental trajectories, repeated measurements, censored reporting, preanalytical differences or mislabeled clinical histories. The 1,000-result design gives the outpatient-only comparator 250 observations, below refineR's 1,000-result caution threshold.

The appropriate number of screened observations still needs to be established for each intended application. The exercise tests errors in estimated reference limits. It does not evaluate classification of individual patients. A clinical study would need an independent reference cohort covering both tails and an evaluation of uncertainty before reporting an interval. No CHOP data were used.

## Credit and license

We acknowledge Tatjana Ammer and colleagues for the published **refineR** method. The released R package by Matthias Beck, Tatjana Ammer, Christopher M. Rank and André Schützenmeister supplies the two baseline estimates. anchoredRI is an independent research project. [Ammer et al., *Scientific Reports* (2021)](https://doi.org/10.1038/s41598-021-95301-2) · [refineR on CRAN](https://cran.r-project.org/package=refineR)

Package code, simulations and website materials are licensed under the [GNU General Public License, version 3](LICENSE). The Hedvig Letters Serif font is distributed under its [SIL Open Font License](assets/FONT-LICENSE.txt).
