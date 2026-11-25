\name{fit_anchored_ri}
\alias{fit_anchored_ri}
\title{Fit a reference component using independently screened observations}
\usage{
fit_anchored_ri(value, setting, anchor,
  transform = c("auto", "identity", "log", "shifted_log"),
  shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975),
  rounding_unit = NULL, fit_seed = 1L, fallback = TRUE,
  fallback_seed = 123L, pit_seed = 1L)
}
\arguments{
  \item{fallback}{Apply the checked anchorRefine fallback after an anchoredRI decline or technical withholding, including a numerical fitting failure. Default \code{TRUE}. The fallback runs only when \code{safeguard = TRUE}.}
  \item{fallback_seed, pit_seed}{Nonnegative integer seeds for the Stage 2 fit and the Stage 1 randomized probability integral transform. They are separate from \code{fit_seed}. The fit preserves the caller's random-number state.}
  \item{rounding_unit}{Assay reporting increment in measurement units. If \code{NULL}, the fitter finds and messages the largest step of the form \eqn{m \times 10^k} that divides every routine and screened value, with \eqn{m} in 1, 2, 2.5 or 5 and integer \eqn{k} from -8 to 3. A positive supplied step must divide all values. If no step fits or you supply zero, only the other component SD bounds apply. Every fit returns the unit used. Bootstrap draws reuse the point fit's unit.}
  \item{fit_seed}{Nonnegative integer for the optimizer starts, default 1. Model choice can occasionally change with the fit seed or computing environment. The bootstrap resampling seed is separate.}
  \item{safeguard}{Logical. Default TRUE applies the anchor likelihood-ratio and observed-minimum checks. FALSE returns raw estimates for research use and records that the checks are disabled.}
  \item{probs}{One probability for a one-sided limit or two increasing probabilities for an interval. Defaults to c(.025, .975). The observed-minimum check is skipped for one-sided requests.}

  \item{value}{Finite numeric routine results for a single analyte, method and comparable age group.}
  \item{setting}{A vector the same length as \code{value}, with one or two observed settings.}
  \item{anchor}{At least two finite values from independently screened reference individuals. Select them on clinical grounds while blinded to the target result. Exclude their records from the routine sample.}
  \item{transform}{\code{identity} fits Gaussian components to measured values; \code{log} uses their logarithms and requires positive values. \code{shifted_log} jointly estimates a common offset and Gaussian components on \code{log(value - offset)}. \code{auto} compares these families using original-scale BIC, omitting log when any value is nonpositive. Automatic selection is the default.}
  \item{shift_bounds}{Optional pair of increasing finite bounds for the offset, in measurement units. Both must be below every routine and screened observation. The default is the minimum of all observations minus five and 0.001 times the routine median absolute deviation (or SD if MAD is zero). Applies to shifted_log and auto only.}
}
\description{Fits a Gaussian mixture to one or two sets of routine results, with screened observations assigned to the reference component. Each setting has its own mixture proportions. The component means and scales are shared across settings.}
\details{The first fit uses the anchoredRI safeguard. A nominal anchor likelihood-ratio test p-value below 0.01 or a two-sided lower limit strictly below every routine value triggers a decline. The likelihood test holds the selected transform and offset fixed. A check that cannot be evaluated withholds the estimate.

When the first fit declines or withholds and \code{fallback = TRUE}, the function pools routine results across settings and tries the anchorRefine Stage 3 post-fit blend. Stage 2 guides the refineR 2.0.0 Box-Cox peak choice with the screened sample and uses the original refineR fit if the guided search cannot return a fit. The Box-Cox search and point quantiles use code from refineR 2.0.0 included in this package. Stage 3 blends the fitted reference mean and log spread toward their screened-sample estimates at weight 0.25. It retains the Stage 2 fit if the blend cannot be estimated. The Stage 1 check uses a randomized probability integral transform within each rounded measurement bin and declines at p < 0.01 or when an anchor's bin is outside the model's support. Fallback requires at least 40 nonnegative screened results, a positive reporting increment, and more than 10 nonnegative routine results. The Stage 1 exact null uses the actual screened count. The Stage 2 vote rules and the blend weight of 0.25 apply to all eligible counts. For fewer than 40 the function explains why it skips the fallback. Counts other than 40 give a warning because Task G evaluated both methods with exactly 40 screened patients. As in refineR, the Box-Cox fit uses the nonnegative routine results. A failed check or unusable limits leaves the interval withheld. If neither method reports, the recommendation is to obtain at least 120 screened patients or verify a suitable published interval with 20.

With more than 5,000 routine results, the point fit gives the sample-size message regardless of its reporting status. RIbench accuracy between 5,000 and 50,000 routine results has not been tested directly.

Heavily mixed populations containing many mildly abnormal patients, including ICU-only populations, can produce wrong intervals that pass both checks.

Each screened observation contributes one reference-component log density and each routine observation contributes one mixture log density. Arm B fits one through five components per eligible family using deterministic and seeded starts. The model with the lowest Bayesian information criterion (BIC) supplies the limits. Automatic selection includes transformation Jacobians and penalizes the additional offset parameter. All components share an offset; their shapes are shared between settings and their proportions may differ. This is a shifted-lognormal family, not a general Box-Cox transformation with an estimated power.

The Task E RIbench comparison before the safeguard used 40 independent screened anchors. At N = 1,000 and 5,000 combined, Arm B had the lowest native mean absolute endpoint z deviation among the methods tested (0.2224 vs 0.3696 for refineR before safeguarding). At N = 50,000, refineR scored 0.2282 and Arm B scored 0.3744.

On RIbench calcium dataset 870, the same fit seed selected five normal components in the archived environment and four shifted-lognormal components in the local environment. The archived fit triggered both safeguards. The local fit passed both.

Numeric bounds regularize component scales, mixture weights and offsets. The shift is estimated within a finite, data-derived support range to avoid collapse onto the smallest observation. In shifted-log fits, component SDs on the log scale range from 0.005 to 5; non-baseline weight logits range from -7 to 7. A selected boundary fit is retained and flagged. Inspect it and consider sensitivity to bounds before using the estimate. BIC is a model-selection heuristic for these constrained mixtures, not proof that the chosen family is correct. Population mismatch, selection on the test result and a poorly represented disease distribution can affect the limits. Use \code{bootstrap_anchored_ri} to assess sampling variation; neither function establishes clinical validity.}
\value{The reported \code{limits}, final \code{status}, and \code{method}. The method is \code{"anchoredRI"}, \code{"anchorRefine"}, or \code{NA} when no interval is returned. \code{anchoredRI_status} and \code{anchoredRI_reason} retain the first decision, even when the fallback reports. The original model fields and \code{safeguard} describe anchoredRI. \code{diagnostics$anchorRefine} contains the Stage 2 fit and peak selection, the blend, the candidate endpoints, and the Stage 1 check statistic, p-value, and decision. It is \code{NULL} when no fallback is needed or enabled. \code{diagnostics$anchored} and \code{diagnostics$anchors_only} retain the available anchoredRI estimates. Diagnostic estimates can be present even when \code{limits} is missing.

A list with named \code{limits} (requested percentiles in measurement units, missing when declined), \code{components} (selected count), \code{weights} (component by setting), \code{transform} (selected family), \code{requested_transform}, \code{shift}, \code{shift_bounds}, \code{bic_original}, \code{selected_bic}, and \code{at_bound}. \code{candidate_limits}, \code{candidate_shift} and \code{candidate_at_bound} give results by component count. The fitter stops if it cannot find a feasible starting fit for a required candidate. \code{selection}, returned for auto, contains all candidate families and counts. \code{rounding_unit} gives the unit used, including zero, and \code{fit_seed} gives the seed. \code{optimizer_status} describes the selected start's exit (\code{converged}, \code{iteration_limit}, \code{line_search_failure} or \code{other_exit}). \code{optimizer_code} and \code{optimizer_message} retain the native exit code and text. A reported fit that hits the iteration limit gives a warning and retains its limits.

\code{runner_up} is a one-row data frame with the next model's \code{transform}, \code{components}, \code{shift}, \code{bic_original} and \code{at_bound}. Its \code{bic_gap} equals the runner-up BIC minus the winning BIC. A zero gap indicates a tie. If only one candidate has a finite BIC, the runner-up fields are missing. \code{shift_at_bound} is also supplied for a selected shifted-log fit. Row 1 of \code{weights} represents the reference component. For compatibility, \code{bic} retains working-scale values for identity/log fits; use \code{bic_original} for comparisons across families.}
\examples{
d <- simulate_ri_study(n_clinic = 100, n_ward = 300, n_anchor = 40)
fit <- fit_anchored_ri(d$routine$value, d$routine$setting, d$anchor)
fit$limits
fit$weights
fit$rounding_unit
fit$optimizer_status
fit$runner_up
}
\seealso{\code{\link{bootstrap_anchored_ri}}, \code{\link{fit_setting_ri}}, \code{\link{compare_ri_methods}}}
