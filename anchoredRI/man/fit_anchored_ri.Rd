\name{fit_anchored_ri}
\alias{fit_anchored_ri}
\title{Fit a reference component using independently screened observations}
\usage{
fit_anchored_ri(value, setting, anchor,
  transform = c("identity", "log", "shifted_log", "auto"),
  shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975))
}
\arguments{
  \item{safeguard}{Logical. Default TRUE applies the anchor likelihood-ratio and observed-minimum checks. FALSE returns raw estimates for research use and records that the checks are disabled.}
  \item{probs}{One probability for a one-sided limit or two increasing probabilities for an interval. Defaults to c(.025, .975). The observed-minimum check is skipped for one-sided requests.}

  \item{value}{Finite numeric routine results for a single analyte, method and comparable age group.}
  \item{setting}{A vector the same length as \code{value}, with one or two observed settings.}
  \item{anchor}{At least two finite values from independently screened reference individuals. Select them on clinical grounds while blinded to the target result. Exclude their records from the routine sample.}
  \item{transform}{\code{identity} fits Gaussian components to measured values; \code{log} uses their logarithms and requires positive values. \code{shifted_log} jointly estimates a common offset and Gaussian components on \code{log(value - offset)}. \code{auto} compares these families using original-scale BIC, omitting log when any value is nonpositive. The identity default is retained for compatibility.}
  \item{shift_bounds}{Optional pair of increasing finite bounds for the offset, in measurement units. Both must be below every routine and screened observation. The default is the minimum of all observations minus five and 0.001 times the routine median absolute deviation (or SD if MAD is zero). Applies to shifted_log and auto only.}
}
\description{Fits a Gaussian mixture to one or two sets of routine results, with screened observations assigned to the reference component. Each setting has its own mixture proportions. The component means and scales are shared across settings.}
\details{With the safeguard enabled, decline when the nominal anchor likelihood-ratio test has p < 0.01 (chi-squared, 2 df) or when a two-sided lower limit is strictly below every routine result. The likelihood test holds the selected transform and offset fixed. Unavailable check statistics cause technical withholding. Reported limits are missing when withheld; raw anchored and anchors-only estimates remain under \code{diagnostics}, never as substitute reported intervals. The output includes \code{status}, \code{reason}, \code{recommendation}, \code{probs}, and \code{safeguard} with both flags, LRT statistic and p-value. A declined result recommends at least 120 screened patients or verification of a suitable published interval with 20. For at least 50,000 routine values it also suggests an established indirect method such as refineR.

Heavily mixed populations containing many mildly abnormal patients, including ICU-only populations, can produce wrong intervals that pass both checks.

Each screened observation contributes one reference-component log density and each routine observation contributes one mixture log density. Models with one, two and three components are considered. The model with the lowest Bayesian information criterion (BIC) supplies the limits. Automatic selection includes transformation Jacobians and penalizes the additional offset parameter. All components share an offset; their shapes are shared between settings and their proportions may differ. This is a shifted-lognormal family, not a general Box-Cox transformation with an estimated power.

Numeric bounds regularize component scales, mixture weights and offsets. The shift is estimated within a finite, data-derived support range to avoid collapse onto the smallest observation. In shifted-log fits, component SDs on the log scale range from 0.005 to 5; non-baseline weight logits range from -7 to 7. A selected boundary fit is retained and flagged. Inspect it and consider sensitivity to bounds before using the estimate. BIC is a model-selection heuristic for these constrained mixtures, not proof that the chosen family is correct. Population mismatch, selection on the test result and a poorly represented disease distribution can affect the limits. Use \code{bootstrap_anchored_ri} to assess sampling variation; neither function establishes clinical validity.}
\value{A list with named \code{limits} (requested percentiles in measurement units, missing when declined), \code{components} (selected count), \code{weights} (component by setting), \code{transform} (selected family), \code{requested_transform}, \code{shift}, \code{shift_bounds}, \code{bic_original}, \code{selected_bic}, and \code{at_bound}. \code{candidate_limits}, \code{candidate_shift} and \code{candidate_at_bound} give results by component count. A failed candidate has infinite BIC and missing limits. \code{selection}, returned for auto, contains all candidate families and counts. \code{shift_at_bound} is also supplied for a selected shifted-log fit. Row 1 of \code{weights} represents the reference component. For compatibility, \code{bic} retains working-scale values for identity/log fits; use \code{bic_original} for comparisons across families.}
\examples{
d <- simulate_ri_study(n_clinic = 100, n_ward = 300, n_anchor = 40)
fit <- fit_anchored_ri(d$routine$value, d$routine$setting, d$anchor)
fit$limits
fit$weights
}
\seealso{\code{\link{bootstrap_anchored_ri}}, \code{\link{fit_setting_ri}}, \code{\link{compare_ri_methods}}}
