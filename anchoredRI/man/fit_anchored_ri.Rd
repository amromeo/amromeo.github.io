\name{fit_anchored_ri}
\alias{fit_anchored_ri}
\title{Fit a reference component using independently screened observations}
\usage{
fit_anchored_ri(value, setting, anchor, transform = c("identity", "log"))
}
\arguments{
  \item{value}{Finite numeric routine results for a single analyte, method and comparable age group.}
  \item{setting}{A vector the same length as \code{value}, with exactly two observed settings.}
  \item{anchor}{At least two finite values from independently screened reference individuals. Select them on clinical grounds while blinded to the target result. Exclude their records from the routine sample.}
  \item{transform}{Scale on which to fit the Gaussian components. Log requires strictly positive values. Reported limits are on the original scale.}
}
\description{Fits a Gaussian mixture to two sets of routine results, with screened observations assigned to the reference component. Each setting has its own mixture proportions. The component means and scales are shared across settings.}
\details{Each screened observation contributes one reference-component log density and each routine observation contributes one mixture log density. Models with one, two and three components are considered. The model with the lowest Bayesian information criterion (BIC) supplies the limits. Numeric bounds regularize component scales and mixture weights. Bound hits, population mismatch, selection on the test result and a poorly represented disease distribution can affect the limits. This research estimator does not quantify sampling uncertainty.}
\value{A list with named \code{limits} (2.5th and 97.5th percentiles), \code{components} (selected count), \code{weights} (component by setting), \code{bic} (on the standardized working scale), and \code{at_bound} for the selected model. \code{candidate_limits} and \code{candidate_at_bound} give the results for each component count. A failed candidate has infinite BIC and missing limits. Row 1 of \code{weights} represents the reference component.}
\examples{
d <- simulate_ri_study(n_clinic = 100, n_ward = 300, n_anchor = 40)
fit <- fit_anchored_ri(d$routine$value, d$routine$setting, d$anchor)
fit$limits
fit$weights
}
\seealso{\code{\link{fit_setting_ri}}, \code{\link{compare_ri_methods}}}
