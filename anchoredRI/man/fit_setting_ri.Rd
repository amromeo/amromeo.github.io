\name{fit_setting_ri}
\alias{fit_setting_ri}
\title{Fit a two-component mixture to routine measurements from two settings}
\usage{
fit_setting_ri(value, setting, enriched_setting,
               transform = c("identity", "log"), mixing = c("by_setting", "pooled"))
}
\arguments{
  \item{value}{Finite numeric results.}
  \item{setting}{A setting for each result, with exactly two distinct observed values.}
  \item{enriched_setting}{The setting expected to contain a higher reference fraction.}
  \item{transform}{Working scale. Choose \code{"log"} only for positive results.}
  \item{mixing}{Fit separate proportions by setting or fit one pooled proportion.}
}
\description{Fits two Gaussian densities on the working scale using multiple starting points and bounded maximum likelihood. After fitting, the candidate reference component is the one with a higher estimated membership fraction in the designated enriched setting. That assumption can fail when case mix is similar in both settings.}
\value{A list with \code{limits}, \code{candidate}, \code{component_summary}, \code{enrichment_gap}, \code{at_bound}, \code{loglik_standardized} and \code{converged_starts}. The enrichment gap is the absolute between-setting difference in estimated membership.}
\examples{
d <- simulate_ri_case(2, 1, n_clinic = 120, n_ward = 360)
fit_setting_ri(d$value, d$setting, "clinic")$limits
}
