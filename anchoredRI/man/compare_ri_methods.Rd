\name{compare_ri_methods}
\alias{compare_ri_methods}
\alias{benchmark_setting_ri}
\alias{summarize_ri_benchmark}
\title{Compare reference interval methods on shared data}
\usage{
compare_ri_methods(value, setting, anchor = numeric(),
                   enriched_setting = "clinic", include_stock = TRUE)
benchmark_setting_ri(reps = 12L, scenarios = 1:4, include_stock = TRUE,
                     n_clinic = 1000L, n_ward = 3000L, n_anchor = 40L)
summarize_ri_benchmark(results)
}
\arguments{
  \item{value, setting, anchor}{Routine results, two setting labels and an optional independent screened sample as in \code{\link{fit_anchored_ri}}.}
  \item{enriched_setting}{Setting expected to have a higher reference fraction, also used for the outpatient-only comparator.}
  \item{include_stock}{Include published \pkg{refineR} fits using all results and the enriched setting alone. Requires that package.}
  \item{reps}{Number of independent replicate identifiers per scenario.}
  \item{scenarios}{Numeric scenario identifiers from 1 through 4.}
  \item{n_clinic, n_ward, n_anchor}{Counts for \code{\link{simulate_ri_study}}. Use zero anchors to compare routine-only methods.}
  \item{results}{A data frame returned by \code{benchmark_setting_ri}.}
}
\description{Compare fitted reference limits on the same routine data. The comparators include setting-specific and pooled two-component fits and, when screened data exist, a BIC-selected anchored model, the anchored model restricted to two components, and screened-only Gaussian limits. Published \pkg{refineR} is invoked with zero bootstrap replicates for its point estimates.}
\value{\code{compare_ri_methods} returns a row per method with lower and upper limits, component count, bound status and enrichment gap where applicable. \code{benchmark_setting_ri} adds scenario, replicate and signed endpoint errors (estimate minus truth). \code{summarize_ri_benchmark} returns mean absolute lower and upper error and the mean of each draw's larger absolute endpoint error.}
\references{Ammer T, Schuetzenmeister A, Prokosch HU, Rauh M, Rank CM, Zierk J (2021). refineR: A Novel Algorithm for Reference Interval Estimation from Real-World Data. \emph{Scientific Reports} 11, 16039. \doi{10.1038/s41598-021-95301-2}.}
\examples{
d <- simulate_ri_study(n_clinic = 100, n_ward = 300)
compare_ri_methods(d$routine$value, d$routine$setting, d$anchor,
                   include_stock = FALSE)
}
