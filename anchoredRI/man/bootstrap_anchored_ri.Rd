\name{bootstrap_anchored_ri}
\alias{bootstrap_anchored_ri}
\title{Bootstrap the complete anchored reference-interval fit}
\usage{
bootstrap_anchored_ri(value, setting, anchor,
  transform = c("auto", "identity", "log", "shifted_log"),
  n_boot = 200L, conf = .95, seed = 123L, workers = 1L,
  shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975))
}
\arguments{
  \item{safeguard}{Logical. Default TRUE applies the anchor likelihood-ratio and observed-minimum checks. FALSE returns raw estimates for research use and records that the checks are disabled.}
  \item{probs}{One probability for a one-sided limit or two increasing probabilities for an interval. Defaults to c(.025, .975). The observed-minimum check is skipped for one-sided requests.}

  \item{value}{Finite numeric routine results from one physiological partition and comparable assay.}
  \item{setting}{One or two service labels, one per routine result.}
  \item{anchor}{Independent clinically screened reference measurements; at least two values.}
  \item{transform}{Family selection passed to \code{fit_anchored_ri}. Defaults to automatic selection.}
  \item{n_boot}{Number of resamples, at least two. The default is 200; very small counts are only useful for software checks.}
  \item{conf}{Confidence level, strictly between zero and one.}
  \item{seed}{Nonnegative integer seed. The caller's random-number state is restored.}
  \item{workers}{Positive integer number of parallel Unix workers. Windows uses one worker and warns if more were requested.}
  \item{shift_bounds}{Optional fixed offset bounds passed to every fit. If NULL, the data-derived bounds are recalculated within each resample.}
}
\description{Resamples routine observations within each setting and screened observations independently. Every resample repeats component selection, family selection when requested, and offset estimation.}
\details{With the safeguard enabled, decline when the nominal anchor likelihood-ratio test has p < 0.01 (chi-squared, 2 df) or when a two-sided lower limit is strictly below every routine result. The likelihood test holds the selected transform and offset fixed. Unavailable check statistics cause technical withholding. Reported limits are missing when withheld; raw anchored and anchors-only estimates remain under \code{diagnostics}, never as substitute reported intervals. The output includes \code{status}, \code{reason}, \code{recommendation}, \code{probs}, and \code{safeguard} with both flags, LRT statistic and p-value. A declined result recommends at least 120 screened patients or verification of a suitable published interval with 20. For at least 50,000 routine values it also suggests an established indirect method such as refineR.

Heavily mixed populations containing many mildly abnormal patients, including ICU-only populations, can produce wrong intervals that pass both checks.

The point estimate comes from the full dataset. Percentile confidence intervals use accepted bootstrap fits and are conditional on passing the safeguard. A declined or technically withheld point fit also withholds confidence intervals. Failed fits and selected numerical boundary fits are recorded. If fewer than 80 percent are accepted, confidence limits are returned as missing with a warning. This threshold is a reporting rule, not evidence of valid coverage. With small screened samples, distribution mismatch or unstable model selection, percentile intervals may have poor coverage. Individuals must contribute one result each; this function does not implement clustered patient resampling. Parallel and sequential runs use the same per-replicate seeds.}
\value{A list containing \code{fit}, \code{limits}, a \code{ci} matrix with rows lower/upper and columns low/high, \code{confidence}, a \code{bootstrap} data frame, and the counts \code{requested}, \code{successful}, \code{failed} and \code{boundary_fits}. The output also contains \code{declined_draws}, \code{withheld_draws}, and \code{ci_conditional_on_safeguard}. The \code{failed} count is reserved for numerical fitting failures. Every draw retains both diagnostic estimates and safeguard results. Each bootstrap row reports limits, selected family, offset, component count, boundary flag, status and error.}
\examples{
\dontrun{
d <- simulate_ri_study()
b <- bootstrap_anchored_ri(d$routine$value, d$routine$setting, d$anchor,
                         n_boot = 200, workers = 2)
b$limits
b$ci
table(b$bootstrap$transform)
}
}
\seealso{\code{\link{fit_anchored_ri}}}
