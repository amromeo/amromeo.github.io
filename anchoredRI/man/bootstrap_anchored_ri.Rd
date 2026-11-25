\name{bootstrap_anchored_ri}
\alias{bootstrap_anchored_ri}
\title{Bootstrap the complete anchored reference-interval fit}
\usage{
bootstrap_anchored_ri(value, setting, anchor,
  transform = c("auto", "identity", "log", "shifted_log"),
  n_boot = 200L, conf = .95, seed = 123L, workers = 1L,
  shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975),
  rounding_unit = NULL, fit_seed = 1L)
}
\arguments{
  \item{rounding_unit}{Assay reporting increment in measurement units. If \code{NULL}, the fitter finds and messages the largest step of the form \eqn{m \times 10^k} that divides every routine and screened value, with \eqn{m} in 1, 2, 2.5 or 5 and integer \eqn{k} from -8 to 3. A positive supplied step must divide all values. If no step fits or you supply zero, only the other component SD bounds apply. Every fit returns the unit used. Bootstrap draws reuse the point fit's unit.}
  \item{fit_seed}{Nonnegative integer for the optimizer starts, default 1. Model choice can occasionally change with the fit seed or computing environment. The bootstrap resampling seed is separate.}
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
\details{With the safeguard enabled, decline when the nominal anchor likelihood-ratio test has p < 0.01 (chi-squared, 2 df) or when a two-sided lower limit is strictly below every routine result. The likelihood test holds the selected transform and offset fixed. Unavailable check statistics cause technical withholding. Reported limits are missing when withheld; raw anchored and anchors-only estimates remain under \code{diagnostics}, never as substitute reported intervals. The output includes \code{status}, \code{reason}, \code{recommendation}, \code{probs}, and \code{safeguard} with both flags, LRT statistic and p-value. A declined result recommends at least 120 screened patients or verification of a suitable published interval with 20. For any call with more than 5,000 routine results, a message explains that RIbench accuracy between 5,000 and 50,000 routine results has not been tested directly. This message applies to all reporting statuses.

Heavily mixed populations containing many mildly abnormal patients, including ICU-only populations, can produce wrong intervals that pass both checks.

The full-data point estimate and resamples use Arm B alone. Use \code{\link{fit_anchored_ri}} for the point fit with the checked fallback. The rounding and size messages apply once per call. A reported point fit that hits the iteration limit gives a warning and retains its interval. A separate warning gives the count of reported bootstrap draws that hit the limit. Those draws remain eligible for the confidence interval. Percentile confidence intervals use accepted bootstrap fits and are conditional on passing the safeguard. A declined or technically withheld point fit also withholds confidence intervals. Failed fits and selected numerical boundary fits are recorded. If fewer than 80 percent are accepted, confidence limits are returned as missing with a warning. This threshold is a reporting rule, not evidence of valid coverage. With small screened samples, distribution mismatch or unstable model selection, percentile intervals may have poor coverage. Individuals must contribute one result each; this function does not implement clustered patient resampling. Parallel and sequential runs use the same per-replicate seeds.}
\value{A list containing \code{fit}, \code{limits}, the reporting \code{method}, a \code{ci} matrix with rows lower/upper and columns low/high, \code{confidence}, a \code{bootstrap} data frame, and the counts \code{requested}, \code{successful}, \code{failed} and \code{boundary_fits}. The output also contains \code{declined_draws}, \code{withheld_draws}, and \code{ci_conditional_on_safeguard}. The \code{failed} count is reserved for numerical fitting failures. Every draw retains both diagnostic estimates and safeguard results. The point fit's \code{rounding_unit}, \code{optimizer_status}, \code{optimizer_code}, \code{optimizer_message} and \code{runner_up} are also available at the top level. Each bootstrap row records limits, selected family, offset, component count, boundary flag, status and error. Fitted draws also give the winning optimizer's status, code and message, the \code{runner_up_transform}, \code{runner_up_components}, \code{runner_up_bic} and \code{bic_gap} (runner-up BIC minus winning BIC).}
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
