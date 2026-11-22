\name{bootstrap_anchored_ri}
\alias{bootstrap_anchored_ri}
\title{Bootstrap the complete anchored reference-interval fit}
\usage{
bootstrap_anchored_ri(value, setting, anchor,
  transform = c("auto", "identity", "log", "shifted_log"),
  n_boot = 200L, conf = .95, seed = 123L, workers = 1L,
  shift_bounds = NULL)
}
\arguments{
  \item{value}{Finite numeric routine results from one physiological partition and comparable assay.}
  \item{setting}{Exactly two service labels, one per routine result.}
  \item{anchor}{Independent clinically screened reference measurements; at least two values.}
  \item{transform}{Family selection passed to \code{fit_anchored_ri}. Defaults to automatic selection.}
  \item{n_boot}{Number of resamples, at least two. The default is 200; very small counts are only useful for software checks.}
  \item{conf}{Confidence level, strictly between zero and one.}
  \item{seed}{Nonnegative integer seed. The caller's random-number state is restored.}
  \item{workers}{Positive integer number of parallel Unix workers. Windows uses one worker and warns if more were requested.}
  \item{shift_bounds}{Optional fixed offset bounds passed to every fit. If NULL, the data-derived bounds are recalculated within each resample.}
}
\description{Resamples routine observations within each setting and screened observations independently. Every resample repeats component selection, family selection when requested, and offset estimation.}
\details{The point estimate comes from the full dataset. Percentile confidence intervals use successful bootstrap fits. Failed fits and selected numerical boundary fits are recorded. If fewer than 80 percent succeed, confidence limits are returned as missing with a warning. This threshold is a reporting rule, not evidence of valid coverage. With small screened samples, distribution mismatch or unstable model selection, percentile intervals may have poor coverage. Individuals must contribute one result each; this function does not implement clustered patient resampling. Parallel and sequential runs use the same per-replicate seeds.}
\value{A list containing \code{fit}, \code{limits}, a \code{ci} matrix with rows lower/upper and columns low/high, \code{confidence}, a \code{bootstrap} data frame, and the counts \code{requested}, \code{successful}, \code{failed} and \code{boundary_fits}. Each bootstrap row reports limits, selected family, offset, component count, boundary flag, status and error.}
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
