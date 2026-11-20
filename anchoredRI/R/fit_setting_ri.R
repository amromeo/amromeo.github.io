# SPDX-License-Identifier: GPL-3.0-only
fit_setting_ri <- function(value, setting, enriched_setting,
                           transform = c("identity", "log"),
                           mixing = c("by_setting", "pooled")) {
  transform <- match.arg(transform)
  mixing <- match.arg(mixing)
  if (!is.numeric(value) || length(value) != length(setting) ||
      any(!is.finite(value)) || anyNA(setting)) {
    stop("Supply one finite numeric result and one setting for every observation.")
  }
  setting <- as.character(setting)
  groups <- unique(setting)
  if (length(groups) != 2L || length(enriched_setting) != 1L ||
      !enriched_setting %in% groups) {
    stop("Use exactly two settings and name the setting expected to be enriched.")
  }
  if (transform == "log" && any(value <= 0)) stop("The log scale needs positive results.")
  groups <- c(enriched_setting, setdiff(groups, enriched_setting))
  group <- match(setting, groups)
  forward <- if (transform == "log") log else identity
  inverse <- if (transform == "log") exp else identity
  y <- forward(value)
  center <- median(y)
  scale <- mad(y, center = center)
  if (!is.finite(scale) || scale <= 0) scale <- sd(y)
  if (!is.finite(scale) || scale <= 0 || length(unique(y)) < 4L) {
    stop("Too few distinct results for a two-component fit.")
  }
  z <- (y - center) / scale
  support <- sort(unique(z))
  n_support <- length(support)
  counts <- matrix(tabulate(match(z, support) + (group - 1L) * n_support,
                            nbins = 2L * n_support), nrow = n_support,
                   dimnames = list(NULL, groups))
  fit_counts <- if (mixing == "by_setting") counts else cbind(rowSums(counts))
  n_group <- colSums(fit_counts)
  n <- sum(n_group)

  # Group identical reported values exactly for speed. Each distinct reported
  # value contributes its point density. Quantization is not explicitly modeled.
  calc <- function(par, gradient = FALSE) {
    mu <- par[1:2]
    sigma <- exp(par[3:4])
    eta <- par[-(1:4)]
    a <- outer(dnorm(support, mu[1], sigma[1], log = TRUE),
               plogis(eta, log.p = TRUE), "+")
    b <- outer(dnorm(support, mu[2], sigma[2], log = TRUE),
               plogis(eta, lower.tail = FALSE, log.p = TRUE), "+")
    top <- pmax(a, b)
    total <- top + log(exp(a - top) + exp(b - top))
    if (!gradient) return(-sum(fit_counts * total) / n)
    first <- fit_counts * exp(a - total)
    w1 <- rowSums(first)
    w2 <- rowSums(fit_counts) - w1
    u1 <- (support - mu[1]) / sigma[1]
    u2 <- (support - mu[2]) / sigma[2]
    score <- c(sum(w1 * u1) / sigma[1], sum(w2 * u2) / sigma[2],
               sum(w1 * (u1^2 - 1)), sum(w2 * (u2^2 - 1)),
               colSums(first) - n_group * plogis(eta))
    -score / n
  }

  # Multiple starts explore narrow and broad components on both sides.
  # The scale floor is a numerical regularizer in robust SD units.
  q <- quantile(z, c(.15, .35, .50, .65, .85), names = FALSE)
  starts <- rbind(c(q[1], q[4], .8, 1.0), c(q[2], q[5], .7, 1.0),
                  c(q[3], q[4], .55, 1.5), c(q[2], q[3], 1.5, .55),
                  c(median(z[group == 1L]), median(z[group == 2L]), .8, 1.2))
  eta0 <- if (mixing == "by_setting") qlogis(c(.65, .35)) else 0
  lo <- c(rep(min(z), 2), rep(log(.05), 2), rep(qlogis(.001), length(eta0)))
  hi <- c(rep(max(z), 2), rep(log(20), 2), rep(qlogis(.999), length(eta0)))
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    init <- c(starts[i, 1:2], log(starts[i, 3:4]), eta0)
    optim(init, calc, function(p) calc(p, TRUE), method = "L-BFGS-B",
          lower = lo, upper = hi, control = list(maxit = 400L, factr = 1e6))
  })
  ok <- vapply(fits, function(x) x$convergence == 0L && is.finite(x$value), logical(1))
  if (!any(ok)) stop("The likelihood optimization did not converge.")
  fits <- fits[ok]
  objective <- vapply(fits, function(x) x$value, numeric(1))
  best <- fits[[which.min(objective)]]
  par <- best$par
  mu <- par[1:2]
  sigma <- exp(par[3:4])
  eta <- rep(par[-(1:4)], length.out = 2L)
  log_ratio <- dnorm(support, mu[1], sigma[1], log = TRUE) -
               dnorm(support, mu[2], sigma[2], log = TRUE)
  membership <- plogis(outer(log_ratio, eta, "+"))
  fraction1 <- colSums(counts * membership) / colSums(counts)
  # In the pooled ablation setting labels enter only *after* the fit to name
  # the component. Use the same enrichment assumption for both models.
  candidate <- if (fraction1[1] >= fraction1[2]) 1L else 2L
  component_limits <- t(vapply(seq_len(2L), function(j) {
    inverse(center + scale * (mu[j] + sigma[j] * qnorm(c(.025, .5, .975))))
  }, numeric(3)))
  colnames(component_limits) <- c("lower", "median", "upper")
  fractions <- rbind(fraction1, 1 - fraction1)
  colnames(fractions) <- groups
  list(limits = component_limits[candidate, c("lower", "upper")],
       candidate = candidate,
       component_summary = data.frame(component = 1:2, component_limits,
                                      fractions, check.names = FALSE),
       enrichment_gap = abs(diff(unname(fraction1))),
       at_bound = any(abs(par - lo) < 1e-5 | abs(par - hi) < 1e-5),
       loglik_standardized = -n * best$value,
       converged_starts = length(fits))
}
