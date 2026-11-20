# SPDX-License-Identifier: GPL-3.0-only
# Joint likelihood for routine results plus independently screened reference
# observations. The latter enter only component 1. Their clinical selection must
# be independent of the target result. Use separate individuals for validation.
# BIC selects one, two or three components. All settings share the component
# shapes, but each setting has its own mixture proportions.
fit_anchored_ri <- function(value, setting, anchor, transform = c("identity", "log")) {
  transform <- match.arg(transform)
  if (!is.numeric(value) || !is.numeric(anchor) || length(anchor) < 2L ||
      length(value) != length(setting) || length(unique(setting)) != 2L ||
      anyNA(setting) || any(!is.finite(c(value, anchor)))) {
    stop("Supply numeric results from two settings and independent reference values.")
  }
  if (transform == "log" && any(c(value, anchor) <= 0)) stop("Log requires positive values.")
  forward <- if (transform == "log") log else identity
  inverse <- if (transform == "log") exp else identity
  y <- forward(value)
  center <- median(y)
  scale <- mad(y)
  if (scale <= 0) scale <- sd(y)
  if (!is.finite(scale) || scale <= 0 || length(unique(y)) < 4L) stop("Insufficient variation.")
  z <- (y - center) / scale
  a <- (forward(anchor) - center) / scale
  x <- sort(unique(z))
  groups <- unique(as.character(setting))
  g <- match(as.character(setting), groups)
  counts <- matrix(tabulate(match(z, x) + (g - 1L) * length(x), 2L * length(x)),
                   ncol = 2L, dimnames = list(NULL, groups))
  n_group <- colSums(counts)
  n <- length(z) + length(a)
  # Component means, log SDs, and non-baseline log odds for each setting.
  fit_k <- function(k) {
    unpack <- function(par) {
      eta <- rbind(matrix(tail(par, 2L * (k - 1L)), ncol = 2L), c(0, 0))
      weights <- exp(sweep(eta, 2, apply(eta, 2, max), "-"))
      list(mu = par[seq_len(k)], sigma = exp(par[k + seq_len(k)]),
           w = sweep(weights, 2, colSums(weights), "/"))
    }
    calc <- function(par, gradient = FALSE) {
      p <- unpack(par)
      log_density <- vapply(seq_len(k), function(j) {
        dnorm(x, p$mu[j], p$sigma[j], log = TRUE)
      }, numeric(length(x)))
      ll <- sum(dnorm(a, p$mu[1], p$sigma[1], log = TRUE))
      contributions <- matrix(0, length(x), k)
      weight_score <- matrix(0, k - 1L, 2L)
      for (h in 1:2) {
        parts <- sweep(log_density, 2, log(p$w[, h]), "+")
        peak <- do.call(pmax, as.data.frame(parts))
        scaled <- exp(parts - peak)
        denom <- rowSums(scaled)
        ll <- ll + sum(counts[, h] * (peak + log(denom)))
        if (gradient) {
          mass <- scaled * (counts[, h] / denom)
          contributions <- contributions + mass
          weight_score[, h] <- (colSums(mass) - n_group[h] * p$w[, h])[-k]
        }
      }
      if (!gradient) return(-ll / n)
      u <- sweep(outer(x, p$mu, "-"), 2, p$sigma, "/")
      mean_score <- colSums(contributions * u) / p$sigma
      scale_score <- colSums(contributions * (u^2 - 1))
      ua <- (a - p$mu[1]) / p$sigma[1]
      mean_score[1] <- mean_score[1] + sum(ua) / p$sigma[1]
      scale_score[1] <- scale_score[1] + sum(ua^2 - 1)
      -c(mean_score, scale_score, weight_score) / n
    }
    slots <- if (k == 2L) list(.15, .5, .85) else list(c(.12,.85), c(.60,.90), c(.10,.40))
    ref_sd <- min(19.9, max(.06, sd(a)))
    lo <- c(rep(min(c(z,a)), k), rep(log(.05), k), rep(-7, 2L * (k-1L)))
    hi <- c(rep(max(c(z,a)), k), rep(log(20), k), rep(7, 2L * (k-1L)))
    fits <- lapply(slots, function(probs) {
      init <- c(mean(a), unname(quantile(z, probs)), log(ref_sd),
                rep(log(1), k-1L), rep(0, 2L * (k-1L)))
      optim(init, calc, function(p) calc(p, TRUE), method = "L-BFGS-B",
            lower = lo, upper = hi, control = list(maxit = 450L, factr = 1e6))
    })
    fits <- Filter(function(f) f$convergence == 0L && is.finite(f$value), fits)
    if (!length(fits)) return(NULL)
    fit <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
    list(p = unpack(fit$par), bic = 2*n*fit$value + length(fit$par)*log(n),
         at_bound = any(abs(fit$par-lo)<1e-5 | abs(fit$par-hi)<1e-5))
  }
  all_z <- c(z, a)
  one_mu <- mean(all_z)
  one_sd <- sqrt(mean((all_z - one_mu)^2))
  one <- list(p = list(mu = one_mu, sigma = one_sd, w = matrix(1, 1L, 2L)),
              bic = -2 * sum(dnorm(all_z, one_mu, one_sd, log = TRUE)) + 2*log(n),
              at_bound = FALSE)
  fits <- list(one, fit_k(2L), fit_k(3L))
  bics <- vapply(fits, function(f) if (is.null(f)) Inf else f$bic, numeric(1))
  k <- which.min(bics)
  p <- fits[[k]]$p
  colnames(p$w) <- groups
  candidate_limits <- t(vapply(fits, function(f) {
    if (is.null(f)) return(c(NA_real_, NA_real_))
    inverse(center + scale * (f$p$mu[1] + qnorm(c(.025, .975)) * f$p$sigma[1]))
  }, numeric(2)))
  dimnames(candidate_limits) <- list(as.character(1:3), c("lower", "upper"))
  limits <- inverse(center + scale * (p$mu[1] + qnorm(c(.025,.975))*p$sigma[1]))
  list(limits = setNames(limits, c("lower","upper")), components = k,
       weights = p$w, bic = setNames(bics, 1:3), at_bound = fits[[k]]$at_bound,
       candidate_limits = candidate_limits,
       candidate_at_bound = vapply(fits, function(f) {
         if (is.null(f)) NA else f$at_bound
       }, logical(1)))
}
