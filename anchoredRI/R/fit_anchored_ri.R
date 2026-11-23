# SPDX-License-Identifier: GPL-3.0-only
# Joint likelihood for routine results plus independently screened reference
# observations. The latter enter only component 1. Their clinical selection must
# be independent of the target result. Use separate individuals for validation.
# BIC selects one, two or three components. All settings share the component
# shapes, but each setting has its own mixture proportions.
.fit_anchored_fixed <- function(value, setting, anchor, transform = c("identity", "log")) {
  transform <- match.arg(transform)
  if (!is.numeric(value) || !is.numeric(anchor) || length(anchor) < 2L ||
      length(value) != length(setting) || !length(unique(setting)) %in% 1:2 ||
      anyNA(setting) || any(!is.finite(c(value, anchor)))) {
    stop("Supply numeric results from one or two settings and independent reference values.")
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
  G <- length(groups)
  g <- match(as.character(setting), groups)
  counts <- matrix(tabulate(match(z, x) + (g - 1L) * length(x), G * length(x)),
                   ncol = G, dimnames = list(NULL, groups))
  n_group <- colSums(counts)
  n <- length(z) + length(a)
  # Component means, log SDs, and non-baseline log odds for each setting.
  fit_k <- function(k) {
    unpack <- function(par) {
      eta <- rbind(matrix(tail(par, G * (k - 1L)), ncol = G), rep(0, G))
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
      weight_score <- matrix(0, k - 1L, G)
      for (h in seq_len(G)) {
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
    lo <- c(rep(min(c(z,a)), k), rep(log(.05), k), rep(-7, G * (k-1L)))
    hi <- c(rep(max(c(z,a)), k), rep(log(20), k), rep(7, G * (k-1L)))
    fits <- lapply(slots, function(probs) {
      init <- c(mean(a), unname(quantile(z, probs)), log(ref_sd),
                rep(log(1), k-1L), rep(0, G * (k-1L)))
      tryCatch(optim(init, calc, function(p) calc(p, TRUE), method = "L-BFGS-B",
            lower = lo, upper = hi, control = list(maxit = 450L, factr = 1e6)),
            error = function(e) NULL)
    })
    fits <- Filter(function(f) !is.null(f) && f$convergence == 0L && is.finite(f$value), fits)
    if (!length(fits)) return(NULL)
    fit <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
    list(p = unpack(fit$par), bic = 2*n*fit$value + length(fit$par)*log(n),
         at_bound = any(abs(fit$par-lo)<1e-5 | abs(fit$par-hi)<1e-5))
  }
  all_z <- c(z, a)
  one_mu <- mean(all_z)
  one_sd <- sqrt(mean((all_z - one_mu)^2))
  one <- list(p = list(mu = one_mu, sigma = one_sd, w = matrix(1, 1L, G)),
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
  # Density conversion back to the measurement scale. This constant cancels
  # across component counts, but is essential when comparing transformations.
  correction <- 2*n*log(scale) + if (transform == "log") 2*sum(log(c(value, anchor))) else 0
  list(limits = setNames(limits, c("lower","upper")), components = k,
       weights = p$w, bic = setNames(bics, 1:3), at_bound = fits[[k]]$at_bound,
       bic_original = setNames(bics + correction, 1:3),
       selected_bic = bics[k] + correction, transform = transform, shift = 0,
       candidate_shift = setNames(rep(0, 3), 1:3), shift_bounds = NULL,
       candidate_limits = candidate_limits,
       candidate_at_bound = vapply(fits, function(f) {
         if (is.null(f)) NA else f$at_bound
       }, logical(1)))
}

# Preserve the existing default; the new automatic selection is opt-in.
fit_anchored_ri <- function(value, setting, anchor,
                            transform = c("identity", "log", "shifted_log", "auto"),
                            shift_bounds = NULL, safeguard = TRUE,
                            probs = c(.025, .975)) {
  transform <- match.arg(transform)
  .validate_ri_reporting(safeguard, probs)
  if (!is.numeric(value) || !is.numeric(anchor) || length(anchor) < 2L ||
      length(value) != length(setting) || !length(unique(setting)) %in% 1:2 ||
      anyNA(setting) || any(!is.finite(c(value, anchor)))) {
    stop("Supply numeric results from one or two settings and independent reference values.")
  }
  if (!is.null(shift_bounds) && !transform %in% c("shifted_log", "auto"))
    stop("shift_bounds applies only to shifted_log or auto.")
  if (transform %in% c("identity", "log")) {
    result <- .fit_anchored_fixed(value, setting, anchor, transform)
  } else if (transform == "shifted_log") {
    result <- .fit_anchored_shifted(value, setting, anchor, shift_bounds)
    if (!is.finite(result$selected_bic)) stop("No shifted-log candidate model converged.")
  } else {
    families <- c("identity", if (all(c(value, anchor) > 0)) "log", "shifted_log")
    fits <- lapply(families, function(family) {
      if (family == "shifted_log") .fit_anchored_shifted(value, setting, anchor, shift_bounds)
      else .fit_anchored_fixed(value, setting, anchor, family)
    })
    scores <- vapply(fits, function(f) f$selected_bic, numeric(1))
    if (!any(is.finite(scores))) stop("No candidate model converged.")
    result <- fits[[which.min(scores)]]
    result$selection <- do.call(rbind, lapply(seq_along(fits), function(i) {
      f <- fits[[i]]
      data.frame(transform = families[i], components = 1:3,
                 shift = unname(f$candidate_shift),
                 bic_original = unname(f$bic_original),
                 at_bound = unname(f$candidate_at_bound))
    }))
  }
  result$requested_transform <- transform
  .apply_ri_safeguard(result, value, anchor, safeguard, probs)
}
