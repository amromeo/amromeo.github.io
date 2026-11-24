# SPDX-License-Identifier: GPL-3.0-only
# Multistart selection and original-scale SD constraints from Task E Arm B.
.log_sd_factor <- function(s) sqrt(exp(s^2) * expm1(s^2))

.fit_mixture_candidate <- function(d, family, k, models, rules, fit_seed) {
  shifted <- family == "shifted_log"
  linear <- family == "identity"
  coords <- .mixture_coordinates(d, family)
  box <- .mixture_bounds(d, family, k, coords)
  spec <- as.double(c(k, d$g, d$n, coords, d$location, d$spread,
                      if (linear) 1 else 0, if (shifted) NA_real_ else 0))
  starts <- .grid_starts(d, family, k, coords)
  if (k > 1L) starts <- c(starts, .nested_starts(models[[paste(family, k-1L)]]))
  if (shifted && d$zero)
    starts[[length(starts)+1L]] <- list(p = models[[paste("log", k)]]$p)

  # Historical drivers can source the R files after installing this checkout.
  native_fit <- get("C_anchored_fit_start", envir = asNamespace("anchoredRI"), inherits = FALSE)
  best <- NULL
  consider <- function(start, bounds) {
    z <- .Call(native_fit, d$x, d$counts, d$anchors, spec, as.double(start),
               as.double(bounds$lo), as.double(bounds$hi), as.double(rules), 1000L, 1L)
    use_start <- is.finite(z[[2]]) && (!is.finite(z[[3]]) || z[[2]] > z[[3]])
    ll <- if (use_start) z[[2]] else z[[3]]
    if (!is.finite(ll) || (!is.null(best) && ll <= best$ll)) return()
    par <- if (use_start) z[[7]] else z[[1]]
    p <- .unpack_mixture(d, family, k, par, coords)
    median <- if (linear) p$mu + 1 else p$shift + exp(p$mu)
    minimum <- pmax(rules[1], rules[2], pmax(median, 0) * rules[3])
    if (shifted) minimum <- pmax(minimum, exp(p$mu) * .log_sd_factor(.005))
    dispersion <- if (linear) p$sigma else exp(p$mu) * .log_sd_factor(p$sigma)
    if (!all(is.finite(c(par, minimum, dispersion))) ||
        any(dispersion < minimum * (1 - 1e-11))) return()
    at_bound <- any(abs(par - bounds$lo) < 1e-5 | abs(par - bounds$hi) < 1e-5) ||
      any(abs(dispersion - minimum) <= 1e-8 * pmax(1, minimum))
    best <<- list(ll = ll, par = par, box = bounds, p = p, at_bound = at_bound)
  }
  for (start in starts) {
    par <- .pack_mixture(d, family, start$p, coords)
    bounds <- box
    if (isTRUE(start$analytic)) bounds$hi <- pmax(bounds$hi, par)
    consider(par, bounds)
  }
  if (is.null(best))
    stop(sprintf("No feasible starting fit for %s with %d component(s).", family, k))
  order <- match(family, c("identity", "log", "shifted_log"))
  seed <- (as.double(fit_seed) + 1000000 * order + 100000 * k) %% (as.double(.Machine$integer.max) + 1)
  set.seed(as.integer(seed), kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  base <- best
  step <- c(.45 * exp(base$par[k + seq_len(k)]), rep(.35, k), rep(.65, d$g * (k-1L)),
            if (shifted) .12 * diff(d$bounds) / d$spread)
  for (j in seq_len(20L)) consider(base$par + rnorm(length(step)) * step, base$box)
  best$k <- k
  best$transform <- family
  best$bic <- -2 * best$ll + (2*k + d$g * (k-1L) + shifted) * log(d$n)
  best
}

.mixture_results <- function(models, d, transform, unit, fit_seed) {
  # Family order, then component count, breaks exact BIC ties as in Task E.
  scores <- vapply(models, function(m) m$bic, numeric(1))
  selected <- models[[which.min(scores)]]
  family <- selected$transform
  k <- selected$k
  candidates <- Filter(function(m) m$transform == family, models)
  central <- function(m) {
    y <- m$p$mu[1] + qnorm(c(.025, .975)) * m$p$sigma[1]
    if (family == "identity") y + 1 else m$p$shift + exp(y)
  }
  limits <- t(vapply(candidates, central, numeric(2)))
  dimnames(limits) <- list(as.character(seq_along(candidates)), c("lower", "upper"))
  shifts <- vapply(candidates, function(m) m$p$shift, numeric(1), USE.NAMES = FALSE)
  at_bound <- vapply(candidates, function(m) m$at_bound, logical(1), USE.NAMES = FALSE)
  bics <- vapply(candidates, function(m) m$bic, numeric(1), USE.NAMES = FALSE)
  names(bics) <- names(shifts) <- rownames(limits)
  p <- selected$p
  colnames(p$w) <- d$groups
  correction <- switch(family, shifted_log = 0,
    identity = 2 * d$n * log(d$spread),
    log = 2 * d$n * log(.routine_scale(log(d$value))) + 2 * sum(log(c(d$value, d$anchor))))
  result <- list(limits = limits[k, ], components = k, weights = p$w,
                 bic = bics - correction, bic_original = bics, selected_bic = unname(bics[k]),
                 candidate_limits = limits, candidate_at_bound = at_bound, candidate_shift = shifts,
                 transform = family, shift = unname(p$shift), at_bound = selected$at_bound,
                 shift_bounds = if (family == "shifted_log") d$bounds else NULL,
                 requested_transform = transform, rounding_unit = unit, fit_seed = fit_seed,
                 reference = list(mu = unname(p$mu[1] + as.integer(family == "identity")),
                                  sigma = unname(p$sigma[1])))
  if (family == "shifted_log") result$shift_at_bound <- any(abs(p$shift - d$bounds) / d$spread < 1e-5)
  if (transform == "auto") {
    result$selection <- do.call(rbind, lapply(unname(models), function(m)
      data.frame(transform = m$transform, components = m$k, shift = unname(m$p$shift),
                 bic_original = m$bic, at_bound = m$at_bound)))
  }
  result
}

.arm_b_fit <- function(value, setting, anchor, transform, shift_bounds, rounding_unit, fit_seed) {
  if (!is.numeric(value) || !is.numeric(anchor) || length(anchor) < 2L ||
      length(value) != length(setting) || !length(unique(setting)) %in% 1:2 ||
      anyNA(setting) || any(!is.finite(c(value, anchor))))
    stop("Supply numeric results from one or two settings and independent reference values.")
  if (transform == "log" && any(c(value, anchor) <= 0)) stop("Log requires positive values.")
  if (!is.null(shift_bounds) && !transform %in% c("auto", "shifted_log"))
    stop("shift_bounds applies only to shifted_log or auto.")
  check_integer(fit_seed, "fit_seed", 0L)
  d <- .mixture_data(value, setting, anchor, shift_bounds)
  unit <- .measurement_unit(d$x, rounding_unit)
  blog <- if (d$positive) max(.005, .05 * .routine_scale(log(value))) else .005
  rules <- c(unit, .05 * d$spread, .log_sd_factor(blog))
  families <- if (transform == "auto") c("identity", if (d$positive) "log", "shifted_log") else transform
  # A feasible zero shift uses the ordinary log fit as one shifted-log start.
  needed <- if (transform == "shifted_log" && d$zero) c("log", "shifted_log") else families
  restore <- preserve_rng()
  on.exit(restore())
  models <- list()
  for (k in 1:5) for (family in needed)
    models[[paste(family, k)]] <- .fit_mixture_candidate(d, family, k, models, rules, fit_seed)
  ordered <- unlist(lapply(families, function(family) paste(family, 1:5)), use.names = FALSE)
  .mixture_results(models[ordered], d, transform, unit, fit_seed)
}
