# SPDX-License-Identifier: GPL-3.0-only
# Study adapter for anchoredRI 0.4.0, pinned at fe8105a in PROTOCOL.md.
# Native likelihood, constraints, gradients, starts, ranking and safeguards come
# from the installed package. Its namespace and on-disk sources are never edited.

h_weighted_candidate <- function(ns, weight) {
  candidate <- ns$.fit_mixture_candidate
  if (weight == 1) return(candidate)

  scope <- new.env(parent = ns)
  scope$.grid_starts <- function(d, family, k, coords) {
    starts <- ns$.grid_starts(d, family, k, coords)
    if (k == 1L && family != "shifted_log") {
      y <- if (family == "identity") d$value - 1 else log(d$value)
      a <- if (family == "identity") d$anchor - 1 else log(d$anchor)
      mu <- (weight * sum(y) + sum(a)) / d$n
      starts[[1L]]$p$mu <- mu
      starts[[1L]]$p$sigma <- sqrt((weight * sum((y - mu)^2) + sum((a - mu)^2)) / d$n)
    }
    starts
  }
  environment(candidate) <- scope
  candidate
}

# The package's anchor-only diagnostic profiles a 1001-point shift grid.
# It depends only on the identical panel shared by the two caps and 100 draws.
h_cached_anchors <- local({
  previous <- NULL
  result <- NULL
  function(anchor) {
    if (!identical(anchor, previous)) {
      result <<- asNamespace("anchoredRI")$.anchors_auto(anchor)
      previous <<- anchor
    }
    result
  }
})

h_missing_fit <- function(error, value, probs, rounding_unit, fit_seed) {
  ns <- asNamespace("anchoredRI")
  reason <- conditionMessage(error)
  limits <- setNames(rep(NA_real_, length(probs)), ns$.ri_limit_names(probs))
  check <- ns$.ri_flags(limits, NA_real_, min(value), reason)
  list(limits = limits, probs = probs, status = "withheld", anchoredRI_status = "withheld",
       reason = reason, anchoredRI_reason = reason, method = NA_character_,
       declined = FALSE, technical_withholding = TRUE, rounding_unit = rounding_unit,
       fit_seed = fit_seed, recommendation = ns$.ri_screening_advice,
       diagnostics = list(anchored = list(error = reason), anchorRefine = NULL),
       safeguard = c(list(enabled = TRUE), check[setdiff(names(check), "limits")]))
}

h_mixture <- function(value, anchor, rounding_unit, probs = c(.025, .975),
                      caps = c(3L, 5L), routine_weight = 1, fit_seed = 1L,
                      fallback = FALSE, fallback_seed = 123L, pit_seed = 1L,
                      keep_candidates = FALSE) {
  ns <- asNamespace("anchoredRI")
  if (!length(caps) || anyNA(caps) || any(!caps %in% 1:5) || anyDuplicated(caps))
    stop("Provide distinct integer component caps between 1 and 5.")
  if (length(routine_weight) != 1L || !is.finite(routine_weight) || routine_weight <= 0)
    stop("Routine weight must be one positive finite number.")
  if (isTRUE(fallback) && (routine_weight != 1 || !5L %in% caps))
    stop("Fallback is available only for the unweighted standard full fit.")
  caps <- as.integer(caps)
  total <- length(value) * routine_weight + length(anchor)
  # Native Objective::n is an int. Remove only floating-point roundoff in the
  # study's integer total (5040 for W, observed n+40 for unweighted fits).
  if (abs(total - round(total)) > 8 * .Machine$double.eps * max(1, total) ||
      total > .Machine$integer.max) stop("The native objective requires an integral total weight.")
  bic_n <- as.integer(round(total))
  ns$.validate_ri_reporting(TRUE, probs)
  setting <- rep("pooled", length(value))
  d <- tryCatch(ns$.arm_b_data(value, setting, anchor, "auto", NULL, rounding_unit, fit_seed),
                error = identity)
  failed <- inherits(d, "error")
  models <- list()
  errors <- setNames(rep("", length(caps)), caps)
  fail_at <- Inf
  failure <- NULL

  if (failed) {
    fail_at <- 1L
    failure <- d
  } else {
    if (routine_weight != 1) {
      # The native routine uses these double counts for mixture responsibilities,
      # the observation Jacobian and every derivative. Anchor counts stay fixed.
      d$counts <- d$counts * routine_weight
      d$n <- bic_n
    }
    blog <- if (d$positive) max(.005, .05 * ns$.routine_scale(log(d$value))) else .005
    rules <- c(d$unit, .05 * d$spread, ns$.log_sd_factor(blog))
    families <- c("identity", if (d$positive) "log", "shifted_log")
    candidate <- h_weighted_candidate(ns, routine_weight)
    restore <- ns$preserve_rng()
    on.exit(restore(), add = TRUE)
    for (k in seq_len(max(caps))) {
      for (family in families) {
        model <- tryCatch(candidate(d, family, k, models, rules, fit_seed), error = identity)
        if (inherits(model, "error")) {
          fail_at <- k
          failure <- model
          break
        }
        models[[paste(family, k)]] <- model
      }
      if (k == fail_at) break
    }
  }

  fits <- setNames(vector("list", length(caps)), caps)
  # Run the pinned reporting and fallback bodies in a private environment. The
  # shared candidates have already been optimized, so supply the ranked result.
  scope <- new.env(parent = ns)
  scope$.arm_b_data <- function(...) d
  scope$.arm_b_optimize <- function(...) {
    if (inherits(raw, "error")) stop(raw)
    raw
  }
  scope$.anchors_auto <- h_cached_anchors
  safeguard <- ns$.apply_ri_safeguard
  environment(safeguard) <- scope
  scope$.apply_ri_safeguard <- safeguard
  report <- ns$fit_anchored_ri
  environment(report) <- scope
  for (cap in caps) {
    key <- as.character(cap)
    raw <- if (cap >= fail_at) failure else tryCatch({
      ordered <- unlist(lapply(families, function(f) paste(f, seq_len(cap))), use.names = FALSE)
      fit <- ns$.mixture_results(models[ordered], d, "auto", d$unit, fit_seed)
      # The optional working-scale BIC in .mixture_results also needs the
      # weighted log Jacobian; original-scale ranking already used model$bic.
      if (routine_weight != 1 && fit$transform == "log")
        fit$bic <- fit$bic - 2 * (routine_weight - 1) * sum(log(d$value))
      fit
    }, error = identity)
    if (inherits(raw, "error")) errors[key] <- conditionMessage(raw)
    result <- if (failed) raw else tryCatch(suppressMessages(
      report(value, setting, anchor, probs = probs, rounding_unit = rounding_unit,
             fit_seed = fit_seed, fallback = fallback && cap == 5L, fallback_seed = fallback_seed,
             pit_seed = pit_seed)), error = identity)
    if (inherits(result, "error")) {
      errors[key] <- conditionMessage(result)
      result <- h_missing_fit(result, value, probs, rounding_unit, fit_seed)
    }
    fits[[key]] <- result
  }
  list(fits = fits, candidates = if (keep_candidates) models else NULL,
       errors = errors, bic_n = bic_n)
}

h_fit_row <- function(fit, bic_n = NA_real_) {
  # The public result masks declined limits. Raw requested limits survive in
  # diagnostics$anchored even when the final package used the fallback.
  first <- fit$diagnostics$anchored
  check <- fit$safeguard
  fallback <- fit$diagnostics$anchorRefine
  num <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_
  flag <- function(x) if (length(x) == 1L) as.logical(x) else NA
  label <- function(x) if (length(x) == 1L) as.character(x) else NA_character_
  reason <- function(x) paste(x, collapse = " | ")
  endpoint <- function(x, side) num(x[side])
  weights <- fit$weights
  valid_weights <- is.matrix(weights) && ncol(weights) == 1L && nrow(weights) >= 1L &&
    all(is.finite(weights)) && all(weights >= 0) && abs(sum(weights) - 1) <= 1e-10
  k <- as.integer(num(fit$components))
  parameters <- 2L * k + (k - 1L) + as.integer(identical(fit$transform, "shifted_log"))
  bic <- num(fit$selected_bic)
  data.frame(raw_lower = endpoint(first$limits, "lower"), raw_upper = endpoint(first$limits, "upper"),
             lower = endpoint(fit$limits, "lower"), upper = endpoint(fit$limits, "upper"),
             status = label(fit$status), first_status = label(fit$anchoredRI_status),
             method = label(fit$method), reason = reason(fit$reason),
             first_reason = reason(fit$anchoredRI_reason),
             fallback_status = label(fallback$status), fallback_reason = reason(fallback$reason),
             share = if (valid_weights) as.numeric(weights[1L, 1L]) else NA_real_,
             weight_sum = if (is.null(weights)) NA_real_ else sum(weights),
             share_reason = if (valid_weights) "" else if (length(first$error)) reason(first$error)
                            else "missing or invalid routine weights",
             components = k, transform = label(fit$transform), shift = num(fit$shift),
             mu = num(fit$reference$mu), sigma = num(fit$reference$sigma),
             bic = bic, bic_n = bic_n, loglik = (parameters * log(bic_n) - bic) / 2,
             p = num(check$lrt_p), lrt_statistic = num(check$lrt_statistic),
             lrt_flag = flag(check$lrt_flag), lower_flag = flag(check$lower_flag),
             lower_check_skipped = flag(check$lower_check_skipped), routine_min = num(check$routine_min),
             fit_error = reason(first$error), check_error = reason(check$error),
             optimizer_status = label(fit$optimizer_status),
             optimizer_code = as.integer(num(fit$optimizer_code)),
             stringsAsFactors = FALSE)
}
