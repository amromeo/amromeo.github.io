# SPDX-License-Identifier: GPL-3.0-only
.fit_anchor_refine <- function(value, anchor, unit, probs, seed, pit_seed) {
  result <- list(status = "withheld", reason = "", weight = .25,
                 rounding_unit = unit, fit_seed = seed, pit_seed = pit_seed,
                 routine_n = sum(value >= 0), anchor_n = length(anchor),
                 parent = NULL, blend = NULL, check = NULL,
                 limits = stats::setNames(rep(NA_real_, length(probs)), .ri_limit_names(probs)))
  result$candidate_limits <- result$limits
  if (result$anchor_n < 40L) {
    result$status <- "skipped"
    result$reason <- "The anchorRefine fallback requires at least 40 screened patients."
    message(result$reason)
    return(result)
  }
  if (any(anchor < 0)) {
    result$reason <- "invalid_anchor_panel"
    return(result)
  }
  if (unit / 2 <= 0) {
    result$reason <- "invalid_rounding_grid"
    return(result)
  }
  if (result$routine_n <= 10L) {
    result$reason <- "insufficient_nonnegative_routine_results"
    return(result)
  }

  restore_rng <- preserve_rng()
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    restore_rng()
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  parent <- tryCatch(stage2_selected_fit(value, anchor, unit, seed), error = identity)
  if (inherits(parent, "error")) {
    result$reason <- conditionMessage(parent)
    return(result)
  }
  # The stored fit parameters suffice for the check and point limits.
  if (!is.null(parent$fit)) parent$fit$Data <- NULL
  result$parent <- parent
  if (!stage2_valid_fit(parent$fit)) {
    result$reason <- "invalid_parent_fit"
    return(result)
  }
  post <- tryCatch(stage3_post(parent, anchor, unit, t = result$weight), error = identity)
  if (inherits(post, "error")) {
    result$reason <- conditionMessage(post)
    return(result)
  }
  result$blend <- list(fit = post$fit, applied = identical(post$reason, "post_adjusted"),
                       reason = post$reason, estimate = post$post_estimate)

  checked <- tryCatch(stage1_check(anchor, post$fit, decimals = NULL, seed = pit_seed,
                                   increment = unit),
                      error = function(e) list(status = "check_failed", reason = conditionMessage(e),
                                               statistic = NA_real_, p = NA_real_))
  result$check <- checked
  model <- post$fit
  endpoints <- tryCatch(.refineR_200$cdfTruncatedBoxCox(model$Mu, model$Sigma, model$Lambda,
                                                      model$Shift, probs), error = identity)
  if (inherits(endpoints, "error")) {
    result$reason <- conditionMessage(endpoints)
    return(result)
  }
  if (length(endpoints) != length(probs)) {
    result$reason <- "unusable_blend_limits"
    return(result)
  }
  result$candidate_limits[] <- endpoints
  result$reason <- checked$reason
  if (checked$status == "decline") {
    result$status <- "declined"
  } else if (checked$status == "keep") {
    if (all(is.finite(endpoints)) &&
        (length(probs) == 1L || endpoints[1L] < endpoints[2L])) {
      result$limits <- result$candidate_limits
      result$status <- "reported"
    } else {
      result$reason <- "unusable_blend_limits"
    }
  }
  result
}
