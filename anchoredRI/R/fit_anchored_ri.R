# SPDX-License-Identifier: GPL-3.0-only
fit_anchored_ri <- function(value, setting, anchor,
                            transform = c("auto", "identity", "log", "shifted_log"),
                            shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975),
                            rounding_unit = NULL, fit_seed = 1L, fallback = TRUE,
                            fallback_seed = 123L, pit_seed = 1L) {
  transform <- match.arg(transform)
  .validate_ri_reporting(safeguard, probs)
  if (!is.logical(fallback) || length(fallback) != 1L || is.na(fallback))
    stop("fallback must be TRUE or FALSE.", call. = FALSE)
  check_integer(fallback_seed, "fallback_seed", 0L)
  check_integer(pit_seed, "pit_seed", 0L)
  if (length(value) > 5000L)
    message("RIbench accuracy between 5,000 and 50,000 routine results has not been tested directly.")

  data <- .arm_b_data(value, setting, anchor, transform, shift_bounds, rounding_unit, fit_seed)
  try_fallback <- safeguard && fallback
  fit <- if (try_fallback) {
    tryCatch(.arm_b_optimize(data, transform, fit_seed), error = identity)
  } else {
    .arm_b_optimize(data, transform, fit_seed)
  }
  if (inherits(fit, "error")) {
    reason <- conditionMessage(fit)
    limits <- stats::setNames(rep(NA_real_, length(probs)), .ri_limit_names(probs))
    check <- .ri_flags(limits, NA_real_, min(value), reason)
    result <- list(limits = limits, status = "withheld", reason = reason,
                   declined = FALSE, technical_withholding = TRUE, probs = probs,
                   rounding_unit = data$unit, fit_seed = fit_seed,
                   recommendation = .ri_screening_advice,
                   diagnostics = list(anchored = list(error = reason)),
                   safeguard = c(list(enabled = TRUE), check[setdiff(names(check), "limits")]))
  } else {
    result <- .apply_ri_safeguard(fit, value, anchor, safeguard, probs)
  }
  .signal_ri_fit(result, inferred_unit = is.null(rounding_unit))
  result$anchoredRI_status <- result$status
  result$anchoredRI_reason <- result$reason
  result$method <- if (result$status == "reported") "anchoredRI" else NA_character_
  result$diagnostics["anchorRefine"] <- list(NULL)
  if (try_fallback && result$status != "reported") {
    second <- .fit_anchor_refine(value, anchor, data$unit, probs, fallback_seed, pit_seed)
    result$diagnostics$anchorRefine <- second
    if (second$status == "reported") {
      result$limits <- second$limits
      result$method <- "anchorRefine"
      result$status <- "reported"
      result$reason <- character()
      result$declined <- result$technical_withholding <- FALSE
      result$recommendation <- ""
    }
  }
  result
}

.fit_anchored_interval <- function(value, setting, anchor, transform, shift_bounds, safeguard,
                                   probs, rounding_unit, fit_seed) {
  fit <- .arm_b_fit(value, setting, anchor, transform, shift_bounds, rounding_unit, fit_seed)
  .apply_ri_safeguard(fit, value, anchor, safeguard, probs)
}

.signal_ri_fit <- function(fit, inferred_unit) {
  if (inferred_unit) {
    unit <- format(fit$rounding_unit, digits = 2L, scientific = FALSE, trim = TRUE)
    message("Inferred rounding unit is ", unit,
            if (fit$rounding_unit == 0) " (no shared increment found)." else ".")
  }
  if (fit$status == "reported" && identical(fit$optimizer_status, "iteration_limit"))
    warning("The optimizer reached its iteration limit for the reported fit.", call. = FALSE)
  invisible(NULL)
}
