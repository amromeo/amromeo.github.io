# SPDX-License-Identifier: GPL-3.0-only
fit_anchored_ri <- function(value, setting, anchor,
                            transform = c("auto", "identity", "log", "shifted_log"),
                            shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975),
                            rounding_unit = NULL, fit_seed = 1L) {
  transform <- match.arg(transform)
  .validate_ri_reporting(safeguard, probs)
  if (length(value) > 5000L)
    message("RIbench accuracy between 5,000 and 50,000 routine results has not been tested. ",
            "At 50,000, refineR was more accurate. ", .ri_large_dataset_advice)
  result <- .fit_anchored_interval(value, setting, anchor, transform, shift_bounds, safeguard,
                                   probs, rounding_unit, fit_seed)
  .signal_ri_fit(result, inferred_unit = is.null(rounding_unit))
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
