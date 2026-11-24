# SPDX-License-Identifier: GPL-3.0-only
fit_anchored_ri <- function(value, setting, anchor,
                            transform = c("auto", "identity", "log", "shifted_log"),
                            shift_bounds = NULL, safeguard = TRUE, probs = c(.025, .975),
                            rounding_unit = NULL, fit_seed = 1L) {
  transform <- match.arg(transform)
  .validate_ri_reporting(safeguard, probs)
  fit <- .arm_b_fit(value, setting, anchor, transform, shift_bounds, rounding_unit, fit_seed)
  .apply_ri_safeguard(fit, value, anchor, safeguard, probs)
}
