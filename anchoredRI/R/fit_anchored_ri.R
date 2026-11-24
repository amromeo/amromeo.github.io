# SPDX-License-Identifier: GPL-3.0-only
fit_anchored_ri <- function(value, setting, anchor,
                            transform = c("auto", "identity", "log", "shifted_log"),
                            shift_bounds = NULL, rounding_unit = NULL, fit_seed = 1L) {
  transform <- match.arg(transform)
  .arm_b_fit(value, setting, anchor, transform, shift_bounds, rounding_unit, fit_seed)
}
