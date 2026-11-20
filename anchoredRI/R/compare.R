# The baseline is the published package, with bootstrap disabled for point estimates.
stock_limits <- function(x) {
  if (!requireNamespace("refineR", quietly = TRUE)) {
    stop("Install the optional package refineR to run its comparators.", call. = FALSE)
  }
  fit <- refineR::findRI(x, NBootstrap = 0L)
  limits <- as.numeric(refineR::getRI(fit)$PointEst)
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("refineR did not return two ordered finite reference limits.", call. = FALSE)
  }
  limits
}

interval_rows <- function(method, limits, components = NA_integer_,
                          at_bound = NA, enrichment_gap = NA_real_) {
  data.frame(method = method, lower = unname(limits[, 1]), upper = unname(limits[, 2]),
             components = components, at_bound = at_bound, enrichment_gap = enrichment_gap,
             row.names = NULL)
}

screened_comparisons <- function(value, setting, anchor) {
  fit <- fit_anchored_ri(value, setting, anchor)
  direct <- mean(anchor) + qnorm(c(.025, .975)) * sd(anchor)
  interval_rows(c("anchored", "anchored_two", "anchor_only"),
                rbind(fit$limits, fit$candidate_limits["2", ], direct),
                c(fit$components, 2L, NA_integer_),
                c(fit$at_bound, fit$candidate_at_bound[2L], NA))
}

compare_ri_methods <- function(value, setting, anchor = numeric(),
                               enriched_setting = "clinic", include_stock = TRUE) {
  grouped <- fit_setting_ri(value, setting, enriched_setting)
  pooled <- fit_setting_ri(value, setting, enriched_setting, mixing = "pooled")
  rows <- interval_rows(c("joint", "pooled_mixture"), rbind(grouped$limits, pooled$limits),
                        components = 2L, at_bound = c(grouped$at_bound, pooled$at_bound),
                        enrichment_gap = c(grouped$enrichment_gap, pooled$enrichment_gap))
  if (length(anchor)) {
    rows <- rbind(rows, screened_comparisons(value, setting, anchor))
  }
  if (include_stock) {
    limits <- rbind(stock_limits(value), stock_limits(value[setting == enriched_setting]))
    rows <- rbind(rows, interval_rows(c("refineR_all", "refineR_clinic"), limits))
  }
  rows
}

benchmark_setting_ri <- function(reps = 12L, scenarios = 1:4, include_stock = TRUE,
                                 n_clinic = 1000L, n_ward = 3000L, n_anchor = 40L) {
  check_integer(reps, "reps", 1L, 1000000L)
  if (!length(scenarios)) stop("Supply at least one scenario.", call. = FALSE)
  results <- list()
  for (s in scenarios) for (k in seq_len(reps)) {
    d <- simulate_ri_study(s, k, n_clinic, n_ward, n_anchor)
    out <- compare_ri_methods(d$routine$value, d$routine$setting, d$anchor,
                              include_stock = include_stock)
    out$scenario <- s
    out$replicate <- k
    out$lower_error <- out$lower - d$truth[1]
    out$upper_error <- out$upper - d$truth[2]
    results[[length(results) + 1L]] <- out
  }
  do.call(rbind, results)
}

summarize_ri_benchmark <- function(results) {
  summary_data <- data.frame(scenario = results$scenario, method = results$method,
                             lower_mae = abs(results$lower_error),
                             upper_mae = abs(results$upper_error))
  summary_data$worst_limit_mae <- pmax(summary_data$lower_mae, summary_data$upper_mae)
  aggregate(cbind(lower_mae, upper_mae, worst_limit_mae) ~ scenario + method,
            data = summary_data, FUN = mean)
}
