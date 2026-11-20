library(anchoredRI)
set.seed(17)
rng_before <- .Random.seed
d <- simulate_ri_study()
with_more <- simulate_ri_study(n_anchor = 80)
restricted <- simulate_ri_study(selection = "in_range")
stopifnot(identical(rng_before, .Random.seed),
          identical(d$routine, with_more$routine),
          identical(d$routine, restricted$routine),
          identical(d$anchor, with_more$anchor[1:40]),
          all(restricted$anchor >= 26 & restricted$anchor <= 34),
          length(d$anchor) == 40L, nrow(d$routine) == 1000L)
x <- d$routine
fit <- fit_anchored_ri(x$value, x$setting, d$anchor)
stopifnot(max(abs(fit$limits - d$truth)) < 3,
          isTRUE(all.equal(colSums(fit$weights), c(clinic = 1, ward = 1))),
          identical(unname(fit$limits), unname(fit$candidate_limits[fit$components, ])))
on_log <- fit_anchored_ri(exp(x$value/10), x$setting, exp(d$anchor/10), "log")
stopifnot(isTRUE(all.equal(unname(on_log$limits), unname(exp(fit$limits/10)), tolerance = 1e-4)))
fails <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
stopifnot(fails(fit_anchored_ri(x$value, rep("a", nrow(x)), d$anchor)),
          fails(fit_anchored_ri(x$value, x$setting, numeric())),
          fails(fit_anchored_ri(x$value, x$setting, -d$anchor, "log")),
          fails(simulate_ri_case(1.5, 1)), fails(simulate_ri_study(n_anchor = -1)))
two <- compare_ri_methods(x$value, x$setting, include_stock = FALSE)
stopifnot(identical(two$method, c("joint", "pooled_mixture")))
bench <- benchmark_setting_ri(1, 2, FALSE, 120, 360, 20)
summary <- summarize_ri_benchmark(bench)
stopifnot(nrow(bench) == 5, nrow(summary) == 5,
          all(summary$worst_limit_mae >= summary$lower_mae))
if (requireNamespace("refineR", quietly = TRUE)) {
  both <- compare_ri_methods(x$value, x$setting, include_stock = TRUE)
  original <- refineR::getRI(refineR::findRI(x$value, NBootstrap = 0L))$PointEst
  row <- both[both$method == "refineR_all", c("lower", "upper")]
  stopifnot(isTRUE(all.equal(as.numeric(row), as.numeric(original))))
}
