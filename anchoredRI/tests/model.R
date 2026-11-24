library(anchoredRI)
# Exercise raw model diagnostics. The reporting contract is tested in safeguard.R.
fit_anchored_ri <- function(...) anchoredRI::fit_anchored_ri(..., safeguard=FALSE)
bootstrap_anchored_ri <- function(...) anchoredRI::bootstrap_anchored_ri(..., safeguard=FALSE)

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
# The single normal component has an analytic fit. The two families have
# different physical SD floors and multistart paths for higher component counts.
identity_fit <- fit_anchored_ri(x$value, x$setting, d$anchor, "identity")
stopifnot(fit$requested_transform == "auto", fit$rounding_unit == .1,
          isTRUE(all.equal(on_log$candidate_limits["1", ],
                           exp(identity_fit$candidate_limits["1", ] / 10), tolerance = 1e-7)),
          max(abs(on_log$limits - exp(on_log$reference$mu +
                    qnorm(c(.025, .975)) * on_log$reference$sigma))) < 1e-10)
fails <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
stopifnot(fails(fit_anchored_ri(x$value, rep(c("a","b","c"), length.out=nrow(x)), d$anchor)),
          fails(fit_anchored_ri(x$value, x$setting, d$anchor, fit_seed = -1L)),
          fails(fit_anchored_ri(x$value, x$setting, d$anchor, rounding_unit = 1)),
          fails(fit_anchored_ri(x$value, x$setting, numeric())),
          fails(fit_anchored_ri(x$value, x$setting, -d$anchor, "log")),
          fails(simulate_ri_case(1.5, 1)), fails(simulate_ri_study(n_anchor = -1)))
two <- compare_ri_methods(x$value, x$setting, include_stock = FALSE)
stopifnot(identical(two$method, c("joint", "pooled_mixture")))
bench <- benchmark_setting_ri(1, 2, FALSE, 120, 360, 20)
summary <- summarize_ri_benchmark(bench)
stopifnot(nrow(bench) == 5, nrow(summary) == 5,
          all(summary$worst_limit_mae >= summary$lower_mae,na.rm=TRUE),
          sum(summary$attempted)==5,all(summary$reported<=summary$attempted))
if (requireNamespace("refineR", quietly = TRUE)) {
  both <- compare_ri_methods(x$value, x$setting, include_stock = TRUE)
  original <- refineR::getRI(refineR::findRI(x$value, NBootstrap = 0L))$PointEst
  row <- both[both$method == "refineR_all", c("lower", "upper")]
  stopifnot(isTRUE(all.equal(as.numeric(row), as.numeric(original))))
}
