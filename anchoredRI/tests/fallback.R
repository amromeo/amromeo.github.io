library(anchoredRI)

# Use all 60 screened results for peak voting, blending and the exact null.
measure <- function(n) round(10 * exp(.15 * qnorm(ppoints(n))), 1)
fallback <- getFromNamespace(".fit_anchor_refine", "anchoredRI")
null_p <- getFromNamespace("stage1_null_p", "anchoredRI")
set.seed(31)
saved_rng <- .Random.seed
result <- fallback(measure(1000L), measure(60L), .1, c(.025, .975), 123L, 1L)
check <- result$check
stopifnot(identical(saved_rng, .Random.seed), identical(result$weight, .25),
          result$status == "reported", isTRUE(result$blend$applied),
          length(result$parent$election$assignments) == 60L,
          check$anchor_n == 60L, check$status == "keep",
          identical(check$p, null_p(check$statistic, 60L)),
          abs(check$p - null_p(check$statistic, 40L)) > 1e-6)
