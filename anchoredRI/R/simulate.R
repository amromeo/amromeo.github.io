# Simulated aPTT-like values in seconds at a single age and on one platform.
# Independent individuals. Target is the continuous pre-rounding interval.
# 1: small fractions, well separated. 2: heavily mixed, shared raised results.
# 3: ward pathology includes an extra low component. 4: equal fractions in both.
simulate_ri_case <- function(scenario, replicate, n_clinic = 1000L, n_ward = 3000L) {
  check_integer(scenario, "scenario", 1L, 4L)
  check_integer(replicate, "replicate", 1L, 1000000L)
  check_integer(n_clinic, "n_clinic", 1L)
  check_integer(n_ward, "n_ward", 1L)
  restore_rng <- preserve_rng()
  on.exit(restore_rng())
  set.seed(26000L + 100L * scenario + replicate)
  setting <- rep(c("clinic", "ward"), c(n_clinic, n_ward))
  n <- length(setting)
  p_clinic <- c(.05, .20, .20, .40)[scenario]
  p_ward <- c(.10, .75, .75, .40)[scenario]
  affected <- runif(n) < ifelse(setting == "clinic", p_clinic, p_ward)
  value <- rnorm(n, 30, 3)
  illness <- rnorm(n, 36, 4.5)
  if (scenario == 1L) illness <- rnorm(n, 42, 4)
  if (scenario == 3L) {
    low <- setting == "ward" & runif(n) < .55
    illness <- rnorm(n, ifelse(low, 24, 36), ifelse(low, 2.5, 4.5))
  }
  value[affected] <- illness[affected]
  data.frame(value = round(value, 1), setting = setting)
}

# Return a cleanup function for deterministic draws without changing the
# caller's subsequent random numbers.
preserve_rng <- function() {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }
}

check_integer <- function(x, label, lower, upper = .Machine$integer.max) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < lower || x > upper) {
    stop(sprintf("%s must be one integer between %d and %d.", label, lower, upper),
         call. = FALSE)
  }
}

simulate_ri_study <- function(scenario = 3L, replicate = 1L, n_clinic = 250L,
                              n_ward = 750L, n_anchor = 40L,
                              selection = c("clinical", "in_range")) {
  selection <- match.arg(selection)
  check_integer(n_anchor, "n_anchor", 0L)
  routine <- simulate_ri_case(scenario, replicate, n_clinic, n_ward)
  restore_rng <- preserve_rng()
  on.exit(restore_rng())
  set.seed(700000L + 100L * scenario + replicate)
  if (selection == "clinical") {
    anchor <- round(rnorm(n_anchor, 30, 3), 1)
  } else {
    # An intentionally biased screen also requires a result of 26 to 34 s.
    anchor <- numeric()
    while (length(anchor) < n_anchor) {
      draw <- round(rnorm(max(256L, n_anchor - length(anchor)), 30, 3), 1)
      anchor <- c(anchor, draw[draw >= 26 & draw <= 34])
    }
    anchor <- head(anchor, n_anchor)
  }
  truth <- setNames(30 + 3 * qnorm(c(.025, .975)), c("lower", "upper"))
  list(routine = routine, anchor = anchor, truth = truth)
}
