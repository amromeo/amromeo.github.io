# Run from anchoredRI: Rscript tools/build-guided-benchmark.R [workers] [reps] [bootstrap]
# New synthetic benchmark following Ammer et al. 2023 and the refineR vignette.
# This is not a reproduction of the six-analyte 2021 paper.
stopifnot(file.exists("DESCRIPTION"), requireNamespace("refineR", quietly = TRUE))
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(path)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args) >= 1L) as.integer(args[1]) else 1L
reps <- if (length(args) >= 2L) as.integer(args[2]) else 1L
B <- if (length(args) >= 3L) as.integer(args[3]) else 200L
stopifnot(workers >= 1L, reps >= 1L, B >= 200L)
dir.create("assets/guided", recursive = TRUE, showWarnings = FALSE)
dir.create("tools/guided-cache", recursive = TRUE, showWarnings = FALSE)

profiles <- data.frame(
  id = c("normal", "overlap", "lognormal", "shifted"),
  title = c("Normal, separated disease", "Normal, overlapping disease",
            "Right-skewed", "Right-skewed and shifted"),
  p_a = c(.05, .20, .05, .05), p_b = c(.15, .30, .15, .15),
  refine_model = c("BoxCox", "BoxCox", "BoxCox", "modBoxCox"),
  anchor_transform = c("identity", "identity", "log", "log"))

reference_quantile <- function(p, profile) {
  switch(profile, normal = qnorm(p, 30, 3), overlap = qnorm(p, 30, 3),
         lognormal = qlnorm(p, log(20), .35), shifted = 10 + qlnorm(p, log(15), .5))
}
reference_cdf <- function(x, profile) {
  switch(profile, normal = pnorm(x, 30, 3), overlap = pnorm(x, 30, 3),
         lognormal = plnorm(x, log(20), .35), shifted = plnorm(pmax(0, x - 10), log(15), .5))
}
make_data <- function(profile, replicate) {
  k <- match(profile, profiles$id)
  set.seed(20261122L + 1000L * k + replicate)
  # One result per person; a single age/sex partition, assay and stable period.
  # Both eligible services share physiology; disease fractions differ by service.
  setting <- rep(c("service_a", "service_b"), each = 2500L)
  affected <- runif(5000L) < rep(c(profiles$p_a[k], profiles$p_b[k]), each = 2500L)
  value <- reference_quantile(runif(5000L), profile)
  path <- switch(profile, normal = rnorm(5000L, 42, 4),
                 overlap = rnorm(5000L, 36, 4.5),
                 lognormal = rlnorm(5000L, log(60), .3),
                 shifted = 10 + rlnorm(5000L, log(65), .25))
  value[affected] <- path[affected]
  value <- round(value, 2)
  # Independent of the routine patients. Eligibility never depends on the result.
  anchor <- round(reference_quantile(runif(400L), profile), 2)
  list(value = value, setting = setting, anchor = anchor,
       truth = reference_quantile(c(.025, .975), profile),
       actual_affected = mean(affected))
}

valid_limits <- function(x) length(x) == 2L && all(is.finite(x)) && x[1] < x[2]
fit_anchor <- function(d, transform) {
  fit_anchored_ri(d$value, d$setting, head(d$anchor, 40L), transform = transform)
}

one_case <- function(i) {
  g <- grid[i, ]
  profile <- profiles[g$profile, ]
  cache <- sprintf("tools/guided-cache/%s-%02d-B%d.rds", profile$id, g$replicate, B)
  if (file.exists(cache)) return(readRDS(cache))
  message(sprintf("Starting %s dataset %d", profile$id, g$replicate))
  d <- make_data(profile$id, g$replicate)
  set.seed(810000L + 1000L * g$profile + g$replicate)
  truth <- d$truth
  output <- list()
  row <- function(method, limits, ci, n_screened, successful, failed = 0L,
                  at_bound = FALSE, np_fraction = NA_real_, model = "", error = "") {
    ok <- valid_limits(limits)
    if (!ok) limits <- c(NA_real_, NA_real_)
    data.frame(profile = profile$id, replicate = g$replicate, method = method,
      n_routine = if (method %in% c("direct120", "direct400")) 0L else 5000L,
      n_screened = n_screened, model = model, lower = limits[1], upper = limits[2],
      lower_ci_low = ci[1, 1], lower_ci_high = ci[2, 1],
      upper_ci_low = ci[1, 2], upper_ci_high = ci[2, 2],
      true_lower = truth[1], true_upper = truth[2],
      lower_error = limits[1] - truth[1], upper_error = limits[2] - truth[2],
      mpe = 100 * mean(abs(limits - truth) / truth),
      worst_error = max(abs(limits - truth)),
      reference_coverage = if (ok) diff(reference_cdf(limits, profile$id)) else NA_real_,
      lower_ci_covers = if (ok) ci[1, 1] <= truth[1] && ci[2, 1] >= truth[1] else NA,
      upper_ci_covers = if (ok) ci[1, 2] <= truth[2] && ci[2, 2] >= truth[2] else NA,
      bootstrap_requested = B, bootstrap_successful = successful,
      bootstrap_failed = failed, at_bound = at_bound, np_fraction = np_fraction,
      expected_affected = mean(c(profile$p_a, profile$p_b)),
      actual_affected = d$actual_affected, status = if (ok) "ok" else "failed",
      error = error, stringsAsFactors = FALSE)
  }
  empty_ci <- matrix(NA_real_, 2, 2)
  # Stock refineR, including its own resampling and median-model selection.
  message(sprintf("  refineR %s: %d bootstrap fits on %d workers", profile$refine_model, B, workers))
  ref_start <- proc.time()[["elapsed"]]
  ref <- tryCatch(refineR::findRI(d$value, model = profile$refine_model,
    NBootstrap = B, seed = 410000L + 1000L * g$profile + g$replicate, NCores = workers),
    error = identity)
  message(sprintf("  refineR finished in %.1f minutes", (proc.time()[["elapsed"]] - ref_start)/60))
  if (inherits(ref, "error")) {
    output[[1]] <- row("refineR", c(NA, NA), empty_ci, 0L, 0L, B,
      model = profile$refine_model, error = conditionMessage(ref))
  } else {
    ri <- refineR::getRI(ref, pointEst = "medianBS", CIprop = .95,
                         UCMargins = FALSE)
    n_success <- length(ref$MuBS)
    output[[1]] <- row("refineR", ri$PointEst, rbind(ri$CILow, ri$CIHigh),
      0L, n_success, B - n_success, np_fraction = ref$PMed, model = profile$refine_model)
    # Store original package plots for every dataset, including the CIs.
    grDevices::svg(sprintf("assets/guided/%s-%02d.svg", profile$id, g$replicate),
                   width = 10, height = 5.5)
    plot(ref, pointEst = "medianBS", showPathol = FALSE, showMargin = TRUE, uncertaintyRegion = "bootstrapCI",
         title = paste(profile$title, "| dataset", g$replicate), xlab = "Synthetic units")
    grDevices::dev.off()
    # Show the default model too for the prespecified shifted distribution.
    if (profile$id == "shifted") {
      initial <- refineR::findRI(d$value, model = "BoxCox", NBootstrap = 0L)
      grDevices::svg(sprintf("assets/guided/shifted-default-%02d.svg", g$replicate),
                     width = 10, height = 5.5)
      plot(initial, showPathol = FALSE, title = "Initial one-parameter Box-Cox fit",
           xlab = "Synthetic units")
      grDevices::dev.off()
    }
  }
  # Paired routine cohort. The additional information is 40 independent labels.
  message(sprintf("  anchoredRI: point estimate and %d bootstrap fits", B))
  af <- tryCatch(fit_anchor(d, profile$anchor_transform), error = identity)
  set.seed(610000L + 1000L * g$profile + g$replicate)
  groups <- split(seq_along(d$value), d$setting)
  anchor_draws <- parallel::mclapply(seq_len(B), function(b) {
    set.seed(610000L + 100000L * g$profile + 1000L * g$replicate + b)
    indices <- unlist(lapply(groups, function(x) sample(x, length(x), replace = TRUE)))
    ad <- list(value = d$value[indices], setting = d$setting[indices],
               anchor = sample(head(d$anchor, 40L), 40L, replace = TRUE))
    fit <- tryCatch(fit_anchor(ad, profile$anchor_transform), error = identity)
    if (inherits(fit, "error") || !valid_limits(fit$limits)) c(NA_real_, NA_real_)
    else unname(fit$limits)
  }, mc.cores = workers, mc.set.seed = FALSE)
  anchor_bs <- do.call(cbind, anchor_draws)
  ok <- apply(anchor_bs, 2, function(x) all(is.finite(x)))
  ci <- if (sum(ok) >= 2L) apply(anchor_bs[, ok, drop = FALSE], 1,
               quantile, c(.025, .975), names = FALSE) else empty_ci
  output[[2]] <- row("anchored40", if (inherits(af, "error")) c(NA, NA) else af$limits,
    ci, 40L, sum(ok), sum(!ok),
    at_bound = if (inherits(af, "error")) NA else af$at_bound,
    model = paste0("Gaussian mixture; ", profile$anchor_transform),
    error = if (inherits(af, "error")) conditionMessage(af) else "")
  # Direct benchmarks match the original paper's percentile method, not mean +/- SD.
  message("  Direct percentile estimates and bootstrap CIs")
  for (n in c(120L, 400L)) {
    a <- head(d$anchor, n)
    limits <- quantile(a, c(.025, .975), type = 7, names = FALSE)
    direct_bs <- replicate(B, quantile(sample(a, n, replace = TRUE),
                                       c(.025, .975), type = 7, names = FALSE))
    ci <- apply(direct_bs, 1, quantile, c(.025, .975), names = FALSE)
    output[[length(output) + 1L]] <- row(paste0("direct", n), limits, ci,
      n, B, model = "Empirical quantiles (R type 7)")
  }
  result <- do.call(rbind, output)
  saveRDS(result, cache)
  message(sprintf("Completed %s dataset %d", profile$id, g$replicate))
  result
}
grid <- expand.grid(profile = seq_len(nrow(profiles)), replicate = seq_len(reps))
results <- lapply(seq_len(nrow(grid)), one_case)
stopifnot(all(vapply(results, is.data.frame, logical(1))))
results <- do.call(rbind, results)
stopifnot(nrow(results) == nrow(profiles) * reps * 4L,
  !anyDuplicated(results[c("profile", "replicate", "method")]))
write.csv(results, "assets/guided/fits.csv", row.names = FALSE, na = "")
summary <- do.call(rbind, lapply(split(results, list(results$profile, results$method), drop = TRUE), function(x) {
  good <- x$status == "ok"
  y <- x[good, ]
  data.frame(profile = x$profile[1], method = x$method[1], attempted = nrow(x),
    successful = sum(good), failed = sum(!good),
    mean_lower = mean(y$lower), mean_upper = mean(y$upper),
    lower_mae = mean(abs(y$lower_error)), upper_mae = mean(abs(y$upper_error)),
    lower_bias = mean(y$lower_error), upper_bias = mean(y$upper_error),
    mean_percentage_error = mean(y$mpe), mean_larger_limit_error = mean(y$worst_error),
    mcse_mpe = sd(y$mpe)/sqrt(nrow(y)),
    mean_reference_coverage = mean(y$reference_coverage),
    lower_ci_coverage = mean(y$lower_ci_covers), upper_ci_coverage = mean(y$upper_ci_covers),
    bootstrap_failed = sum(x$bootstrap_failed), boundary_fits = sum(y$at_bound, na.rm = TRUE))
}))
write.csv(summary, "assets/guided/summary.csv", row.names = FALSE)
write.csv(profiles, "assets/guided/profiles.csv", row.names = FALSE)
writeLines(c(R.version.string, paste("refineR", packageVersion("refineR")),
  paste("future", packageVersion("future")), paste("future.apply", packageVersion("future.apply")),
  paste("ash", packageVersion("ash")),
  paste("Replicates per profile:", reps), paste("Bootstrap resamples:", B),
  "5,000 routine results from two eligible services; same routine inputs for both models.",
  "40/120/400 nested independent screened reference individuals, never result-filtered.",
  "Original refineR fit plots, medianBS, 95% bootstrap CIs; no clinical data."),
  "assets/guided/versions.txt")
print(summary, row.names = FALSE)
