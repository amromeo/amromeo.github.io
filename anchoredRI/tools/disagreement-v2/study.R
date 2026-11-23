# Run from the anchoredRI package directory.
study_dir <- "tools/disagreement-v2"
results_dir <- file.path(study_dir, "results")
work_dir <- Sys.getenv("B2_WORK", file.path(study_dir, "work"))
thresholds <- c(.50, .55, .60, .65, .70, .75, .80)

frozen_sources <- c("R/fit_anchored_ri.R", "R/fit_anchored_shifted.R",
                    "tools/stress-design.R", "tools/disagreement/run.R",
                    "tools/anchors-only/estimators.R")
frozen_md5 <- c("cc19e27c8929097fe0fce7b0a23e6fa8", "61b1697e8a55e640fdf4c4d65ae22f10",
                "360d92a5d7479c975b222a042ddb5891", "dcc0a34ebc6ba9797fcc96bfde6ad478",
                "44461b63fb3308f25646a849034565d0")
stopifnot(identical(unname(tools::md5sum(frozen_sources)), frozen_md5))
source("tools/stress-design.R")

frozen_fit <- new.env(parent = globalenv())
sys.source("R/fit_anchored_ri.R", envir = frozen_fit)
sys.source("R/fit_anchored_shifted.R", envir = frozen_fit)
component_fit <- new.env(parent = frozen_fit)
sys.source(file.path(study_dir, "component-fits.R"), envir = component_fit)
component_fit$fit_anchored_ri <- frozen_fit$fit_anchored_ri
environment(component_fit$fit_anchored_ri) <- component_fit

# The historical Step B script runs other analyses at top level. Load only its
# LRT and the likelihood helper, directly from the pinned sources.
read_definition <- function(path, name, env) {
  code <- parse(path, keep.source = FALSE)
  matches <- vapply(code, function(expr) is.call(expr) &&
    identical(expr[[1L]], as.name("<-")) && identical(expr[[2L]], as.name(name)), logical(1))
  stopifnot(sum(matches) == 1L)
  eval(code[[which(matches)]], envir = env)
}
lrt_env <- new.env(parent = globalenv())
read_definition("tools/anchors-only/estimators.R", "ao_normal", lrt_env)
read_definition("tools/disagreement/run.R", "lrt", lrt_env)

check_equivalence <- function(frozen, local) {
  check <- function(field, a, b, categorical = FALSE) {
    equal <- if (categorical) identical(unname(a), unname(b)) else
      length(a) == 1L && length(b) == 1L && is.finite(a) && is.finite(b) && abs(a - b) <= 1e-10
    if (!isTRUE(equal)) stop(sprintf("%s differs: frozen %s, local %s",
      field, paste(format(a, digits = 17), collapse = ","), paste(format(b, digits = 17), collapse = ",")))
  }
  check("lower limit", frozen$limits[1L], local$limits[1L])
  check("upper limit", frozen$limits[2L], local$limits[2L])
  check("transform", frozen$transform, local$transform, TRUE)
  check("offset", frozen$shift, local$shift)
  check("component count", frozen$components, local$components, TRUE)
}

component_summary <- function(fit, setting) {
  k <- unname(fit$components)
  stopifnot(length(fit$component_mu) == k, length(fit$component_sd) == k,
            all(is.finite(fit$component_mu)), all(is.finite(fit$component_sd)), all(fit$component_sd > 0))
  w <- fit$weights[, c("service_a", "service_b"), drop = FALSE]
  n <- as.numeric(table(factor(setting, levels = colnames(w))))
  stopifnot(sum(n) == length(setting), all(n > 0), all(is.finite(w)), all(w > 0))
  pooled <- drop(w %*% (n / sum(n)))
  z <- fit$component_mu[1L] + qnorm(c(.025, .975)) * fit$component_sd[1L]
  inside <- pnorm(z[2L], fit$component_mu, fit$component_sd) -
            pnorm(z[1L], fit$component_mu, fit$component_sd)
  disease <- seq_len(k) > 1L
  overlap <- if (k == 1L) NA_real_ else weighted.mean(inside[disease], pooled[disease])
  if (k > 1L) stopifnot(is.finite(overlap), overlap >= 0, overlap <= 1)
  list(O = overlap, components = data.frame(component = seq_len(k), mu_internal = fit$component_mu,
    sd_internal = fit$component_sd, weight_a = w[, 1L], weight_b = w[, 2L],
    pooled_weight = pooled, probability_inside = inside))
}

evaluate_fits <- function(value, setting, anchor) {
  frozen <- frozen_fit$fit_anchored_ri(value, setting, anchor, transform = "auto")
  local <- component_fit$fit_anchored_ri(value, setting, anchor, transform = "auto")
  check_equivalence(frozen, local)
  s <- component_summary(local, setting)
  test <- lrt_env$lrt(anchor, frozen$limits, frozen$transform, frozen$shift)
  row <- data.frame(diagnostic_lower = unname(frozen$limits[1L]),
    diagnostic_upper = unname(frozen$limits[2L]), transform = frozen$transform,
    offset = unname(frozen$shift), components = unname(frozen$components), O = s$O,
    lrt_statistic = unname(test[1L]), lrt_p = unname(test[2L]),
    lower_abs_difference = unname(abs(frozen$limits[1L] - local$limits[1L])),
    upper_abs_difference = unname(abs(frozen$limits[2L] - local$limits[2L])),
    offset_abs_difference = unname(abs(frozen$shift - local$shift)))
  list(row = row, components = s$components)
}

# Vectorized for the paired candidate comparisons. An absent O is expected only for K = 1.
decide <- function(row, rule, threshold = NA_real_) {
  stopifnot(rule %in% c("a", "b", "lrt"),
            all((row$components == 1L) == is.na(row$O)))
  use_O <- rule != "lrt" & row$components > 1L
  use_lrt <- rule != "a" | row$components == 1L
  if (rule != "lrt") stopifnot(length(threshold) == 1L, threshold %in% thresholds)
  separation_flag <- use_O & row$O > threshold
  lrt_available <- is.finite(row$lrt_p)
  positive_lrt <- lrt_available & row$lrt_p < .01
  lrt_flag <- use_lrt & positive_lrt
  technical <- use_lrt & !lrt_available & !separation_flag
  declined <- separation_flag | lrt_flag | technical
  message <- rep("", nrow(row))
  message[technical] <- "LRT could not be evaluated"
  message[lrt_flag] <- "routine data disagree with the screened sample"
  message[separation_flag] <- "routine results are too mixed to separate healthy from affected patients"
  data.frame(separation_flag = separation_flag, lrt_flag = lrt_flag,
    lrt_unavailable = !lrt_available, technical_withholding = technical, declined = declined,
    reported_lower = ifelse(declined, NA_real_, row$diagnostic_lower),
    reported_upper = ifelse(declined, NA_real_, row$diagnostic_upper), message = message,
    recommendation = ifelse(declined,
      "Establish an interval with at least 120 screened patients or verify a suitable published interval with 20 screened patients.", ""))
}

keys_for_phase <- function(phase) {
  replicates <- switch(phase, development = 1:20, confirmation = 21:40, stop("Unknown phase"))
  expand.grid(case = stress_design$case, replicate = replicates, KEEP.OUT.ATTRS = FALSE)
}

checkpoint <- function(dir, case, replicate) file.path(dir, sprintf("%02d-%02d.rds", case, replicate))

save_atomic <- function(object, path) {
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) stop("Could not save ", path)
}

selection_path <- file.path(results_dir, "selected-rule.csv")
locked_choice <- function() {
  if (!file.exists(selection_path)) stop("Commit the development selection before confirmation.")
  sha <- system2("git", c("log", "-1", "--format=%H", "--", selection_path), stdout = TRUE)
  if (length(sha) != 1L || !grepl("^[0-9a-f]{40}$", sha)) stop("Missing selection commit.")
  files <- file.path(study_dir, c("PROTOCOL.md", "component-fits.R", "study.R", "run.R", "report.R", "test.R"))
  needed <- c(files, selection_path)
  if (system2("git", c("ls-files", "--error-unmatch", "--", needed), stdout = FALSE) != 0L ||
      system2("git", c("diff", "--quiet", sha, "--", needed)) != 0L)
    stop("The study code, selection and protocol must match the selection commit.")
  choice <- read.csv(selection_path, stringsAsFactors = FALSE)
  if (nrow(choice) != 1L || choice$status != "selected") stop("Development produced no eligible rule.")
  stopifnot(choice$rule %in% c("a", "b"), choice$threshold %in% thresholds)
  list(rule = choice$rule, threshold = choice$threshold, commit = sha)
}

development_index <- function() {
  path <- "assets/comparison-study/summary/comparison-results.zip"
  prior <- read.csv(unz(path, "all-results.csv"), stringsAsFactors = FALSE)
  prior <- prior[prior$method == "anchored40", c("case", "replicate", "seed")]
  key <- function(x) paste(x$case, x$replicate, sep = "/")
  expected <- key(keys_for_phase("development"))
  stopifnot(nrow(prior) == length(expected), !anyDuplicated(key(prior)), setequal(key(prior), expected))
  setNames(prior$seed, key(prior))
}
