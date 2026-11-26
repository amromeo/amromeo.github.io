# SPDX-License-Identifier: GPL-3.0-only
# Run from the study repository root. The input adapter verifies the Stage 2
# archives and gives the estimator only its measurements, panel and precision.
source("anchoredRI/tools/subsample/pilot.R")
confirmation_probs <- c(.025, .975)

confirmation_route <- function(full, m, votes, majority) {
  usable <- h_usable(m$lower, m$upper, confirmation_probs)
  lower <- upper <- NA_real_
  reason <- ""
  if (full$first_status == "reported") {
    step <- "a"
    status <- if (usable) "reported" else "withheld"
    if (usable) {
      lower <- m$lower
      upper <- m$upper
    } else reason <- "unavailable subsample median"
  } else if (full$status == "reported") {
    step <- "b"
    if (h_usable(full$lower, full$upper, confirmation_probs)) {
      status <- "reported"
      lower <- full$lower
      upper <- full$upper
    } else {
      status <- "withheld"
      reason <- "unusable package interval"
    }
  } else if (majority == "reported" && usable) {
    step <- "c"
    status <- "reported"
    lower <- m$lower
    upper <- m$upper
  } else {
    step <- "d"
    status <- if (full$status == "declined" && majority == "declined" && usable) "declined" else "withheld"
    reason <- paste(c(if (!usable) "unavailable subsample median",
                      if (majority == "declined") "majority declined",
                      if (majority == "withheld") "unresolved majority",
                      if (nzchar(full$reason)) full$reason,
                      if (nzchar(full$fallback_reason)) full$fallback_reason), collapse = " | ")
  }
  data.frame(step, status, reason, lower, upper, median_lower = m$lower, median_upper = m$upper,
             contributors = m$contributors, as.data.frame(as.list(votes)), majority,
             first_status = full$first_status, first_reason = full$first_reason,
             package_status = full$status, package_method = full$method, package_reason = full$reason,
             package_lower = full$lower, package_upper = full$upper, fallback_status = full$fallback_status,
             fallback_reason = full$fallback_reason, stringsAsFactors = FALSE)
}

confirmation_fit <- function(value, anchor, unit, j, fallback) {
  fit <- function() suppressMessages(anchoredRI::fit_anchored_ri(
    value, rep("pooled", length(value)), anchor, transform = "auto", safeguard = TRUE,
    probs = confirmation_probs, rounding_unit = unit, fit_seed = 1L, fallback = fallback,
    fallback_seed = 123L, pit_seed = 38000000L + j))
  # The public full call handles first-stage optimization errors and exercises
  # the actual package fallback. An unexpected full-call error stops the job.
  if (fallback) return(fit())
  tryCatch(fit(), error = function(error) h_missing_fit(error, value, confirmation_probs, unit, 1L))
}

confirmation_row <- function(fit, j, n, bank = NA_character_, draw = 0L,
                             sample_seed = NA_integer_, positions_hash = NA_character_) {
  key <- data.frame(j = j, bank = bank, draw = draw, fit_type = if (draw == 0) "full" else "subsample",
                    sample_seed = sample_seed, positions_hash = positions_hash)
  p <- fit$diagnostics$anchorRefine$check$p
  cbind(key, h_fit_row(fit, n + 40L), fallback_p = if (length(p) == 1L) p else NA_real_)
}

confirmation_write_manifest <- function(dir) {
  files <- c("fits.csv", "models.rds", "positions.rds", "receipt.csv", "routes.csv")
  sums <- vapply(file.path(dir, files), digest::digest, character(1), algo = "sha256", file = TRUE,
                 USE.NAMES = FALSE)
  writeLines(paste(sums, files, sep = "  "), file.path(dir, "SHA256SUMS"))
}

confirmation_worker <- function(input_dir, output_dir, j) {
  stopifnot(as.character(getRversion()) == "4.5.0", as.character(packageVersion("anchoredRI")) == "0.4.0")
  target <- file.path(output_dir, j)
  if (dir.exists(target)) stop("Dataset directory already exists: ", target)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create ", output_dir)
  staging <- tempfile(pattern = paste0(".working-", j, "-"), tmpdir = output_dir)
  if (!dir.create(staging)) stop("Cannot create staging directory: ", staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  x <- confirmation_read_case(input_dir, j)
  n <- length(x$value)
  stopifnot(n == 50000L, length(x$anchor) == 40L, length(x$rounding_unit) == 1L,
            x$rounding_unit > 0, is.finite(x$rounding_unit), all(is.finite(x$value)),
            all(is.finite(x$anchor)), identical(h_vector_hash(x$value), x$routine_hash),
            identical(h_vector_hash(x$anchor), x$anchor_hash))

  models <- positions <- rows <- routes <- list()
  models$full <- confirmation_fit(x$value, x$anchor, x$rounding_unit, j, TRUE)
  rows$full <- confirmation_row(models$full, j, n)
  for (bank in c("A", "B")) {
    keys <- sprintf("%s%02d", bank, 1:50)
    for (draw in seq_along(keys)) {
      key <- keys[draw]
      idx <- h_sample_positions(j, n, bank, draw)
      positions[[key]] <- idx
      fit <- confirmation_fit(x$value[idx], x$anchor, x$rounding_unit, j, FALSE)
      models[[key]] <- fit
      rows[[key]] <- confirmation_row(fit, j, length(idx), bank, draw, h_sample_seed(j, bank, draw),
                                      h_vector_hash(idx))
    }
    draws <- do.call(rbind, rows[keys])
    m <- h_median(draws, confirmation_probs)
    route <- confirmation_route(rows$full, m, h_votes(draws$first_status), h_majority(draws$first_status))
    routes[[bank]] <- cbind(j = j, bank = bank, route)
    message("Completed ", j, " bank ", bank)
  }
  rownames_table <- function(rows) {
    x <- do.call(rbind, rows)
    row.names(x) <- NULL
    x
  }
  h_csv(rownames_table(rows), file.path(staging, "fits.csv"))
  h_csv(rownames_table(routes), file.path(staging, "routes.csv"))
  h_csv(data.frame(j, routine_n = n, anchor_n = length(x$anchor), rounding_unit = x$rounding_unit,
                   routine_hash = x$routine_hash, anchor_hash = x$anchor_hash), file.path(staging, "receipt.csv"))
  saveRDS(models, file.path(staging, "models.rds"), version = 2, compress = "gzip")
  saveRDS(positions, file.path(staging, "positions.rds"), version = 2, compress = "gzip")
  confirmation_write_manifest(staging)
  if (!file.rename(staging, target)) stop("Cannot publish completed case ", j)
  message("Completed and preserved 101 fits for j = ", j)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) < 3 || any(!grepl("^[1-9][0-9]*$", args[-(1:2)]))) {
    stop("Usage: worker.R INPUT_DIR OUTPUT_DIR j [j ...]")
  }
  source("anchoredRI/tools/subsample/confirmation/input-case.R")
  future::plan(future::sequential)
  for (j in as.integer(args[-(1:2)])) confirmation_worker(args[1], args[2], j)
}
