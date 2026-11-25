# Run in the Task G container with the Stage 3 dependency snapshot.
library(anchoredRI)
stopifnot(packageVersion("anchoredRI") == "0.4.0", packageVersion("refineR") == "2.0.0",
          getRversion() == "4.5.0")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
fixtures <- read.csv(file.path(args[1L], "fixtures.csv"), na.strings = "NA",
                     colClasses = c(method = "character", first_reason = "character"))
sha <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
compare <- function(observed, expected, label) {
  if (!isTRUE(all.equal(observed, expected, check.attributes = FALSE, tolerance = 0)))
    stop(label, ": expected ", toString(expected), ", found ", toString(observed))
}
rows <- list()
for (shard_name in unique(fixtures$rds)) {
  shard <- readRDS(file.path(args[1L], shard_name))
  for (i in which(fixtures$rds == shard_name)) {
    f <- fixtures[i, ]
    case <- shard$cases[[f$position]]
    test_case <- function() {
      compare(case$j, f$j, "Dataset ID")
      compare(sha(case$routine), f$routine_hash, "Routine vector")
      compare(sha(case$panels$guidance), f$anchor_hash, "Guidance panel")
      compare(length(case$routine), f$N, "Routine sample size")
      fit <- fit_anchored_ri(case$routine, rep("routine", f$N), case$panels$guidance,
                             rounding_unit = f$h, fit_seed = 1L,
                             fallback_seed = f$fallback_seed, pit_seed = f$pit_seed)
      compare(fit$anchoredRI_status, f$first_status, "First fit decision")
      compare(paste(fit$anchoredRI_reason, collapse = " | "), f$first_reason, "First fit reason")
      compare(fit$transform, f$first_transform, "First fit transformation")
      compare(fit$components, f$first_components, "First fit component count")
      compare(if (is.na(fit$method)) "" else fit$method, f$method, "Release route")
      compare(fit$status == "reported", nzchar(f$method), "Final release status")
      fallback <- fit$diagnostics$anchorRefine
      if (f$first_status == "reported") {
        stopifnot(is.null(fallback))
      } else {
        compare(fallback$check$status, f$check_status, "Stage 1 decision")
        compare(fallback$check$reason, f$check_reason, "Stage 1 reason")
      }
      expected <- as.numeric(c(f$lower, f$upper))
      compare(is.na(fit$limits), is.na(expected), "Endpoint missingness")
      error <- if (nzchar(f$method)) max(abs(fit$limits - expected)) else 0
      if (!is.finite(error) || error > 1e-10)
        stop("Expected limits ", toString(expected), ", found ", toString(fit$limits),
             ", max absolute difference ", format(error, digits = 17))
      if (is.na(fit$method)) {
        stopifnot(grepl("120 screened", fit$recommendation, fixed = TRUE),
                  grepl("20 screened", fit$recommendation, fixed = TRUE))
      }
      data.frame(j = f$j, method = f$method, first_status = f$first_status,
                 lower = unname(fit$limits[1L]), upper = unname(fit$limits[2L]),
                 max_endpoint_error = error, error = "")
    }
    rows[[length(rows) + 1L]] <- tryCatch(test_case(), error = function(e) {
      message("Dataset ", f$j, ": ", conditionMessage(e))
      data.frame(j = f$j, method = f$method, first_status = f$first_status,
                 lower = NA_real_, upper = NA_real_, max_endpoint_error = NA_real_,
                 error = conditionMessage(e))
    })
  }
}
results <- do.call(rbind, rows)
write.csv(results, args[2L], row.names = FALSE)
print(table(results$method, nzchar(results$error), useNA = "always"))
stopifnot(!any(nzchar(results$error)))
