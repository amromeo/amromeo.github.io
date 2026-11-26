# SPDX-License-Identifier: GPL-3.0-only
# Independent readback, using only archived inputs and saved model objects.
# Rscript --vanilla verify-records.R INPUT_DIR RECORDS_DIR VERIFICATION_DIR j [j ...]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 4L)
input <- normalizePath(args[[1L]])
records <- normalizePath(args[[2L]])
out <- args[[3L]]
ids <- as.integer(args[-(1:3)])
stopifnot(!anyNA(ids), !anyDuplicated(ids))

# Source the original Stage 2 design and data generator without launching any fitter.
for (file in c("generator.R", "confirmation-design.R", "confirmation-generator.R"))
  source(file.path(input, file))
planned <- stage2_confirmation_design()
design <- read.csv(file.path(input, "confirmation-design.csv"), stringsAsFactors = FALSE,
                   colClasses = vapply(planned, typeof, character(1)))
stopifnot(identical(design$j, planned$j), nrow(design) == 2880L)
design <- design[design$N == 50000L, ]
stopifnot(nrow(design) == 1440L, !anyDuplicated(design$j), all(ids %in% design$j))
references <- readRDS(file.path(input, "references.rds"))
archive <- read.csv(file.path(input, "reference-hashes.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(archive) == 1440L, !anyDuplicated(archive$j), setequal(archive$j, design$j))

vector_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)
file_hash <- function(x) digest::digest(file = x, algo = "sha256")
read_tab <- function(d, name) read.csv(file.path(d, name), stringsAsFactors = FALSE,
                                        na.strings = c("", "NA"), check.names = FALSE)
read_vector <- function(x) if (!length(x)) NA else unname(x)
same <- function(x, y, tolerance = 5e-13) {
  x <- read_vector(x)
  y <- read_vector(y)
  if (length(x) != length(y) || any(is.na(x) != is.na(y))) return(FALSE)
  good <- !is.na(x)
  if (!any(good)) return(TRUE)
  x <- x[good]
  y <- y[good]
  if (is.numeric(x) && is.numeric(y)) {
    finite <- is.finite(x) & is.finite(y)
    return(all(x[!finite] == y[!finite]) &&
           all(abs(x[finite] - y[finite]) <= tolerance * pmax(1, abs(y[finite]))))
  }
  all(as.character(x) == as.character(y))
}
context <- ""
require_check <- function(ok, message) {
  if (!isTRUE(ok)) stop(context, " | ", message, call. = FALSE)
}
compare_fields <- function(row, expected) {
  for (field in names(expected)) {
    require_check(field %in% names(row), paste("missing field", field))
    require_check(same(row[[field]], expected[[field]]), paste("field differs:", field))
  }
}
reason <- function(x) {
  message <- paste(x, collapse = " | ")
  if (nzchar(message)) message else NA_character_
}
valid <- function(limits) length(limits) == 2L && all(is.finite(limits)) && limits[1L] < limits[2L]
files <- c("fits.csv", "models.rds", "positions.rds", "receipt.csv", "routes.csv")
draw_names <- c(sprintf("A%02d", 1:50), sprintf("B%02d", 1:50))
model_names <- c("full", draw_names)

for (j in ids) {
  context <- paste("j", j)
  d <- file.path(records, j)
  require_check(setequal(list.files(d, all.files = TRUE, no.. = TRUE), c(files, "SHA256SUMS")),
                "record file inventory")
  hashes <- vapply(files, function(f) file_hash(file.path(d, f)), character(1))
  received <- readLines(file.path(d, "SHA256SUMS"))
  expected <- paste(hashes, files, sep = "  ")
  require_check(length(received) == length(expected) && setequal(received, expected), "SHA256SUMS")

  meta <- design[design$j == j, , drop = FALSE]
  case <- stage2_confirmation_generate(meta)
  source_row <- archive[archive$j == j, , drop = FALSE]
  ref <- references[[as.character(j)]]
  rh <- vector_hash(case$x)
  ah <- vector_hash(case$anchors)
  require_check(length(case$x) == 50000L && length(case$anchors) == 40L &&
    identical(case$anchors, ref$anchor) && same(rh, source_row$routine_hash) &&
    same(ah, source_row$anchor_hash) && same(ref$routine_hash, rh) && same(ref$anchor_hash, ah),
    "generated input identities vs archived guidance and routine hashes")
  receipt <- read_tab(d, "receipt.csv")
  require_check(nrow(receipt) == 1L, "case receipt")
  compare_fields(receipt, list(j=j, routine_n=50000L, anchor_n=40L, rounding_unit=case$h,
                               routine_hash=rh, anchor_hash=ah))
  models <- readRDS(file.path(d, "models.rds"))
  positions <- readRDS(file.path(d, "positions.rds"))
  fits <- read_tab(d, "fits.csv")
  routes <- read_tab(d, "routes.csv")
  keys <- ifelse(fits$fit_type == "full", "full", sprintf("%s%02d", fits$bank, fits$draw))
  require_check(length(models) == 101L && setequal(names(models), model_names) &&
                length(positions) == 100L && setequal(names(positions), draw_names) &&
                nrow(fits) == 101L && setequal(keys, model_names) && all(fits$j == j) &&
                nrow(routes) == 2L && setequal(routes$bank, c("A","B")) && all(routes$j == j),
                "model, draw and route keys")
  max_p_diff <- 0
  max_limit_diff <- 0
  max_weight_error <- 0
  for (key in model_names) {
    context <- paste("j", j, key)
    row <- fits[match(key, keys), ]
    fit <- models[[key]]
    raw <- fit$diagnostics$anchored$limits
    final <- fit$limits
    safeguard <- fit$safeguard
    fallback <- fit$diagnostics$anchorRefine
    full <- key == "full"
    if (!full) {
      bank <- substring(key, 1L, 1L)
      b <- as.integer(substring(key, 2L))
      seed <- if (bank == "A") 100000000L else 200000000L
      seed <- seed + 100L*j + b
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      ix <- sort.int(sample.int(50000L, 5000L, replace = FALSE))
      require_check(identical(ix, positions[[key]]), "5,000 prespecified sample positions")
      compare_fields(row, list(bank=bank, draw=b, sample_seed=seed, positions_hash=vector_hash(ix),
                               fit_type="subsample"))
      observed <- case$x[ix]
    } else {
      compare_fields(row, list(bank=NA, draw=0, sample_seed=NA, positions_hash=NA, fit_type="full"))
      observed <- case$x
    }
    compare_fields(row, list(first_status=fit$anchoredRI_status, status=fit$status, method=fit$method,
      first_reason=reason(fit$anchoredRI_reason), reason=reason(fit$reason),
      raw_lower=raw["lower"], raw_upper=raw["upper"], lower=final["lower"], upper=final["upper"],
      components=fit$components, transform=fit$transform, shift=fit$shift,
      mu=fit$reference$mu, sigma=fit$reference$sigma, bic=fit$selected_bic,
      optimizer_status=fit$optimizer_status, optimizer_code=fit$optimizer_code,
      p=safeguard$lrt_p, lrt_statistic=safeguard$lrt_statistic, lrt_flag=safeguard$lrt_flag,
      lower_flag=safeguard$lower_flag, lower_check_skipped=FALSE, routine_min=min(observed),
      fit_error=reason(fit$diagnostics$anchored$error), check_error=reason(safeguard$error),
      fallback_status=fallback$status, fallback_reason=reason(fallback$reason),
      fallback_p=fallback$check$p, bic_n=length(observed) + 40L))
    require_check(fit$fit_seed == 1L && same(fit$rounding_unit, case$h) &&
                  same(fit$probs, c(.025,.975)) && isTRUE(safeguard$enabled), "public settings")
    require_check(identical(full && fit$anchoredRI_status != "reported", !is.null(fallback)),
                  "fallback call rule")
    if (!is.null(fallback)) {
      require_check(fallback$fit_seed == 123L && fallback$pit_seed == 38000000L + j &&
                    same(fallback$weight, .25) && same(fallback$rounding_unit, case$h),
                    "fallback seeds and weight")
      checked <- fallback$check
      if (!is.null(checked) && is.finite(checked$p)) {
        require_check(identical(checked$status, if (checked$p < .01) "decline" else "keep"),
                      "fallback check decision")
      }
    }
    chosen <- if (identical(fallback$status, "reported")) "anchorRefine" else
      if (identical(fit$anchoredRI_status, "reported")) "anchoredRI" else NA_character_
    require_check(same(fit$method, chosen) &&
      same(fit$status, if (!is.na(chosen)) "reported" else fit$anchoredRI_status),
      "package reporting method and status")
    if (!is.na(chosen)) {
      require_check(valid(final) && same(final, if (chosen == "anchorRefine") fallback$limits else raw),
                    "released package endpoints")
    } else require_check(length(final) == 2L && all(is.na(final)), "package withholding")

    if (length(fit$weights)) {
      s <- fit$selection
      candidates <- s[is.finite(s$bic_original), ]
      candidate_order <- order(candidates$bic_original,
                               match(candidates$transform,c("identity","log","shifted_log")),
                               candidates$components)
      winner <- candidates[candidate_order[1], ]
      require_check(nrow(candidates) > 0L && all(s$components %in% 1:5) &&
                    same(winner$bic_original, fit$selected_bic) &&
                    same(winner$transform, fit$transform) && same(winner$components, fit$components),
                    "selected original-scale BIC and tie order")
      weights <- as.vector(fit$weights)
      error <- abs(sum(weights) - 1)
      require_check(length(weights) == fit$components && all(is.finite(weights) & weights >= 0) &&
                    error <= 1e-10 && same(row$weight_sum, sum(weights)) && same(row$share, weights[1]),
                    "selected routine mixture weights")
      max_weight_error <- max(max_weight_error, error)
      kpars <- 3L * fit$components - 1L + (fit$transform == "shifted_log")
      require_check(same(row$loglik, (kpars * log(length(observed) + 40) - fit$selected_bic)/2),
                    "selected candidate log likelihood")
      limits <- fit$reference$mu + stats::qnorm(c(.025, .975)) * fit$reference$sigma
      if (!identical(fit$transform, "identity")) limits <- fit$shift + exp(limits)
      require_check(same(limits, raw, 1e-10), "raw limits from component 1")
      if (all(is.finite(limits))) max_limit_diff <- max(max_limit_diff, abs(limits - raw))
    } else require_check(is.na(row$share) && row$first_status == "withheld", "technical fit status")
    if (fit$anchoredRI_status != "withheld") {
      y <- switch(fit$transform, identity=case$anchors, log=log(case$anchors),
                  shifted_log=log(case$anchors - fit$shift))
      sd0 <- sqrt(mean((y - mean(y))^2))
      statistic <- max(0, 2 * sum(stats::dnorm(y, mean(y), sd0, log = TRUE) -
                                    stats::dnorm(y, fit$reference$mu, fit$reference$sigma, log = TRUE)))
      p <- stats::pchisq(statistic, 2, lower.tail=FALSE)
      low <- raw[1L] < min(observed)
      require_check(abs(p - row$p) <= 1e-9 && same(statistic, row$lrt_statistic, 1e-10) &&
        same(p < .01, row$lrt_flag) && same(low, row$lower_flag) &&
        same(fit$anchoredRI_status, if (p < .01 || low) "declined" else "reported"),
        "independent first-stage anchor likelihood ratio, minimum and decision")
      max_p_diff <- max(max_p_diff, abs(p - row$p))
    } else require_check(isTRUE(safeguard$technical) || length(fit$diagnostics$anchored$error) > 0L,
                         "technical first-stage withholding")
  }

  full <- models$full
  for (bank in c("A","B")) {
    context <- paste("j", j, bank, "routes")
    r <- routes[routes$bank == bank, ]
    bank_models <- models[sprintf("%s%02d", bank, 1:50)]
    statuses <- vapply(bank_models, function(m) m$anchoredRI_status, character(1))
    require_check(all(statuses %in% c("declined","reported","withheld")), "50 first-stage statuses")
    pairs <- lapply(bank_models, function(m) m$diagnostics$anchored$limits)
    pairs <- Filter(valid, pairs)
    n <- length(pairs)
    med <- if (n) apply(do.call(rbind, pairs), 2L, stats::median) else c(NA_real_, NA_real_)
    pass <- sum(statuses == "reported")
    decline <- sum(statuses == "declined")
    technical <- 50L - pass - decline
    majority <- if (decline >= 26L) "declined" else if (decline + technical <= 25L) "reported" else "withheld"
    step <- if (full$anchoredRI_status == "reported") "a" else if (full$status == "reported") "b" else
              if (valid(med) && majority == "reported") "c" else "d"
    limits <- switch(step, a=med, b=full$limits, c=med, d=c(NA_real_, NA_real_))
    status <- if (valid(limits)) "reported" else if (step == "d" && full$status == "declined" &&
                  majority == "declined" && valid(med)) "declined" else "withheld"
    compare_fields(r, list(step=step, status=status, lower=limits[1], upper=limits[2],
      median_lower=med[1], median_upper=med[2], contributors=n, votes_pass=pass, votes_decline=decline,
      votes_technical=technical, majority=majority, first_status=full$anchoredRI_status,
      package_status=full$status, package_method=full$method, package_lower=full$limits[1],
      package_upper=full$limits[2], first_reason=reason(full$anchoredRI_reason),
      package_reason=reason(full$reason), fallback_status=full$diagnostics$anchorRefine$status,
      fallback_reason=reason(full$diagnostics$anchorRefine$reason)))
  }

  dest <- file.path(out, j)
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  all_files <- c(files, "SHA256SUMS")
  all_paths <- file.path(d, all_files)
  hash_table <- data.frame(j=j, file=all_files, bytes=file.size(all_paths),
                           sha256=unname(vapply(all_paths, file_hash, character(1))))
  checks <- data.frame(j=j, routine_hash=rh, anchor_hash=ah, n_models=101L, n_draws=100L,
    n_positions=sum(lengths(positions)), n_routes=2L, technical_fits=sum(fits$first_status == "withheld"),
    max_abs_lrt_p_difference=max_p_diff, max_abs_limit_difference=max_limit_diff,
    max_abs_component_sum_error=max_weight_error, verified=TRUE)
  write.csv(hash_table, file.path(dest, "hashes.csv"), row.names=FALSE)
  write.csv(checks, file.path(dest, "records.csv"), row.names=FALSE)
  message("Verified all 101 fits, 500,000 stored positions and two routes for j=", j)
}
