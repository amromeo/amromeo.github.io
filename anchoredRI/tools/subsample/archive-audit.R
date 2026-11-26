# SPDX-License-Identifier: GPL-3.0-only
# Compare the unweighted full-data cap-3 first stage with Task E.
# No subsampling, fallback or weighted fits are reachable from this script.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 5L)) stop("Usage: archive-audit.R INPUT_DIR OUT SHARD [PART PARTS]")
input_dir <- normalizePath(args[[1L]])
output_dir <- args[[2L]]
shard <- as.integer(args[[3L]])
part <- if (length(args) == 5L) as.integer(args[[4L]]) else 1L
parts <- if (length(args) == 5L) as.integer(args[[5L]]) else 1L
stopifnot(shard %in% 1:48, part >= 1L, part <= parts, parts <= 30L)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
study <- dirname(normalizePath(script))
source(file.path(study, "inputs.R"))
source(file.path(study, "mixture.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
read_input <- function(name) read.csv(file.path(input_dir, name), stringsAsFactors = FALSE,
                                      na.strings = c("", "NA"), check.names = FALSE)
roster <- read_input("task-g-index.csv")
baseline <- read_input("baseline-fits.csv")
weights <- read_input("selected-component-parameters.csv")
weights <- weights[weights$method == "anchoredRI" & weights$setting == "pooled", ]
weight_sums <- tapply(weights$routine_weight, weights$dataset_id, sum)
healthy <- weights[weights$component == 1L & weights$screened_component, ]
stopifnot(nrow(roster) == 1440L, !anyDuplicated(roster$Index), nrow(baseline) == 1440L,
          !anyDuplicated(baseline$Index), all(baseline$method == "anchoredRI"),
          all(baseline$max_components == 3L), nrow(healthy) == 1440L,
          !anyDuplicated(healthy$dataset_id), setequal(roster$Index, baseline$Index),
          setequal(healthy$dataset_id, baseline$dataset_id), all(abs(weight_sums - 1) <= 1e-10))
ids <- roster$Index[roster$shard == shard]
stopifnot(length(ids) == 30L)
ids <- ids[(seq_along(ids) - 1L) %% parts == part - 1L]
fit_fields <- c("transform", "components", "lrt_flag", "lower_flag", "status", "raw_lower",
                "raw_upper", "lower", "upper", "p", "share", "shift", "loglik", "bic",
                "routine_min", "lrt_statistic")
endpoints <- c("raw_lower", "raw_upper", "lower", "upper")

archive_values <- function(b, a) {
  stopifnot(b$transform == a$transform, b$shift == a$shift)
  data.frame(transform = b$transform, components = b$components, lrt_flag = b$lrt_p < .01,
             lower_flag = b$lower_flag, status = b$report_status, raw_lower = b$lower,
             raw_upper = b$upper, lower = b$reported_lower, upper = b$reported_upper,
             p = b$lrt_p, share = a$routine_weight, shift = b$shift, loglik = b$loglik,
             bic = b$bic, routine_min = b$routine_min, lrt_statistic = b$lrt_statistic)
}

# Unavailable required quantities remain exceptions; jointly missing release values
# are expected when both methods decline. CRP never requests a lower endpoint.
compare <- function(observed, ref, two_sided, width) {
  reason <- character()
  for (field in c("transform", "components", "lrt_flag", if (two_sided) "lower_flag", "status")) {
    if (!isTRUE(observed[[field]] == ref[[field]])) reason <- c(reason, field)
  }
  for (field in endpoints) {
    new <- observed[[field]]
    old <- ref[[field]]
    if (is.na(new) != is.na(old)) {
      reason <- c(reason, paste0(field, "_missingness"))
    } else if (!is.na(old) && (!is.finite(new) || !is.finite(old) || abs(new - old) > .001 * width)) {
      reason <- c(reason, field)
    }
  }
  if (!is.finite(observed$p) || !is.finite(ref$p) || abs(observed$p - ref$p) > .0001)
    reason <- c(reason, "p")
  reason
}

run_fit <- function(case, seed, label, ref, b, width) {
  start <- proc.time()
  warnings <- character()
  result <- withCallingHandlers(
    h_mixture(case$value, case$anchor, case$rounding_unit, probs = case$probs, caps = 3L,
              routine_weight = 1, fit_seed = as.integer(seed), fallback = FALSE,
              pit_seed = 38000000L + case$Index),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  timing <- proc.time() - start
  fit <- result$fits[["3"]]
  stopifnot(!is.null(fit))
  observed <- h_fit_row(fit, bic_n = result$bic_n)
  stopifnot(nrow(observed) == 1L, observed$status == observed$first_status,
            all(fit_fields %in% names(observed)))
  reasons <- compare(observed, ref, length(case$probs) == 2L, width)
  key <- data.frame(Index = b$Index, Analyte = b$Analyte, fraction = b$pathological_fraction,
                    shard = shard, part = part, attempt = label, fit_seed = seed,
                    archive_seed = b$Index, anchor_seed = b$anchor_seed,
                    routine_hash = b$routine_hash, anchor_hash = b$anchor_hash,
                    routine_n = b$n_routine, rounding_step = b$rounding_step,
                    lower_probability = b$lower_probability, upper_probability = b$upper_probability,
                    archive_width = width)
  both <- cbind(setNames(ref[fit_fields], paste0("archive_", fit_fields)),
                setNames(observed[fit_fields], paste0("observed_", fit_fields)))
  differences <- data.frame(share_difference_pp = 100 * (observed$share - ref$share),
                            elapsed_s = unname(timing["elapsed"]),
                            cpu_s = sum(timing[c("user.self", "sys.self", "user.child", "sys.child")]),
                            pass = !length(reasons), differences = paste(reasons, collapse = "|"),
                            fit_reason = observed$first_reason, check_error = observed$check_error,
                            share_reason = observed$share_reason, weight_sum = observed$weight_sum,
                            warnings = paste(unique(warnings), collapse = "|"))
  list(row = cbind(key, both, differences), fit = fit, engine_errors = result$errors)
}

records <- list()
for (id in ids) {
  case <- h_read_case(input_dir, id)
  b <- baseline[match(id, baseline$Index), ]
  a <- healthy[match(b$dataset_id, healthy$dataset_id), ]
  m <- case$meta
  requested <- if (is.na(b$lower_probability)) b$upper_probability else c(b$lower_probability, b$upper_probability)
  stopifnot(id == case$Index, id == m$Index, b$dataset_id == paste0("r-", id), b$N == 50000L,
            b$n_routine == length(case$value), b$n_anchor == length(case$anchor),
            b$n_anchor == 40L, b$rounding_step == case$rounding_unit, b$anchor_seed == 17000000L + id,
            identical(as.numeric(case$probs), requested), b$routine_hash == m$routine_hash,
            b$anchor_hash == m$anchor_hash, b$ribench_hashcode == m$ribench_hashcode,
            (length(requested) == 1L) == (b$Analyte == "CRP"))
  width <- if (length(requested) == 2L) b$upper - b$lower else b$upper - 0
  stopifnot(is.finite(width), width > 0)
  reference <- archive_values(b, a)
  primary <- run_fit(case, 1L, "study", reference, b, width)
  records[[length(records) + 1L]] <- primary
  if (!primary$row$pass)
    records[[length(records) + 1L]] <- run_fit(case, id, "archived_seed", reference, b, width)
  cat(sprintf("Index %d complete; archived seed rerun %s\n", id, !primary$row$pass))
}

stem <- file.path(output_dir, sprintf("cap3-shard-%02d-part-%02d", shard, part))
rows <- do.call(rbind, lapply(records, "[[", "row"))
# Base write.csv ordinarily emits 15 significant digits. Preserve comparisons
# at tolerance boundaries and each component share as a round-trippable double.
rows[] <- lapply(rows, function(value) {
  if (!is.double(value)) return(value)
  ifelse(is.na(value), NA_character_, sprintf("%.17g", value))
})
saveRDS(records, paste0(stem, ".rds.tmp"), version = 3)
write.csv(rows, paste0(stem, ".csv.tmp"), row.names = FALSE, na = "")
stopifnot(file.rename(paste0(stem, ".rds.tmp"), paste0(stem, ".rds")),
          file.rename(paste0(stem, ".csv.tmp"), paste0(stem, ".csv")))
