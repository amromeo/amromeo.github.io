# Rscript replay-task-e.R TASK_E_PACKAGE SMALL_RIBENCH_DIR LARGE_INPUT_SHARDS
# Replay archived Arm B fits through the installed package, without re-optimizing.
library(anchoredRI)
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
study <- normalizePath(args[1])
small_data <- normalizePath(args[2])
shards <- normalizePath(args[3])
study_tools <- file.path(study, "tools/boxcox")
source(file.path(study_tools, "inputs.R"))
source(file.path(study, "tools/stress-design.R"))
small <- file.path(study_tools, "accuracy/results")
large <- file.path(study_tools, "accuracy/large/results")
read_parts <- function(dir) {
  files <- sort(list.files(dir, "^part-[0-9]+[.]csv$", full.names = TRUE))
  stopifnot(length(files) > 0L)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

# Component 1 is shared by both settings. Keep one identical parameter record per fit.
params <- function(x) unique(x[x$method == "armB" & x$screened_component,
  c("dataset_id", "transform", "shift", "transformed_mean", "transformed_sd")])
ref <- rbind(params(read_parts(file.path(small, "selected-component-parameters"))),
             params(read.csv(file.path(large, "selected-component-parameters.csv"))))
ds <- read_parts(file.path(small, "diagnostics"))
ds <- ds[ds$method == "armB", ]
lg <- read.csv(file.path(large, "armB-fits.csv"))
stopifnot(all(lg$method == "armB"), all(lg$max_components == 5))
lg$stage <- "50000"
lg$reported <- lg$report_status == "reported"
lg$code <- lg$convergence_code
cols <- c("dataset_id", "stage", "transform", "components", "shift", "lower", "upper",
          "lrt_p", "lower_flag", "declined", "technical", "reported", "code")
expected <- rbind(ds[cols], lg[cols])
meta <- read.csv(file.path(study_tools, "anchorRefine-comparator/anchoredRI-comparator.csv"))
stopifnot(nrow(expected) == 5220L, !anyDuplicated(expected$dataset_id),
          !anyDuplicated(ref$dataset_id), !anyDuplicated(meta$dataset_id),
          setequal(expected$dataset_id, ref$dataset_id), setequal(expected$dataset_id, meta$dataset_id))
ids <- expected$dataset_id
ref <- ref[match(ids, ref$dataset_id), ]
meta <- meta[match(ids, meta$dataset_id), ]
stopifnot(all(expected$transform == ref$transform), all(expected$shift == ref$shift),
          all(expected$stage == meta$stage))
index <- setNames(seq_along(ids), ids)
results <- vector("list", length(ids))
apply_safeguard <- getFromNamespace(".apply_ri_safeguard", "anchoredRI")
infer_unit <- getFromNamespace(".measurement_unit", "anchoredRI")
rank_models <- getFromNamespace(".rank_ri_candidates", "anchoredRI")

candidate_cols <- c("dataset_id", "transform", "k", "bic", "code")
choices <- rbind(read_parts(file.path(small, "candidates"))[candidate_cols],
                 read.csv(file.path(large, "candidates.csv"))[candidate_cols])
names(choices) <- c("dataset_id", "transform", "components", "bic_original", "code")
by_id <- split(choices, choices$dataset_id)
gap_cols <- c("dataset_id", "method", "second", "second_k", "bic_gap")
gaps <- rbind(read_parts(file.path(small, "runner-up-gaps"))[gap_cols],
              read.csv(file.path(large, "runner-up-gaps.csv"))[gap_cols])
gaps <- gaps[gaps$method == "armB", ]
stopifnot(setequal(ids, names(by_id)), !anyDuplicated(gaps$dataset_id), setequal(ids, gaps$dataset_id))
gaps <- gaps[match(ids, gaps$dataset_id), ]

replay <- function(id, value, anchor, probs) {
  i <- unname(index[[id]])
  row <- expected[i, ]
  p <- ref[i, ]
  stopifnot(is.null(results[[i]]),
            identical(digest::digest(value, algo = "sha256"), meta$routine_hash[i]),
            identical(digest::digest(anchor, algo = "sha256"), meta$anchor_hash[i]),
            min(value) == meta$routine_min[i])
  y <- p$transformed_mean + qnorm(c(.025, .975)) * p$transformed_sd
  central <- if (p$transform == "identity") y else p$shift + exp(y)
  unit <- infer_unit(unique(c(value, anchor)), NULL)
  candidates <- by_id[[id]]
  ranking <- rank_models(candidates)
  winner <- candidates[ranking$winner, ]
  second <- ranking$runner_up
  # The archive stores decimal text for candidate BICs. Allow roundoff in their difference.
  gap_tolerance <- 128 * .Machine$double.eps * max(1, abs(winner$bic_original), abs(second$bic_original))
  fit <- list(limits = setNames(central, c("lower", "upper")), transform = p$transform,
              shift = p$shift, components = row$components, rounding_unit = unit, runner_up = second,
              optimizer_code = winner$code,
              reference = list(mu = p$transformed_mean, sigma = p$transformed_sd))
  result <- apply_safeguard(fit, value, anchor, probs = probs)
  check <- result$safeguard
  lrt <- row$lrt_p < .01
  status <- if (row$technical) "withheld" else if (row$declined) "declined" else "reported"
  reasons <- c(if (lrt) "routine data disagree with the screened sample",
               if (row$lower_flag) "estimated lower limit is below every observed result")
  if (is.null(reasons)) reasons <- character()
  raw <- if (length(probs) == 1L) row$upper else c(row$lower, row$upper)
  close <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-10))
  matches <- c(
    rounding_unit = result$rounding_unit == meta$rounding_step[i],
    winner = identical(winner$transform, row$transform) && winner$components == row$components,
    optimizer_code = identical(result$optimizer_code, row$code),
    runner_up = identical(second$transform, gaps$second[i]) && second$components == gaps$second_k[i],
    bic_gap = abs(second$bic_gap - gaps$bic_gap[i]) <= gap_tolerance,
    lrt = identical(check$lrt_flag, lrt),
    lower = identical(check$lower_flag, row$lower_flag),
    skipped = identical(check$lower_check_skipped, length(probs) == 1L),
    technical = identical(result$technical_withholding, row$technical),
    declined = identical(result$declined, row$declined),
    status = identical(result$status, status),
    reported = identical(result$status == "reported", row$reported),
    reason = identical(as.character(result$reason), as.character(reasons)),
    p = close(check$lrt_p, row$lrt_p),
    raw = close(result$diagnostics$anchored$limits, raw),
    released = if (row$reported) close(result$limits, raw) else all(is.na(result$limits))
  )
  if (!all(matches)) stop(id, " differs on ", paste(names(matches)[!matches], collapse = ", "))
  results[[i]] <<- data.frame(stage = row$stage, lrt = lrt, lower = row$lower_flag,
                              both = lrt && row$lower_flag, declined = result$declined,
                              reported = row$reported, one_sided = length(probs) == 1L)
}
defs <- as.data.frame(RIbench::loadTestsetDefinition())
for (i in which(meta$stage == "small")) {
  s <- defs[match(meta$Index[i], defs$Index), ]
  d <- bc_ribench_input(s, small_data, meta$anchor_hash[i])
  replay(ids[i], d$value, d$anchor, d$probs)
}
for (i in which(meta$stage == "grid")) {
  d <- stress_data(meta$case[i], meta$replicate[i])
  replay(ids[i], d$value, d$anchor, c(.025, .975))
}
for (file in sort(list.files(shards, "^inputs-[0-9]+[.]rds$", full.names = TRUE))) {
  for (d in readRDS(file)) replay(paste0("r-", d$s$Index), d$value, d$anchor, d$probs)
}
stopifnot(all(vapply(results, is.data.frame, logical(1))))
rows <- do.call(rbind, results)
counts <- aggregate(cbind(fits = 1L + 0L * lrt, lrt, lower, both, declined, reported, one_sided) ~ stage,
                    rows, sum)
counts <- counts[match(c("small", "grid", "50000"), counts$stage), ]
stopifnot(identical(as.integer(counts$fits), c(2880L, 900L, 1440L)),
          identical(as.integer(counts$declined), c(99L, 91L, 305L)))
print(counts, row.names = FALSE)
cat("All 5,220 Arm B decisions matched, including both flags, one-sided skips and reporting status.\n")
cat("All 5,220 rounding units, model rankings, BIC gaps and winning optimizer codes matched Task E.\n")
