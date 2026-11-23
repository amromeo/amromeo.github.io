source("tools/disagreement-v2/study.R")

run_phase <- function(phase, workers) {
  keys <- keys_for_phase(phase)
  stopifnot(length(workers) == 1L, is.finite(workers), workers >= 1)
  choice <- if (phase == "confirmation") locked_choice() else NULL
  old_seeds <- if (phase == "development") development_index() else NULL
  dir <- file.path(work_dir, phase)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(dir, "manifest.rds")
  manifest <- list(phase = phase, sources = tools::md5sum(c(frozen_sources,
    file.path(study_dir, c("component-fits.R", "study.R", "run.R")))),
    runtime = R.version.string, platform = R.version$platform,
    blas = unname(extSoftVersion()["BLAS"]), choice = choice)
  if (file.exists(manifest_path) && !identical(readRDS(manifest_path), manifest))
    stop("Checkpoints have different inputs. Use a fresh B2_WORK directory.")
  save_atomic(manifest, manifest_path)
  halt_path <- file.path(dir, "HALT.txt")
  if (file.exists(halt_path)) stop(paste(readLines(halt_path), collapse = "\n"))
  paths <- checkpoint(dir, keys$case, keys$replicate)
  pending <- which(!file.exists(paths))

  one <- function(i) {
    key <- keys[i, , drop = FALSE]
    case <- key$case
    replicate <- key$replicate
    tryCatch({
      d <- stress_data(case, replicate)
      if (!is.null(old_seeds)) stopifnot(d$seed == old_seeds[[paste(case, replicate, sep = "/")]])
      fit <- evaluate_fits(d$value, d$setting, d$anchor)
      fit$row <- cbind(key, seed = d$seed, fit$row)
      fit$components <- cbind(case = case, replicate = replicate, fit$components)
      # On confirmation the two decisions are saved without reading d$truth or any labels.
      if (!is.null(choice)) {
        fit$decisions <- list(selected = decide(fit$row, choice$rule, choice$threshold),
                              lrt = decide(fit$row, "lrt"))
      }
      if (file.exists(halt_path)) stop("Another dataset failed.")
      save_atomic(fit, paths[i])
      TRUE
    }, error = function(e) {
      msg <- sprintf("case %d replicate %d: %s", case, replicate, conditionMessage(e))
      if (!file.exists(halt_path)) writeLines(msg, halt_path)
      stop(msg, call. = FALSE)
    })
  }

  done <- nrow(keys) - length(pending)
  cat(sprintf("%s: %d/%d datasets verified; %d workers\n", phase, done, nrow(keys), workers))
  if (!length(pending)) return(invisible(NULL))
  if (workers == 1L || .Platform$OS.type == "windows") {
    for (i in pending) {
      one(i)
      done <- done + 1L
      if (done %% 25L == 0L) cat(sprintf("%s %d/%d\n", phase, done, nrow(keys)))
    }
    return(invisible(NULL))
  }
  jobs <- list()
  on.exit(for (j in jobs) try(tools::pskill(j$pid), silent = TRUE), add = TRUE)
  while (length(pending) || length(jobs)) {
    if (file.exists(halt_path)) stop(paste(readLines(halt_path), collapse = "\n"))
    while (length(pending) && length(jobs) < workers) {
      i <- pending[1L]
      pending <- pending[-1L]
      child <- parallel::mcparallel(one(i), mc.set.seed = FALSE, silent = TRUE)
      jobs[[as.character(child$pid)]] <- child
    }
    finished <- parallel::mccollect(jobs, wait = FALSE, timeout = 1)
    for (pid in names(finished)) {
      jobs[[pid]] <- NULL
      if (!identical(finished[[pid]], TRUE)) stop("Worker failed: ", paste(finished[[pid]], collapse = " "))
      done <- done + 1L
      if (done %% 25L == 0L) cat(sprintf("%s %d/%d\n", phase, done, nrow(keys)))
    }
  }
  stopifnot(all(file.exists(paths)))
}

args <- commandArgs(TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Usage: Rscript run.R development|confirmation [workers]")
run_phase(args[1L], if (length(args) == 2L) as.integer(args[2L]) else 1L)
