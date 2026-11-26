# SPDX-License-Identifier: GPL-3.0-only
# From the repository root: Rscript anchoredRI/tools/subsample/pilot.R INPUT_DIR OUTPUT_DIR INDEX
# Each process fits one member of the frozen twenty-dataset pilot.
.h_dir <- "anchoredRI/tools/subsample"
source(file.path(.h_dir, "inputs.R"))
source(file.path(.h_dir, "mixture.R"))
source(file.path(.h_dir, "aggregate.R"))

h_clock <- function() {
  t <- proc.time()
  c(wall_seconds = unname(t[3]), cpu_seconds = unname(sum(t[c(1, 2, 4, 5)])))
}

h_measure <- function(work) {
  start <- h_clock()
  value <- force(work)
  list(value = value, time = h_clock() - start)
}

# The fallback installs its own on.exit RNG cleanup. Measure its entry, then
# the public return, which also includes assignment of the reporting fields.
h_observe_fallback <- function(public_call, where) {
  timer <- new.env(parent = emptyenv())
  old <- options(taskH.fallback.timer = timer)
  on.exit(options(old), add = TRUE)
  suppressMessages(trace(".fit_anchor_refine", where = where, print = FALSE,
    tracer = quote(assign("start", proc.time(), envir = getOption("taskH.fallback.timer")))))
  on.exit(suppressMessages(untrace(".fit_anchor_refine", where = where)), add = TRUE)
  fit <- force(public_call)
  returned <- proc.time()
  spent <- if (is.null(timer$start)) numeric() else {
    delta <- returned - timer$start
    c(wall_seconds = unname(delta[3]), cpu_seconds = unname(sum(delta[c(1, 2, 4, 5)])))
  }
  stopifnot((length(spent) > 0L) == !is.null(fit$diagnostics$anchorRefine))
  list(fits = list("5" = fit), fallback_time = spent)
}

h_public <- function(case, pit_seed) {
  h_observe_fallback(
    anchoredRI::fit_anchored_ri(case$value, rep("pooled", length(case$value)), case$anchor,
      transform = "auto", safeguard = TRUE, probs = case$probs, rounding_unit = case$rounding_unit,
      fit_seed = 1L, fallback = TRUE, fallback_seed = 123L, pit_seed = pit_seed),
    asNamespace("anchoredRI"))
}

h_csv <- function(x, path) {
  # Keep IEEE doubles round-trippable in text exports as well as in the RDS.
  for (name in names(x)[vapply(x, is.double, logical(1))]) {
    z <- x[[name]]
    x[[name]] <- ifelse(is.na(z), NA_character_, sprintf("%.17g", z))
  }
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

h_pilot <- function(input_dir, output_dir, Index, approved_ids = NULL) {
  begin <- h_clock()
  preflight <- h_measure({
    stopifnot(as.character(getRversion()) == "4.5.0",
              as.character(utils::packageVersion("anchoredRI")) == "0.4.0")
    if (is.null(approved_ids)) {
      roster <- read.csv(file.path(.h_dir, "pilot-roster.csv"))
      stopifnot(nrow(roster) == 20L, !anyDuplicated(roster$Index))
      approved_ids <- roster$Index
    }
    stopifnot(Index %in% approved_ids)
    case <- h_read_case(input_dir, Index)
    n <- length(case$value)
    stopifnot(n >= 5000L, length(case$anchor) == 40L)
    pit_seed <- 38000000L + Index
    routine_hash <- h_vector_hash(case$value)
    anchor_hash <- h_vector_hash(case$anchor)
    target <- file.path(output_dir, as.character(Index))
    checkpoint_dir <- file.path(target, "attempts")
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    sources <- file.path(.h_dir, c("pilot.R", "aggregate.R", "inputs.R", "mixture.R"))
    version_id <- list(Index = Index, routine_hash = routine_hash, anchor_hash = anchor_hash,
                       code = vapply(sources, function(file) digest::digest(file = file, algo = "sha256"),
                                     character(1)), runtime = R.version.string,
                       anchoredRI = as.character(utils::packageVersion("anchoredRI")),
                       container = Sys.getenv("TASK_H_CONTAINER_DIGEST"))
    identity_file <- file.path(target, "input-id.rds")
    if (file.exists(identity_file)) stopifnot(identical(version_id, readRDS(identity_file))) else
      saveRDS(version_id, identity_file)
    if (file.exists(file.path(target, "COMPLETE"))) stop("Pilot output already complete for Index ", Index)
    TRUE
  })

  # Publish only complete checkpoints. A retry reconstructs a missing draw
  # from its own frozen seed and compares saved work with the input/source ID above.
  restored <- 0L
  retained_time <- c(wall_seconds = 0, cpu_seconds = 0)
  attempt <- function(key, work) {
    path <- file.path(checkpoint_dir, paste0(key, ".rds"))
    if (file.exists(path)) {
      item <- readRDS(path)
      stopifnot(identical(item$key, key))
      restored <<- restored + 1L
      retained_time <<- retained_time + item$time
      return(item)
    }
    item <- h_measure(work)
    item$key <- key
    tmp <- tempfile(tmpdir = checkpoint_dir)
    on.exit(unlink(tmp))
    saveRDS(item, tmp, compress = FALSE)
    stopifnot(file.rename(tmp, path))
    item
  }

  timing <- list()
  add_time <- function(phase, time, level = "phase", bank = NA_character_, draw = NA_integer_) {
    timing[[length(timing) + 1L]] <<- data.frame(Index = Index, level = level, phase = phase,
      seed_set = bank, draw = draw, wall_seconds = unname(time[1]), cpu_seconds = unname(time[2]))
  }
  add_time("preflight", preflight$time)
  fits <- models <- draws <- positions <- list()
  stash <- function(fitted, type, weight, bank = NA_character_, draw = NA_integer_,
                    seed = NA_integer_, pos_hash = NA_character_) {
    observed <- if (type == "subsample") 5000L else n
    for (cap_name in names(fitted)) {
      fit <- fitted[[cap_name]]
      row <- h_fit_row(fit, bic_n = round(weight * observed + length(case$anchor)))
      key <- data.frame(Index = Index, cap = as.integer(cap_name), fit_type = type,
                        seed_set = bank, draw = draw, fit_seed = 1L, sampling_seed = seed,
                        routine_hash = routine_hash, anchor_hash = anchor_hash, position_hash = pos_hash,
                        routine_n = observed, routine_weight = weight, rounding_unit = case$rounding_unit)
      fits[[length(fits) + 1L]] <<- cbind(key, row)
      models[[length(models) + 1L]] <<- list(key = key, fit = fit)
    }
  }

  std <- attempt("full_public", h_public(case, pit_seed))
  add_time("full_public", std$time)
  if (length(std$value$fallback_time)) add_time("fallback", std$value$fallback_time, "detail")
  stash(std$value$fits, "full", 1)
  search <- function(value, weight, caps) {
    h_mixture(value, case$anchor, case$rounding_unit, probs = case$probs, caps = caps,
              routine_weight = weight, fit_seed = 1L, fallback = FALSE, pit_seed = pit_seed)$fits
  }
  cap3 <- attempt("full_cap3", search(case$value, 1, 3L))
  add_time("full_cap3", cap3$time)
  stash(cap3$value, "full", 1)
  weight <- 5000 / n
  w <- attempt("W", search(case$value, weight, c(3L, 5L)))
  add_time("W", w$time)
  stash(w$value, "W", weight)

  for (bank in c("A", "B")) {
    message("Index ", Index, ", starting seed set ", bank)
    bank_start <- h_clock()
    previous <- retained_time
    for (b in seq_len(50L)) {
      label <- sprintf("%s%02d", bank, b)
      item <- attempt(label, {
        sample <- h_measure({
          indices <- h_sample_positions(Index, n, bank, b)
          list(indices = indices, hash = h_vector_hash(indices), seed = h_sample_seed(Index, bank, b))
        })
        model <- h_measure(search(case$value[sample$value$indices], 1, c(3L, 5L)))
        list(sample = sample$value, fits = model$value, sample_time = sample$time, fit_time = model$time)
      })
      z <- item$value
      add_time("sample", z$sample_time, "detail", bank, b)
      add_time("subsample_fit", z$fit_time, "detail", bank, b)
      positions[[label]] <- z$sample$indices
      draws[[label]] <- data.frame(Index = Index, seed_set = bank, draw = b, sampling_seed = z$sample$seed,
                                   position_hash = z$sample$hash, positions_n = length(z$sample$indices))
      stash(z$fits, "subsample", 1, bank, b, z$sample$seed, z$sample$hash)
      if (b %% 10L == 0L) message("Index ", Index, ", set ", bank, ": ", b, "/50 draws complete")
    }
    # On resume, include saved fit time once alongside the current bank I/O.
    add_time(paste0("bank_", bank), h_clock() - bank_start + retained_time - previous)
  }
  agg <- h_measure({
    fit_rows <- do.call(rbind, fits)
    draw_rows <- do.call(rbind, draws)
    stopifnot(nrow(fit_rows) == 204L, nrow(draw_rows) == 100L,
              !anyDuplicated(fit_rows[c("cap", "fit_type", "seed_set", "draw")]),
              !anyDuplicated(draw_rows[c("seed_set", "draw")]))
    h_summarize(fit_rows, case$probs)
  })
  add_time("aggregation", agg$time)
  output <- h_measure({
    saveRDS(models, file.path(target, "models.rds"), compress = "gzip")
    saveRDS(positions, file.path(target, "positions.rds"), compress = "gzip")
    h_csv(fit_rows, file.path(target, "fits.csv"))
    h_csv(draw_rows, file.path(target, "draws.csv"))
    h_csv(agg$value$medians, file.path(target, "medians.csv"))
    h_csv(agg$value$decisions, file.path(target, "decisions.csv"))
  })
  add_time("output", output$time)
  add_time("dataset", h_clock() - begin + retained_time, "total")
  times <- do.call(rbind, timing)
  times$restored_checkpoints <- restored
  h_csv(times, file.path(target, "timing.csv"))
  # The final bundle contains the complete fits and positions, so drop the duplicates.
  unlink(checkpoint_dir, recursive = TRUE)
  writeLines(sprintf("Index %d; 204 fits, 100 draws", Index), file.path(target, "COMPLETE"))
  invisible(target)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: pilot.R INPUT_DIR OUTPUT_DIR INDEX")
  h_pilot(args[[1L]], args[[2L]], as.integer(args[[3L]]))
}
