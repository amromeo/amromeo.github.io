# SPDX-License-Identifier: GPL-3.0-only
# One authorized slot of a full-run block, executed serially by its workflow.
source("anchoredRI/tools/subsample/pilot.R")

h_full_slot <- function(input_dir, block, slot) {
  stopifnot(length(block) == 1L, block %in% 1:178, length(slot) == 1L, slot %in% 1:8)
  frozen <- read.csv(file.path(.h_dir, "full-roster.csv"))
  pilot <- read.csv(file.path(.h_dir, "pilot-roster.csv"))
  original <- read.csv(file.path(input_dir, "task-g-index.csv"))
  stopifnot(nrow(original) == 1440L, nrow(pilot) == 20L, all(pilot$Index %in% original$Index),
            !anyDuplicated(original$Index), !anyDuplicated(pilot$Index))
  remaining <- sort(setdiff(original$Index, pilot$Index))
  pos <- seq_along(remaining) - 1L
  source_rows <- original[match(remaining, original$Index), ]
  stopifnot(length(remaining) == 1420L, identical(frozen$Index, remaining),
            all(frozen$block == pos %/% 8L + 1L), all(frozen$slot == pos %% 8L + 1L),
            identical(frozen$routine_hash, source_rows$routine_hash),
            identical(frozen$anchor_hash, source_rows$anchor_hash),
            identical(frozen$shard, source_rows$shard), identical(frozen$routine_n, source_rows$routine_n))
  row <- frozen[frozen$block == block & frozen$slot == slot, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  row
}

h_full_receipt <- function(path, Index) {
  receipt <- read.csv(file.path(path, "case-receipt.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(names(receipt), c("Index", "state", "path", "bytes", "sha256")),
            nrow(receipt) >= 1L, !anyNA(receipt), all(receipt$Index == Index),
            !anyDuplicated(receipt$path), length(unique(receipt$state)) == 1L)
  complete <- file.exists(file.path(path, "COMPLETE"))
  if (complete) {
    expected <- c("COMPLETE", "decisions.csv", "draws.csv", "fits.csv", "input-id.rds",
                  "medians.csv", "models.rds", "positions.rds", "timing.csv")
    stopifnot(identical(readLines(file.path(path, "COMPLETE")), sprintf("Index %d; 204 fits, 100 draws", Index)))
  } else {
    attempts <- list.files(file.path(path, "attempts"), pattern = "[.]rds$")
    stopifnot(all(grepl("^(full_public|full_cap3|W|[AB](0[1-9]|[1-4][0-9]|50))[.]rds$", attempts)))
    expected <- c("input-id.rds", if (length(attempts)) paste0("attempts/", attempts))
  }
  stopifnot(setequal(receipt$path, expected), all(receipt$state == if (complete) "complete" else "partial"))
  files <- file.path(path, receipt$path)
  stopifnot(all(file.info(files)$size == receipt$bytes), all(receipt$bytes > 0),
            all(vapply(files, function(x) digest::digest(file = x, algo = "sha256"),
                       character(1), USE.NAMES = FALSE) == receipt$sha256))
  identity <- readRDS(file.path(path, "input-id.rds"))
  stopifnot(identity$Index == Index)
  invisible(list(complete = complete, identity = identity))
}

h_full_worker <- function(input_dir, output_dir, block, slot) {
  row <- h_full_slot(input_dir, block, slot)
  Index <- row$Index
  path <- file.path(output_dir, as.character(Index))
  receipt <- file.path(path, "case-receipt.csv")
  if (file.exists(receipt)) {
    saved <- h_full_receipt(path, Index)
    id <- saved$identity
    stopifnot(identical(id$routine_hash, row$routine_hash), identical(id$anchor_hash, row$anchor_hash))
    if (saved$complete) {
      sources <- file.path(.h_dir, c("pilot.R", "aggregate.R", "inputs.R", "mixture.R"))
      current <- vapply(sources, function(x) digest::digest(file = x, algo = "sha256"), character(1))
      stopifnot(identical(id$code, current), identical(id$runtime, R.version.string),
                identical(id$anchoredRI, as.character(utils::packageVersion("anchoredRI"))),
                identical(id$container, Sys.getenv("TASK_H_CONTAINER_DIGEST")))
      message("Reusing verified completed Index ", Index)
      return(invisible(path))
    }
    message("Resuming verified attempts for Index ", Index)
  } else if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("Restored payload has no verified receipt for Index ", Index)
  }
  h_pilot(input_dir, output_dir, Index, approved_ids = row$Index)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) == 3L && identical(args[[1L]], "--verify")) {
    h_full_receipt(args[[2L]], as.integer(args[[3L]]))
  } else {
    if (length(args) != 4L || any(!grepl("^[1-9][0-9]*$", args[3:4])))
      stop("Usage: full-worker.R INPUT_DIR OUTPUT_DIR BLOCK SLOT, or --verify CASE_DIR INDEX")
    h_full_worker(args[[1L]], args[[2L]], as.integer(args[[3L]]), as.integer(args[[4L]]))
  }
}
