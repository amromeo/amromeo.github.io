# SPDX-License-Identifier: GPL-3.0-only
# Only these measurements and frozen hashes are exposed to the fitting worker.
confirmation_read_case <- local({
    bundle_path <- NULL
    source_env <- NULL
    design <- NULL
    reference <- NULL
    function(input_dir, j) {
        path <- normalizePath(input_dir, mustWork = TRUE)
        if (!identical(path, bundle_path)) {
            stopifnot(as.character(getRversion()) == "4.5.0")
            manifest <- read.table(file.path(path, "INPUT-SHA256SUMS"), col.names = c("hash", "file"),
                                   colClasses = "character", stringsAsFactors = FALSE)
            expected <- c("archive-SHA256SUMS", "confirmation-design.R", "confirmation-design.csv",
                          "confirmation-generator.R", "generator.R", "input-sources.json", "reference-hashes.csv",
                          "references.rds")
            stopifnot(length(manifest$file) == length(expected), setequal(manifest$file, expected))
            for (i in seq_len(nrow(manifest))) {
                stopifnot(identical(digest::digest(file = file.path(path, manifest$file[i]), algo = "sha256"),
                                    manifest$hash[i]))
            }
            source_env <<- new.env(parent = parent.env(globalenv()))
            for (file in c("generator.R", "confirmation-design.R", "confirmation-generator.R")) {
                sys.source(file.path(path, file), envir = source_env)
            }
            # Read the original 17-digit CSV with its original R column classes.
            template <- source_env$stage2_confirmation_design()
            original <- read.csv(file.path(path, "confirmation-design.csv"),
                                 stringsAsFactors = FALSE, colClasses = vapply(template, typeof, character(1)))
            stopifnot(identical(original$j, seq_len(2880L)),
                      identical(names(original), names(template)))
            design <<- original[original$N == 50000L, , drop = FALSE]
            reference <<- readRDS(file.path(path, "references.rds"))
            stopifnot(identical(as.character(design$j), names(reference)), nrow(design) == 1440L)
            bundle_path <<- path
        }
        j <- as.integer(j)
        stopifnot(length(j) == 1L, !is.na(j), j %in% design$j)
        row <- design[match(j, design$j), , drop = FALSE]
        archived <- reference[[as.character(j)]]
        regenerated <- source_env$stage2_confirmation_generate(row)
        routine_hash <- digest::digest(regenerated$x, algo = "sha256", serializeVersion = 2)
        anchor_hash <- digest::digest(regenerated$anchors, algo = "sha256", serializeVersion = 2)
        stopifnot(identical(archived$j, j), identical(routine_hash, archived$routine_hash),
                  identical(anchor_hash, archived$anchor_hash), identical(regenerated$anchors, archived$anchor),
                  length(regenerated$x) == 50000L, length(archived$anchor) == 40L,
                  is.finite(row$h), row$h > 0, identical(row$h, regenerated$h))
        list(value = regenerated$x, anchor = archived$anchor, rounding_unit = row$h,
             routine_hash = routine_hash, anchor_hash = anchor_hash)
    }
})
