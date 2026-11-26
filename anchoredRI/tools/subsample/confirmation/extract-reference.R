# SPDX-License-Identifier: GPL-3.0-only
# Export only the Stage 2 routine hash, guidance hash and original 40 values.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, as.character(getRversion()) == "4.5.0")
shards <- args[1]
out <- args[2]
design <- read.csv(file.path(out, "confirmation-design.csv"))
ids <- design$j[design$N == 50000L]
stopifnot(length(ids) == 1440L, !anyDuplicated(ids))
records <- setNames(vector("list", length(ids)), as.character(ids))
for (shard in seq_len(48L)) {
    archived <- readRDS(file.path(shards, sprintf("%d.rds", shard)))
    expected <- seq.int(shard, 2880L, by = 48L)
    # Do not expose scores, parameters or decisions to the input adapter.
    id_fields <- unique(archived$fits[c("j", "routine_hash", "anchor_hash")])
    stored_ids <- vapply(archived$cases, "[[", integer(1), "j")
    stopifnot(identical(stored_ids, expected), identical(as.integer(id_fields$j), expected))
    for (i in seq_along(expected)) {
        j <- expected[i]
        if (!j %in% ids) next
        case <- archived$cases[[i]]
        panel <- case$panels$guidance
        routine_hash <- id_fields$routine_hash[i]
        anchor_hash <- id_fields$anchor_hash[i]
        stopifnot(identical(case$routine_hash, routine_hash),
                  grepl("^[0-9a-f]{64}$", routine_hash), length(panel) == 40L,
                  is.double(panel), all(is.finite(panel)), is.null(records[[as.character(j)]]),
                  identical(digest::digest(panel, algo = "sha256", serializeVersion = 2), anchor_hash))
        records[[as.character(j)]] <- list(j = j, routine_hash = routine_hash,
                                           anchor_hash = anchor_hash, anchor = panel)
    }
}
stopifnot(all(vapply(records, length, integer(1)) == 4L))
saveRDS(records, file.path(out, "references.rds"), version = 2)
hashes <- data.frame(j = ids, routine_hash = unname(vapply(records, "[[", character(1), "routine_hash")),
                     anchor_hash = unname(vapply(records, "[[", character(1), "anchor_hash")))
write.table(hashes, file.path(out, "reference-hashes.csv"), sep = ",", quote = FALSE, row.names = FALSE)
