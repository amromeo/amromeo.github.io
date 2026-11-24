# Run from anchoredRI with the same RIbench release and reader as the accepted stage.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
archive_repo <- normalizePath(args[1])
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("tools/boxcox/reliability-v3/common.R")
stopifnot(getRversion() == "4.5.0", packageVersion("RIbench") == "1.0.2",
          packageVersion("digest") == "0.6.39", packageVersion("data.table") == "1.16.4")

large <- file.path(archive_repo, "stage1/large")
stage1 <- file.path(archive_repo, "results/stage1-45276858699/artifacts/large-inputs.zip")
scored <- file.path(large, "ribench-scored-independent-50000.zip")
stopifnot(digest::digest(file = scored, algo = "sha256") ==
              "15a5122d5e792aa947daea1bcb5fd8bf57e46f2b6f683c71c23f826fb71a2210",
          digest::digest(file = stage1, algo = "sha256") ==
              "57477e21f7fd4eefd4c42da064e570b16e344517449ce0270b17be600a7a85af")
hashes <- read.csv(unz(stage1, "large-routines.csv"))
anchors <- read.csv(unz(scored, "anchor-sample-hashes.csv"))
defs <- subset(loadTestsetDefinition(), N == 50000)
spec <- read.csv(file.path(large, "testsets.csv"))
stopifnot(nrow(defs) == 1440L, !anyDuplicated(defs$Index), all(table(defs$Analyte) == 144L),
          identical(defs$Index, spec$Index), identical(defs$HashCode, spec$HashCode),
          identical(defs$Index, hashes$Index), identical(defs$HashCode, hashes$routine_md5),
          setequal(defs$Index, anchors$Index), !anyDuplicated(anchors$Index))
anchors <- anchors[match(defs$Index, anchors$Index), ]
stopifnot(all(anchors$seed == 17000000L + defs$Index))
started <- Sys.time()
data_dir <- tempfile("task-e-routines-")
dir.create(data_dir)
RIbench::generateBiomarkerTestSets(workingDir = data_dir, subset = defs, rounding = TRUE, verbose = FALSE)
shards <- 24L
for (shard in seq_len(shards)) {
    positions <- seq(shard, nrow(defs), by = shards)
    data <- lapply(positions, function(i) {
        d <- bc_ribench_input(defs[i, ], data_dir, anchors$sha256[i])
        stopifnot(length(d$anchor) == 40L, length(d$value) == hashes$n_routine[i],
                  digest::digest(d$value, algo = "sha256") == hashes$routine_sha256[i])
        list(value = d$value, anchor = d$anchor, s = defs[i, ], probs = d$probs,
             routine_hash = hashes$routine_sha256[i], anchor_hash = anchors$sha256[i])
    })
    saveRDS(data, file.path(out, sprintf("inputs-%02d.rds", shard)))
}
write.csv(data.frame(started = as.numeric(started), finished = as.numeric(Sys.time()),
                     datasets = nrow(defs), shards), file.path(out, "generation.csv"), row.names = FALSE)
