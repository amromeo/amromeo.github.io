# Ordered vectors, rounding and requested percentiles from the primary Task G records.
h_vector_hash <- function(x) digest::digest(x, algo = "sha256", serializeVersion = 2L)

h_read_case <- local({
  last_dir <- NULL
  last_shard <- NULL
  stored <- NULL
  index <- NULL
  design <- NULL

  function(input_dir, Index) {
    dir <- normalizePath(input_dir, mustWork = TRUE)
    if (!identical(last_dir, dir)) {
      read_table <- function(name) read.csv(file.path(dir, name), stringsAsFactors = FALSE)
      index <<- read_table("task-g-index.csv")
      design <<- read_table("design-large.csv")
      stopifnot(!anyDuplicated(index$Index), !anyDuplicated(design$Index),
                setequal(index$Index, design$Index))
      last_dir <<- dir
      last_shard <<- NULL
    }
    stopifnot(length(Index) == 1L, !is.na(Index), Index %in% index$Index)
    entry <- index[match(Index, index$Index), , drop = FALSE]
    meta <- design[match(Index, design$Index), , drop = FALSE]
    if (!identical(last_shard, entry$shard)) {
      stored <<- readRDS(file.path(dir, sprintf("taskg-shard-%d.rds", entry$shard)))$cases
      last_shard <<- entry$shard
    }
    case <- stored[[entry$rds_position]]
    value <- case$routine
    anchor <- case$panels$guidance
    stopifnot(case$Index == Index, case$j == Index, length(value) == meta$routine_n,
              length(anchor) == 40L, is.finite(meta$h), meta$h > 0, meta$anchor_seed == 17000000L + Index,
              identical(h_vector_hash(value), meta$routine_hash),
              identical(h_vector_hash(anchor), meta$anchor_hash),
              identical(meta$routine_hash, entry$routine_hash), identical(meta$anchor_hash, entry$anchor_hash),
              isTRUE(meta$requested_lower == (meta$Analyte != "CRP")))
    probs <- if (meta$requested_lower) c(.025, .975) else .95
    stopifnot(tail(probs, 1L) == meta$upper_probability)
    list(Index = as.integer(Index), value = value, anchor = anchor, rounding_unit = meta$h,
         probs = probs, meta = meta)
  }
})

h_sample_seed <- function(Index, bank, draw) {
  stopifnot(length(Index) == 1L, !is.na(Index), Index > 0, Index == floor(Index),
            length(bank) == 1L, bank %in% c("A", "B"),
            length(draw) == 1L, !is.na(draw), draw %in% seq_len(50L))
  as.integer(if (bank == "A") 100000000 + 100 * Index + draw else 200000000 + 100 * Index + draw)
}

h_sample_positions <- function(Index, n, bank, draw) {
  stopifnot(length(n) == 1L, !is.na(n), n >= 5000L, n == floor(n), n <= .Machine$integer.max)
  seed <- h_sample_seed(Index, bank, draw)
  # Preserve the caller's RNG state, including the absence of a seed.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  kinds <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kinds))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  sort(sample.int(as.integer(n), 5000L, replace = FALSE))
}
