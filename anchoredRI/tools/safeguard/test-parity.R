# Compare installed public fits with the verified Task E fitter on the same R runtime.
# Arguments are the Task E anchoredRI directory, the small-N RIbench directory
# the shared object built with tools/boxcox/reliability-v3/build.R and the large input shards.
library(anchoredRI)
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
study <- normalizePath(args[1])
small_data <- normalizePath(args[2])
object <- normalizePath(args[3])
shards <- normalizePath(args[4])
setwd(study)
source("tools/boxcox/reliability-v3/common.R")
source("tools/stress-design.R")
dyn.load(object)
meta <- read.csv("tools/boxcox/anchorRefine-comparator/anchoredRI-comparator.csv")
defs <- loadTestsetDefinition()
keys <- c("r-1", "r-7", "r-150", "r-4609", "r-4755", "r-870", "s-01-01")
cols <- c("dataset_id", "transform", "components", "declined")
stored <- function(x) x[x$method == "armB" & x$dataset_id %in% keys, cols]
parts <- list.files("tools/boxcox/accuracy/results/diagnostics", "^part-", full.names = TRUE)
historical <- do.call(rbind, lapply(parts, function(f) stored(read.csv(f))))
historical <- rbind(historical, stored(read.csv("tools/boxcox/accuracy/large/results/armB-fits.csv")))
stopifnot(!anyDuplicated(historical$dataset_id), setequal(historical$dataset_id, keys))
for (key in keys) {
  row <- meta[meta$dataset_id == key, ]
  stopifnot(nrow(row) == 1L)
  if (row$stage == "grid") {
    d <- stress_data(row$case, row$replicate)
    seed <- d$seed
    digits <- 2L
    probs <- c(.025, .975)
  } else {
    s <- defs[defs$Index == row$Index, ]
    if (s$N == 50000L) {
      position <- match(s$Index, defs$Index[defs$N == 50000L])
      file <- file.path(shards, sprintf("inputs-%02d.rds", (position - 1L) %% 24L + 1L))
      inputs <- Filter(function(x) x$s$Index == s$Index, readRDS(file))
      stopifnot(length(inputs) == 1L)
      d <- inputs[[1L]]
      d$setting <- rep("pooled", length(d$value))
    } else {
      d <- bc_ribench_input(s, small_data, row$anchor_hash)
    }
    seed <- s$Index
    digits <- s$decimals
    probs <- d$probs
  }
  stopifnot(identical(digest::digest(d$value, algo = "sha256"), row$routine_hash),
            identical(digest::digest(d$anchor, algo = "sha256"), row$anchor_hash))
  data <- v2_data(d$value, d$setting, d$anchor)
  models <- r3_baseline(data, seed, digits, 5L)
  table <- v2_table(models, probs)
  winner <- table[table$selected, ]
  old <- models[[paste(winner$transform, winner$k)]]
  p <- old$p
  if (old$transform == "identity") p$mu <- p$mu + 1
  candidate <- vector("list", old$k)
  candidate[[old$k]] <- p
  archived <- bc_safeguard(list(transform = old$transform, components = old$k,
       shift = p$shift, candidate_parameters = candidate), d$value, d$anchor, probs)
  before <- .Random.seed
  new <- fit_anchored_ri(d$value, d$setting, d$anchor, fit_seed = seed, probs = probs)
  close <- function(x, y) max(abs(as.numeric(x) - as.numeric(y))) <= 1e-10
  stopifnot(identical(before, .Random.seed), new$requested_transform == "auto",
            new$rounding_unit == row$rounding_step, identical(new$transform, winner$transform),
            new$components == winner$k, close(new$weights, p$w), close(new$shift, p$shift),
            identical(new$selection$transform, table$transform),
            identical(new$selection$components, table$k),
            close(new$selection$bic_original, table$bic),
            close(new$diagnostics$anchored$limits, archived$limits),
            close(new$safeguard$lrt_statistic, archived$lrt_statistic),
            identical(new$safeguard$lower_flag, archived$lower_flag),
            identical(new$safeguard$lrt_flag, archived$lrt_p < .01),
            identical(new$declined, archived$declined),
            identical(new$technical_withholding, archived$technical))
  h <- historical[historical$dataset_id == key, ]
  cat(sprintf("%s: package and study source agree (%s, K=%d, %s); archive (%s, K=%d, %s)\n",
              key, new$transform, new$components, new$status, h$transform, h$components,
              if (h$declined) "declined" else "reported"))
}
