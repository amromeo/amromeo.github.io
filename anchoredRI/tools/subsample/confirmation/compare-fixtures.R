# SPDX-License-Identifier: GPL-3.0-only
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, args[1] %in% c("same", "cross"))
same <- args[1] == "same"
left <- args[2]
right <- args[3]
read <- function(dir, name) read.csv(file.path(dir, name), stringsAsFactors = FALSE,
                                     na.strings = "", check.names = FALSE)
ah <- read(left, "host.csv")
bh <- read(right, "host.csv")
stopifnot(nrow(ah) == 1L, nrow(bh) == 1L, !is.na(ah$boot), !is.na(bh$boot),
          (ah$boot == bh$boot) == same, ah$invocation != bh$invocation,
          ah$R == bh$R, ah$anchoredRI == bh$anchoredRI)
files <- c("positions.rds", "routes.csv")
if (same) files <- c(files, "fits.csv", "science.rds")
for (f in files) stopifnot(identical(digest::digest(file = file.path(left, f), algo = "sha256"),
                                      digest::digest(file = file.path(right, f), algo = "sha256")))
a <- read(left, "fits.csv")
b <- read(right, "fits.csv")
stopifnot(identical(a$fixture, c("full", "A01", "full_on_A01")), identical(a$fixture, b$fixture))
decisions <- c("components", "transform", "status", "first_status", "method", "reason", "first_reason",
               "lrt_flag", "lower_flag", "lower_check_skipped", "fallback_status", "fallback_reason")
stopifnot(identical(a[decisions], b[decisions]))
max_end <- max_p <- 0
for (field in c("raw_lower", "raw_upper", "lower", "upper", "p", "fallback_p")) {
  x <- a[[field]]
  y <- b[[field]]
  stopifnot(identical(is.na(x), is.na(y)))
  valid <- which(!is.na(x))
  if (!length(valid)) next
  delta <- abs(x[valid] - y[valid])
  if (field %in% c("p", "fallback_p")) {
    stopifnot(all(is.finite(delta)), all(delta <= 0.0001))
    max_p <- max(max_p, delta)
  } else {
    if (field %in% c("lower", "upper")) width <- a$upper[valid] - a$lower[valid]
    else width <- a$raw_upper[valid] - a$raw_lower[valid]
    stopifnot(all(is.finite(width)), all(width > 0))
    fraction <- delta / width
    stopifnot(all(is.finite(fraction)), all(fraction <= 0.001))
    max_end <- max(max_end, fraction)
  }
}
cat(sprintf("%s host comparison passed: boot IDs %s, %s; max endpoint difference/width %.17g; max p difference %.17g\n",
            args[1], ah$boot, bh$boot, max_end, max_p))
