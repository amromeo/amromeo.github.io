# SPDX-License-Identifier: GPL-3.0-only
# Score final endpoints, including out-of-support values, with the pinned Stage 2 function.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, as.character(getRversion()) == "4.5.0")
source(args[[1]])
x <- read.csv(args[[2]], stringsAsFactors = FALSE, na.strings = "")
x$z_lower <- x$z_upper <- rep(NA_real_, nrow(x))
for (fam in unique(x$family)) {
    i <- which(x$family == fam)
    pars <- unique(x[i, c("lambda", "mu", "sigma")])
    stopifnot(nrow(pars) == 1L)
    eval_z <- function(v) stage2_healthy_z(v, pars$lambda, pars$mu, pars$sigma)
    truth <- unique(x[i, c("truth_lower", "truth_upper")])
    stopifnot(nrow(truth) == 1L, max(abs(eval_z(unlist(truth)) - qnorm(c(.025,.975)))) < 1e-12)
    x$z_lower[i] <- qnorm(.025) - eval_z(x$lower[i])
    x$z_upper[i] <- qnorm(.975) - eval_z(x$upper[i])
}
x$d <- (abs(x$z_lower) + abs(x$z_upper)) / 2
for (col in c("z_lower", "z_upper", "d")) {
    v <- x[[col]]
    x[[col]] <- ifelse(is.na(v), NA_character_, sprintf("%.17g", v))
}
write.csv(x[c("j", "bank", "estimator", "z_lower", "z_upper", "d")], args[[3]], row.names = FALSE, na = "")
