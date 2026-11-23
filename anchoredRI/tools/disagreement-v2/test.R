source("tools/disagreement-v2/report.R")

# Known design overlaps and unequal service sizes. The healthy weight is excluded.
q <- qnorm(c(.025, .975))
inside <- function(mu, sigma) diff(pnorm(q, mu, sigma))
mixed <- list(components = 3L, component_mu = c(0, 4, .75), component_sd = c(1, 1.25, 1.25),
  weights = matrix(c(.8, .15, .05, .2, .30, .50), nrow = 3L,
                   dimnames = list(NULL, c("service_a", "service_b"))))
setting <- c(rep("service_b", 3), "service_a")
s <- component_summary(mixed, setting)
w <- .25 * mixed$weights[, 1L] + .75 * mixed$weights[, 2L]
expected <- (w[2L] * inside(4, 1.25) + w[3L] * inside(.75, 1.25)) / sum(w[2:3])
stopifnot(abs(s$O - expected) < 1e-14, max(abs(s$components$pooled_weight - w)) < 1e-14,
          abs(s$components$probability_inside[1L] - .95) < 1e-14)
for (i in seq_len(nrow(stress_design))) {
  design <- stress_design[i, ]
  stopifnot(abs(inside(design$delta, design$disease_sd_z) - design$disease_inside_true_ri) < 1e-14)
}
one <- list(components = 1L, component_mu = 2, component_sd = 1,
  weights = matrix(1, nrow = 1L, ncol = 2L, dimnames = list(NULL, c("service_a", "service_b"))))
stopifnot(is.na(component_summary(one, setting)$O))

# A mismatch beyond the absolute tolerance stops the dataset.
f <- list(limits = c(10, 20), components = 2L, shift = 0, transform = "log")
check_equivalence(f, f)
for (field in c("lower", "upper", "shift", "transform", "components")) {
  g <- f
  if (field == "lower") g$limits[1L] <- g$limits[1L] + 2e-10
  if (field == "upper") g$limits[2L] <- g$limits[2L] + 2e-10
  if (field == "shift") g$shift <- 2e-10
  if (field == "transform") g$transform <- "identity"
  if (field == "components") g$components <- 1L
  stopifnot(inherits(try(check_equivalence(f, g), silent = TRUE), "try-error"))
}
a <- seq_len(40)
own <- lrt_env$ao_normal(a)
stopifnot(lrt_env$lrt(a, own$mu + q * own$mle_sd, "identity", 0)[2L] > .9999999)

# Threshold equality, both flags, unavailable LRT, and K = 1 under both candidates.
rows <- data.frame(components = c(2L,2L,2L,2L,1L,1L,1L,2L),
  O = c(.60,.61,.59,.61,NA,NA,NA,.59), lrt_p = c(.01,.009,.009,NA,.01,.009,NA,NA),
  diagnostic_lower = 10, diagnostic_upper = 20)
aa <- decide(rows, "a", .60)
bb <- decide(rows, "b", .60)
ll <- decide(rows, "lrt")
stopifnot(identical(aa$declined, c(FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE)),
          identical(bb$declined, c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE)),
          identical(ll$declined, c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE)),
          bb$separation_flag[2L], bb$lrt_flag[2L],
          bb$message[2L] == "routine results are too mixed to separate healthy from affected patients",
          bb$message[3L] == "routine data disagree with the screened sample",
          !bb$technical_withholding[4L], all(bb$technical_withholding[c(7L,8L)]),
          identical(aa$declined[5:7], ll$declined[5:7]))
for (d in list(aa, bb, ll)) {
  stopifnot(identical(is.na(d$reported_lower), d$declined),
            identical(is.na(d$reported_upper), d$declined),
            identical(nzchar(d$recommendation), d$declined))
}
m <- metrics(c(5,10,25,26,0,11,4,100), aa)
stopifnot(m$n_under5 == 2L, m$declined_under5 == 1L, m$n_over10 == 4L,
          m$declined_over10 == 2L, m$n_over25 == 2L, m$declined_over25 == 1L,
          m$reported == 4L, m$reported_over10 == 2L,
          m$reported_p90 == unname(quantile(c(5,25,0,100), .9, type = 7)))
cand <- data.frame(rule = c("b","a","b","a"), threshold = c(.55,.70,.70,.80),
  n_under5 = 100L, declined_under5 = c(5L,5L,5L,6L), declined_over10 = c(20L,20L,20L,50L))
stopifnot(choose(cand)$rule == "a", choose(cand)$threshold == .70)
cand$declined_under5 <- 6L
stopifnot(nrow(choose(cand)) == 0L)
criterion <- data.frame(declined_over10 = 2L, n_over10 = 3L, declined_under5 = 1L, n_under5 = 20L)
stopifnot(passes(criterion))
criterion$n_under5 <- 19L
stopifnot(!passes(criterion))
criterion$n_over10 <- 0L
stopifnot(is.na(passes(criterion)))
cat("All separation, equivalence, decision, selection and reporting checks passed.\n")
