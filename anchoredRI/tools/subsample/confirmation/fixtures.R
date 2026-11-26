# SPDX-License-Identifier: GPL-3.0-only
# Fixed numerical and routing fixtures. Run in a separate R process on each host.
source("anchoredRI/tools/subsample/confirmation/worker.R")
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, as.character(getRversion()) == "4.5.0",
          as.character(packageVersion("anchoredRI")) == "0.4.0")
out <- args[1]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Fixture output exists")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
future::plan(future::sequential)
host <- data.frame(boot = readLines("/proc/sys/kernel/random/boot_id", n = 1L),
                   invocation = readLines("/proc/sys/kernel/random/uuid", n = 1L),
                   runner = Sys.getenv("RUNNER_NAME"), pid = Sys.getpid(), R = R.version.string,
                   anchoredRI = as.character(packageVersion("anchoredRI")))
h_csv(host, file.path(out, "host.csv"))

# Reuse the deterministic routine and panel from the Task H timing fixtures.
qn <- function(n, mu, sigma) mu + sigma * qnorm((seq_len(n) - .5) / n)
routine <- rep(round(4 * c(qn(85L, 12, 1.4), qn(40L, 18, 2))) / 4, 50L)
anchors <- round(4 * qn(40L, 12, 1.4)) / 4
j <- 17L
positions <- list()
for (bank in c("A", "B")) {
  for (draw in c(1L, 50L)) {
    key <- sprintf("%s%02d", bank, draw)
    position <- h_sample_positions(j, length(routine), bank, draw)
    seed <- (if (bank == "A") 100000000L else 200000000L) + 100L * j + draw
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    stopifnot(identical(seed, h_sample_seed(j, bank, draw)),
              identical(position, sort(sample.int(length(routine), 5000L, replace = FALSE))))
    positions[[key]] <- position
  }
}
saveRDS(positions, file.path(out, "positions.rds"), version = 2)
subset_x <- routine[positions$A01]
models <- list(full = confirmation_fit(routine, anchors, .25, j, TRUE),
               A01 = confirmation_fit(subset_x, anchors, .25, j, FALSE),
               full_on_A01 = confirmation_fit(subset_x, anchors, .25, j, TRUE))
rows <- list(full = confirmation_row(models$full, j, length(routine)),
             A01 = confirmation_row(models$A01, j, length(subset_x), "A", 1L,
                                    h_sample_seed(j, "A", 1L), h_vector_hash(positions$A01)),
             full_on_A01 = confirmation_row(models$full_on_A01, j, length(subset_x)))
stopifnot(rows$full$first_status == "reported", rows$A01$first_status == "declined",
          rows$full_on_A01$first_status == rows$A01$first_status,
          !is.na(rows$full_on_A01$fallback_status), is.null(models$A01$diagnostics$anchorRefine))
for (key in names(rows)) {
  r <- rows[[key]]
  stopifnot(h_usable(r$raw_lower, r$raw_upper, confirmation_probs), is.finite(r$p))
}
fit_table <- do.call(rbind, rows)
fit_table$fixture <- names(rows)
h_csv(fit_table, file.path(out, "fits.csv"))
# Include the full-precision first-stage selection and safeguard records, plus
# fallback fit/check outputs that determine the final interval.
stable_model <- function(fit) {
  fallback <- fit$diagnostics$anchorRefine
  list(fit = fit[c("reference", "weights", "selection", "bic_original", "candidate_limits",
                   "candidate_shift", "candidate_at_bound", "limits", "runner_up", "transform",
                   "components", "reason", "status", "anchoredRI_reason", "anchoredRI_status",
                   "method", "safeguard", "optimizer_status", "optimizer_code")],
       fallback = fallback[c("weight", "status", "reason", "candidate_limits", "limits", "check")])
}
saveRDS(lapply(models, stable_model), file.path(out, "science.rds"), version = 2)

routing_test <- function(name, first, package, v, t, step, status, expected,
                         usable = TRUE, omit_first = FALSE) {
  full <- rows$full
  full$first_status <- first
  full$status <- package
  full$reason <- full$fallback_reason <- ""
  full$lower <- if (package == "reported") 4 else NA_real_
  full$upper <- if (package == "reported") 12 else NA_real_
  draws <- data.frame(first_status = rep(c("declined", "withheld", "reported"), c(v,t,50L-v-t)),
                      raw_lower = as.double(1:50), raw_upper = 100 + as.double(1:50), share = 0.6)
  if (!usable) draws$raw_upper <- NA_real_
  if (omit_first) draws$raw_upper[1L] <- NA_real_
  m <- h_median(draws, confirmation_probs)
  votes <- h_votes(draws$first_status)
  route <- confirmation_route(full, m, votes, h_majority(draws$first_status))
  stopifnot(identical(as.numeric(votes), c(v, t, 50-v-t)), route$step == step,
            route$status == status, route$contributors == if (usable) 50L - omit_first else 0L,
            identical(c(route$lower, route$upper), as.numeric(expected)))
  cbind(fixture = name, route)
}
missing <- c(NA_real_, NA_real_)
m <- c(25.5, 125.5)
routes <- rbind(
  routing_test("a over contrary vote", "reported", "reported", 26, 0, "a", "reported", m),
  routing_test("b over approving vote", "declined", "reported", 0, 0, "b", "reported", c(4,12)),
  routing_test("b over declining vote", "declined", "reported", 50, 0, "b", "reported", c(4,12)),
  routing_test("c 25 versus 25", "declined", "declined", 25, 0, "c", "reported", m),
  routing_test("c robust technical vote", "withheld", "withheld", 24, 1, "c", "reported", m),
  routing_test("d 26 declines", "declined", "declined", 26, 0, "d", "declined", missing),
  routing_test("d unresolved majority", "declined", "declined", 25, 1, "d", "withheld", missing),
  routing_test("a unusable median", "reported", "reported", 0, 0, "a", "withheld", missing, FALSE),
  routing_test("d unusable median", "declined", "declined", 0, 0, "d", "withheld", missing, FALSE),
  routing_test("paired contributor omission", "reported", "reported", 10, 0, "a", "reported", c(26,126),
               omit_first = TRUE))
h_csv(routes, file.path(out, "routes.csv"))
cat("Verified 10 routing cases, three public fits and four draw position fixtures.\n")
