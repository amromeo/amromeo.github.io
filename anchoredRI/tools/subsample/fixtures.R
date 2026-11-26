# SPDX-License-Identifier: GPL-3.0-only
# From the repository root, after the study freeze: Rscript --vanilla
# anchoredRI/tools/subsample/fixtures.R <fresh-output-directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
  stop("Fixture output directory must be empty. Keep prior attempts.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
study <- "anchoredRI/tools/subsample"
options(digits = 17)
for (file in c("inputs.R", "mixture.R", "aggregate.R")) source(file.path(study, file))
pkg <- function(name) getFromNamespace(name, "anchoredRI")
host <- data.frame(boot_id = readLines("/proc/sys/kernel/random/boot_id", n = 1L),
                   host = Sys.info()[["nodename"]], pid = Sys.getpid(),
                   invocation = readLines("/proc/sys/kernel/random/uuid", n = 1L),
                   job = Sys.getenv("GITHUB_JOB"), run = Sys.getenv("GITHUB_RUN_ID"),
                   attempt = Sys.getenv("GITHUB_RUN_ATTEMPT"), R = R.version.string,
                   package = as.character(packageVersion("anchoredRI")))
write.csv(host, file.path(out, "host.csv"), row.names = FALSE)
checks <- fits <- list()
science <- list()
check <- function(name, ok, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(test = name, pass = isTRUE(ok),
                                              detail = as.character(detail))
}
section <- function(label, expr) {
  tryCatch(force(expr), error = function(e) check(label, FALSE, conditionMessage(e)))
}
# Include all candidate scores, fitted parameters, safeguards and estimates in
# the exact-repeat evidence. Timing and host metadata are separate.
keep_fit <- function(label, fit, probs) {
  sg <- fit$safeguard
  raw <- fit$diagnostics$anchored$limits
  pick <- function(x, key) if (key %in% names(x)) unname(x[[key]]) else NA_real_
  available <- length(probs) == 1L || is.finite(pick(raw, "lower"))
  check(paste(label, "usable first stage"), available && is.finite(pick(raw, "upper")) &&
          is.finite(sg$lrt_p) && fit$anchoredRI_status %in% c("reported", "declined"))
  row <- data.frame(case = label, shape = if (length(probs) == 1L) "upper" else "two",
                    transform = fit$transform, components = fit$components,
                    first_status = fit$anchoredRI_status, final_status = fit$status,
                    method = if (is.null(fit$method)) NA_character_ else fit$method,
                    lrt_flag = sg$lrt_flag, lower_flag = sg$lower_flag,
                    lower_skipped = sg$lower_check_skipped, lrt_p = sg$lrt_p,
                    share = fit$weights[1L, 1L], shift = fit$shift,
                    raw_lower = pick(raw, "lower"), raw_upper = pick(raw, "upper"),
                    lower = pick(fit$limits, "lower"), upper = pick(fit$limits, "upper"))
  exported <- h_fit_row(fit)
  for (field in c("raw_lower", "raw_upper", "lower", "upper", "share", "lrt_flag", "lower_flag")) {
    check(paste(label, "export", field), identical(exported[[field]], row[[field]]))
  }
  check(paste(label, "export first-stage p"), identical(exported$p, row$lrt_p))
  fits[[label]] <<- row
  science[[label]] <<- fit[c("reference", "shift", "weights", "selection", "candidate_limits",
                              "candidate_at_bound", "candidate_shift", "bic_original", "limits",
                              "transform", "components", "reason", "status", "anchoredRI_reason",
                              "anchoredRI_status", "method", "runner_up", "safeguard")]
}
same_fit <- function(tag, actual, expected) {
  # Inspect the pre- and post-fallback public records field by field. Include
  # details that are absent from h_fit_row, including nonwinning candidate BICs.
  fields <- c("reference", "transform", "components", "weights", "shift", "selection", "bic_original",
              "candidate_limits", "candidate_shift", "candidate_at_bound", "selected_bic",
              "runner_up", "limits", "safeguard", "status", "reason")
  for (field in fields)
    check(paste(tag, field), identical(actual[[field]], expected[[field]]))
  check(paste(tag, "raw limits"), identical(actual$diagnostics$anchored$limits,
                                            expected$diagnostics$anchored$limits))
}
normal_values <- function(n, mu, sigma) mu + sigma * qnorm((seq_len(n) - .5) / n)
# 125 distinct design positions repeated in order create ties on the reporting
# grid and 6,250 routine records (so W has a fractional routine count).
x <- rep(round(4 * c(normal_values(85L, 12, 1.4), normal_values(40L, 18, 2))) / 4, 50L)
a <- round(4 * normal_values(40L, 12, 1.4)) / 4
two <- c(.025, .975)
weight <- 5000 / length(x)
index <- 17L
pit <- 38000000L + index
call_fit <- function(v, probs, q, caps, fallback = FALSE)
  h_mixture(v, a, .25, probs = probs, caps = caps, routine_weight = q,
            fit_seed = 1L, fallback = fallback, pit_seed = pit, keep_candidates = TRUE)
public <- function(v, probs, fallback)
  anchoredRI::fit_anchored_ri(v, rep("pooled", length(v)), a, transform = "auto",
    rounding_unit = .25, safeguard = TRUE, fit_seed = 1L, fallback = fallback,
    fallback_seed = 123L, pit_seed = pit, probs = probs)

# Refit using the installed ordinary package candidate machinery and stop at
# three. This independent oracle never calls h_mixture or filters its results.
ordinary_capped <- function(v, probs) {
  d <- pkg(".arm_b_data")(v, rep("pooled", length(v)), a, "auto", NULL, .25, 1L)
  rules <- c(d$unit, .05 * d$spread, pkg(".log_sd_factor")(max(.005, .05 *
                pkg(".routine_scale")(log(d$value)))))
  families <- c("identity", "log", "shifted_log")
  models <- list()
  for (k in 1:3) for (family in families)
    models[[paste(family, k)]] <- pkg(".fit_mixture_candidate")(d, family, k, models, rules, 1L)
  order <- unlist(lapply(families, paste, 1:3), use.names = FALSE)
  fit <- pkg(".mixture_results")(models[order], d, "auto", d$unit, 1L)
  list(fit = pkg(".apply_ri_safeguard")(fit, v, a, TRUE, probs), candidates = models)
}
candidate_record <- function(m) m[c("ll", "bic", "p")]
compare_candidates <- function(tag, actual, expected) {
  check(paste(tag, "candidate IDs"), setequal(names(actual), names(expected)) && length(actual) == 9L)
  for (id in names(expected)) {
    check(paste(tag, id), identical(candidate_record(actual[[id]]), candidate_record(expected[[id]])))
  }
}
# Evaluate directly on each uncompressed observation, on the measurement scale.
# Rowwise log-sum-exp avoids underflow for very narrow candidate components.
density <- function(v, p, transform, k) {
  coord <- if (transform == "identity") v - 1 else log(v - p$shift)
  dnorm(coord, p$mu[k], p$sigma[k], log = TRUE) -
    if (transform == "identity") 0 else log(v - p$shift)
}
loglike <- function(v, m, q) {
  p <- m$p
  terms <- vapply(seq_len(m$k), function(k)
    log(p$w[k, 1L]) + density(v, p, m$transform, k), numeric(length(v)))
  highest <- apply(terms, 1L, max)
  q * sum(highest + log(rowSums(exp(terms - highest)))) + sum(density(a, p, m$transform, 1L))
}
direct_scores <- function(tag, result, v, q) {
  mass <- q * length(v) + length(a)
  check(paste(tag, "BIC total weight"), identical(as.double(result$bic_n), as.double(mass)))
  ids <- names(result$candidates)
  check(paste(tag, "fifteen candidates"), length(ids) == 15L)
  for (id in ids) {
    m <- result$candidates[[id]]
    ll <- loglike(v, m, q)
    bic <- -2 * ll + (3 * m$k - 1 + (m$transform == "shifted_log")) * log(mass)
    for (field in c("ll", "bic")) {
      diff <- abs(m[[field]] - if (field == "ll") ll else bic)
      check(paste(tag, id, field, "direct"), is.finite(diff) && diff <= 1e-7, diff)
    }
  }
  science[[paste(tag, "candidate parameters")]] <<- lapply(result$candidates, candidate_record)
}

section("sampler fixture", {
  samples <- list()
  for (bank in c("A", "B")) for (draw in c(1L, 50L)) {
    base <- if (bank == "A") 100000000L else 200000000L
    seed <- base + 100L * index + draw
    pos <- h_sample_positions(index, length(x), bank, draw)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    expected <- sort(sample.int(length(x), 5000L, replace = FALSE))
    key <- paste(bank, draw)
    check(paste(key, "sampler"), identical(as.integer(h_sample_seed(index, bank, draw)), seed) &&
            identical(pos, expected) && length(unique(pos)) == 5000L)
    samples[[key]] <- data.frame(bank = bank, draw = draw, seed = seed, position = pos)
  }
  write.csv(do.call(rbind, samples), file.path(out, "positions.csv"), row.names = FALSE)
})
section("two-sided unit-weight fits", {
  result <- call_fit(x, two, 1, c(3L, 5L), TRUE)
  ref <- public(x, two, TRUE)
  same_fit("public fit vs cap 5", result$fits[["5"]], ref)
  for (field in c("method", "anchoredRI_status", "anchoredRI_reason"))
    check(paste("public reporting", field), identical(result$fits[["5"]][[field]], ref[[field]]))
  direct_scores("unit weights", result, x, 1)
  # Public selection retains every original-scale BIC, including nonwinners.
  # Recover the ordinary full-search likelihood from it for all 15 candidates.
  for (i in seq_len(nrow(ref$selection))) {
    row <- ref$selection[i, ]
    id <- paste(row$transform, row$components)
    p <- 3L * row$components - 1L + (row$transform == "shifted_log")
    ref_ll <- (p * log(length(x) + 40L) - row$bic_original) / 2
    diff <- abs(result$candidates[[id]]$ll - ref_ll)
    check(paste("ordinary likelihood", id), is.finite(diff) && diff <= 1e-9, diff)
  }
  capped <- ordinary_capped(x, two)
  ids <- names(capped$candidates)
  compare_candidates("unit weights capped reuse", result$candidates[ids], capped$candidates)
  same_fit("ordinary independent cap 3", result$fits[["3"]], capped$fit)
  for (cap in names(result$fits)) keep_fit(paste("full", cap), result$fits[[cap]], two)
})
section("weighted fits", {
  result <- call_fit(x, two, weight, c(3L, 5L))
  direct_scores("W", result, x, weight)
  combined <- c(x, a)
  weights <- c(rep(weight, length(x)), rep(1, length(a)))
  mu <- sum(weights * combined) / sum(weights)
  sigma <- sqrt(sum(weights * (combined - mu)^2) / sum(weights))
  expected <- c(mu, sigma)
  candidate <- result$candidates[["identity 1"]]$p
  diff <- abs(c(candidate$mu + 1, candidate$sigma) - expected)
  for (i in seq_along(diff)) check(paste("Gaussian analytic", c("mean", "SD")[i]),
                                   is.finite(diff[i]) && diff[i] <= 1e-6, diff[i])
  # Independent capped run uses the same W objective on fresh candidate starts.
  capped <- call_fit(x, two, weight, 3L)
  ids <- names(capped$candidates)
  compare_candidates("weighted capped reuse", result$candidates[ids], capped$candidates)
  same_fit("W independent cap 3", result$fits[["3"]], capped$fits[["3"]])
  for (cap in names(result$fits)) keep_fit(paste("W", cap), result$fits[[cap]], two)
})
section("upper-only subsample fit", {
  v <- x[h_sample_positions(index, length(x), "A", 1L)]
  result <- call_fit(v, .95, 1, c(3L, 5L))
  ref <- public(v, .95, FALSE)
  same_fit("public CRP first stage", result$fits[["5"]], ref)
  for (cap in names(result$fits)) {
    fit <- result$fits[[cap]]
    check(paste("CRP skips lower check", cap), isTRUE(fit$safeguard$lower_check_skipped) &&
            identical(fit$safeguard$lower_flag, FALSE) && identical(names(fit$limits), "upper"))
    keep_fit(paste("CRP subsample", cap), fit, .95)
  }
})

# Saved rows exercise release policies independently of whether an optimizer
# happens to generate borderline vote counts on these particular fixtures.
section("aggregation and safeguards", {
  draws <- data.frame(raw_lower = 1:50, raw_upper = (1:50) + 100,
                      share = seq(.4, .89, by = .01), components = 2L,
                      first_status = rep(c("reported", "declined"), each = 25))
  all <- h_median(draws, two)
  pass <- h_median(draws, two, passing_only = TRUE)
  check("median of 50", all$contributors == 50 && all$lower == 25.5 && all$upper == 125.5)
  check("passing-only median", pass$contributors == 25 && pass$lower == 13 && pass$upper == 113)
  # An invalid member drops its paired endpoint as well.
  bad <- draws
  bad$raw_lower[50] <- 999
  m <- h_median(bad, two)
  check("invalid paired endpoints", m$contributors == 49 && m$lower == 25 && m$upper == 125)
  crp <- draws
  crp$raw_lower <- NA_real_
  crp$raw_upper[1L] <- NA_real_
  m <- h_median(crp, .95)
  check("upper-only median", m$contributors == 49 && is.na(m$lower) && m$upper == 126)
  draws$first_status <- "declined"
  none <- h_median(draws, two, passing_only = TRUE)
  check("zero passing fits", none$contributors == 0 && is.na(none$lower) && is.na(none$upper))
  check("unavailable median overrides statistical decline",
        identical(h_gate("declined", none$lower, none$upper, two), "withheld"))
  majority <- function(v, u) h_majority(c(rep("declined", v), rep("withheld", u),
                                          rep("reported", 50 - v - u)))
  for (row in list(c(25, 0, "reported"), c(26, 0, "declined"), c(25, 1, "withheld"),
                   c(26, 1, "declined"), c(24, 1, "reported"), c(0, 50, "withheld"))) {
    check(paste("majority declines/unresolved", row[1], row[2]),
          identical(majority(as.integer(row[1]), as.integer(row[2])), row[3]))
  }
  sg <- pkg(".ri_flags")(c(lower = 2, upper = 4), 20, 3)
  check("both statistical flags can coexist", isTRUE(sg$lower_flag) && isTRUE(sg$lrt_flag))
  saved <- data.frame(status = "withheld", method = NA_character_, first_status = "declined",
                      lower = NA_real_, upper = NA_real_, reason = "withheld fixture")
  route <- h_route("declined", 7, 11, saved, two)
  check("lookup withheld package despite usable raw median",
        route$route == "full_nonreturn" && route$status == "withheld" && is.na(route$lower) && is.na(route$upper) &&
          route$reason == "withheld fixture")
  # Actual per-Index summaries must always look up the standard full package.
  # On cap 3 a deliberately different saved full record would release limits.
  key_rows <- function(cap) {
    frame <- data.frame(Index = 17L, cap = cap,
                        fit_type = c("full", "W", rep("subsample", 100)),
                        seed_set = c("", "", rep(c("A", "B"), each = 50)),
                        draw = c(NA_integer_, NA_integer_, rep(1:50, 2)), components = 2L,
                        raw_lower = 7, raw_upper = 11, share = .6,
                        first_status = "declined", first_reason = "flag", status = "declined",
                        method = NA_character_, lower = NA_real_, upper = NA_real_, reason = "flag")
    if (cap == 5L) {
      frame$status[1L] <- "withheld"
      frame$reason[1L] <- "withheld fixture"
    } else {
      frame$status[1L] <- "reported"
      frame$method[1L] <- "anchoredRI"
      frame$lower[1L] <- 3
      frame$upper[1L] <- 4
    }
    frame
  }
  policies <- h_summarize(rbind(key_rows(5L), key_rows(3L)), two)
  routed <- subset(policies$decisions, rule == "majority_or_package")
  check("both caps and both definitions look up the withheld standard full record",
        nrow(routed) == 8L && all(routed$route == "full_nonreturn") &&
          all(routed$status == "withheld") && all(routed$reason == "withheld fixture") &&
          all(is.na(routed$lower)) && all(is.na(routed$upper)))
  science$aggregation <- list(all = all, passing = pass, none = none, route = route, routed = routed)
})

summary <- do.call(rbind, checks)
write.csv(summary, file.path(out, "checks.csv"), row.names = FALSE, na = "")
if (length(fits)) write.csv(do.call(rbind, fits), file.path(out, "fits.csv"), row.names = FALSE, na = "")
saveRDS(science, file.path(out, "science.rds"), version = 2, compress = FALSE)
if (any(!summary$pass)) {
  print(summary[!summary$pass, ], row.names = FALSE)
  stop(sum(!summary$pass), " fixture checks failed. See checks.csv.")
}
