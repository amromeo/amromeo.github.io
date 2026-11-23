source("tools/disagreement-v2/study.R")

# Decision denominators always use the original diagnostic estimate, including after withholding.
metrics <- function(error, decisions) {
  declined <- decisions$declined
  reported <- error[!declined]
  over10 <- error > 10
  over25 <- error > 25
  under5 <- error < 5
  stopifnot(all(is.finite(error)), !anyNA(declined))
  quant <- function(p) if (length(reported)) unname(quantile(reported, p, type = 7)) else NA_real_
  data.frame(attempted = length(error), declined = sum(declined), reported = length(reported),
    technical_withheld = sum(decisions$technical_withholding), lrt_unavailable = sum(decisions$lrt_unavailable),
    declined_over10 = sum(declined & over10), n_over10 = sum(over10),
    declined_over25 = sum(declined & over25), n_over25 = sum(over25),
    declined_under5 = sum(declined & under5), n_under5 = sum(under5),
    reported_mean = if (length(reported)) mean(reported) else NA_real_,
    reported_median = quant(.5), reported_p90 = quant(.9),
    reported_max = if (length(reported)) max(reported) else NA_real_,
    reported_over10 = sum(reported > 10))
}

choose <- function(candidates) {
  eligible <- candidates[candidates$n_under5 > 0L &
    20L * candidates$declined_under5 <= candidates$n_under5, , drop = FALSE]
  if (!nrow(eligible)) return(eligible)
  order_index <- order(-eligible$declined_over10, -eligible$threshold, eligible$rule)
  eligible[order_index[1L], , drop = FALSE]
}

passes <- function(m) {
  if (m$n_over10 == 0L || m$n_under5 == 0L) return(NA)
  3L * m$declined_over10 >= 2L * m$n_over10 && 20L * m$declined_under5 <= m$n_under5
}

# Join the design and endpoint truth after decisions have been computed.
score <- function(rows) {
  s <- stress_design[match(rows$case, stress_design$case), ]
  truth <- t(vapply(s$family, function(f) stress_forward(qnorm(c(.025, .975)), f), numeric(2L)))
  e <- vapply(seq_len(nrow(rows)), function(i) {
    stress_metrics(c(rows$diagnostic_lower[i], rows$diagnostic_upper[i]), truth[i, ], s$family[i])$width_error
  }, numeric(1L))
  cbind(rows, s[, c("family", "overlap", "affected"), drop = FALSE],
        true_overlap = s$disease_inside_true_ri, true_lower = truth[, 1L], true_upper = truth[, 2L], error = e)
}

read_phase <- function(phase) {
  dir <- file.path(work_dir, phase)
  keys <- keys_for_phase(phase)
  if (file.exists(file.path(dir, "HALT.txt"))) stop(paste(readLines(file.path(dir, "HALT.txt")), collapse = "\n"))
  paths <- checkpoint(dir, keys$case, keys$replicate)
  if (!all(file.exists(paths))) stop(sprintf("%s has %d/%d completed fits", phase, sum(file.exists(paths)), length(paths)))
  items <- lapply(paths, readRDS)
  rows <- do.call(rbind, lapply(items, "[[", "row"))
  par <- do.call(rbind, lapply(items, "[[", "components"))
  stopifnot(identical(as.integer(rows$case), keys$case),
            identical(as.integer(rows$replicate), keys$replicate), nrow(par) == sum(rows$components))
  list(items = items, rows = rows, components = par)
}

write_csv <- function(x, name) {
  path <- file.path(results_dir, name)
  if (grepl("-(datasets|components|decisions-before-scoring)\\.csv$", name)) {
    connection <- gzfile(paste0(path, ".gz"), open = "wt", compression = 9L)
    on.exit(close(connection), add = TRUE)
    write.csv(x, connection, row.names = FALSE, na = "")
  } else write.csv(x, path, row.names = FALSE, na = "")
}
read_csv <- function(name) read.csv(file.path(results_dir, name), stringsAsFactors = FALSE)
format_number <- function(x, digits = 1L) ifelse(is.finite(x),
  formatC(x, format = "f", digits = digits), "undefined")
fraction <- function(n, d) {
  ifelse(d == 0L, paste0(n, "/", d, " (undefined)"),
         sprintf("%d/%d (%.1f%%)", n, d, 100 * n / d))
}
report_table <- function(headers, rows) {
  render <- function(items) paste0("| ", paste(items, collapse = " | "), " |")
  c(render(headers), render(rep("---", length(headers))),
    apply(as.matrix(rows), 1L, render), "")
}

update_protocol <- function(selected) {
  protocol_path <- file.path(study_dir, "PROTOCOL.md")
  p <- readLines(protocol_path)
  status <- grep("^\\*\\*Status\\.", p)
  heading <- match("**Development choice to commit before confirmation**", p)
  stopifnot(length(status) == 1L, !is.na(heading), startsWith(p[heading + 4L], "| "))
  if (selected$status == "selected") {
    p[status] <- "**Status.** Development complete. The rule below is fixed for confirmation."
    p[heading + 4L] <- sprintf("| %s | %.2f | %d/%d | %d/%d |", selected$rule, selected$threshold,
      selected$declined_over10, selected$n_over10, selected$declined_under5, selected$n_under5)
  } else {
    p[status] <- "**Status.** Development complete. No candidate met the prespecified 5% cap."
    p[heading + 4L] <- "| No eligible rule | None | See candidate table | All 14 exceeded the cap |"
  }
  writeLines(p, protocol_path)
}

summarize_phase <- function(phase) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  selection <- NULL
  if (phase == "confirmation") selection <- locked_choice()
  fits <- read_phase(phase)
  rows <- fits$rows
  candidates <- NULL
  if (phase == "development") {
    candidate_ids <- expand.grid(rule = c("a", "b"), threshold = thresholds,
                                 stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    all_decisions <- lapply(seq_len(nrow(candidate_ids)), function(i)
      decide(rows, candidate_ids$rule[i], candidate_ids$threshold[i]))
    decisions <- list(lrt = decide(rows, "lrt"))
    scored <- score(rows)
    historical <- read.csv("assets/disagreement/datasets.csv")
    key <- function(x) paste(x$case, x$replicate, sep = "/")
    idx <- match(key(rows), key(historical))
    stopifnot(nrow(historical) == nrow(rows), !anyDuplicated(key(historical)), !anyNA(idx))
    prior <- historical[idx, ]
    # Same LRT function. Original Step B used archived estimates from a separate R runtime.
    write_csv(data.frame(datasets = nrow(rows),
      decisions_agree = sum(prior$secondary_declined == decisions$lrt$declined),
      max_abs_p_difference = max(abs(prior$lrt_p - rows$lrt_p), na.rm = TRUE),
      max_abs_lower_difference = max(abs(prior$diagnostic_anchored_lower - rows$diagnostic_lower)),
      max_abs_upper_difference = max(abs(prior$diagnostic_anchored_upper - rows$diagnostic_upper))),
      "development-archived-step-b-check.csv")
    candidates <- cbind(candidate_ids, do.call(rbind, lapply(all_decisions, function(x) metrics(scored$error, x))))
    candidates$eligible <- candidates$n_under5 > 0L & 20L * candidates$declined_under5 <= candidates$n_under5
    winner <- choose(candidates)
    if (nrow(winner)) {
      selected <- data.frame(status = "selected", winner[, c("rule", "threshold", "declined_over10",
        "n_over10", "declined_under5", "n_under5")])
      j <- which(candidate_ids$rule == selected$rule & candidate_ids$threshold == selected$threshold)
      decisions <- c(list(selected = all_decisions[[j]]), decisions)
    } else {
      selected <- data.frame(status = "no_eligible_rule", rule = NA_character_, threshold = NA_real_,
        declined_over10 = NA_integer_, n_over10 = candidates$n_over10[1L],
        declined_under5 = NA_integer_, n_under5 = candidates$n_under5[1L])
    }
    write_csv(candidates, "development-candidates.csv")
    write_csv(selected, "selected-rule.csv")
    update_protocol(selected)
  } else {
    decisions <- lapply(c("selected", "lrt"), function(method)
      do.call(rbind, lapply(fits$items, function(item) item$decisions[[method]])))
    names(decisions) <- c("selected", "lrt")
    stopifnot(all(vapply(decisions, nrow, integer(1L)) == nrow(rows)))
    # Save the complete paired decisions before joining endpoint truth.
    blind <- do.call(rbind, lapply(names(decisions), function(method)
      cbind(method = method, rows, decisions[[method]])))
    write_csv(blind, "confirmation-decisions-before-scoring.csv")
    scored <- score(rows)
  }
  output <- function(suffix) paste0(phase, "-", suffix, ".csv")
  methods <- names(decisions)
  for (method in methods) {
    d <- decisions[[method]]
    stopifnot(all(is.na(d$reported_lower) == d$declined), all(is.na(d$reported_upper) == d$declined),
      all(d$reported_lower[!d$declined] == rows$diagnostic_lower[!d$declined]),
      all(d$reported_upper[!d$declined] == rows$diagnostic_upper[!d$declined]))
  }
  full <- do.call(rbind, lapply(methods, function(m) cbind(method = m, scored, decisions[[m]])))
  summary <- do.call(rbind, lapply(methods, function(m) cbind(method = m,
    one_component = sum(rows$components == 1L), fit_failures = 0L,
    max_abs_parity_difference = max(rows$lower_abs_difference, rows$upper_abs_difference, rows$offset_abs_difference),
    metrics(scored$error, decisions[[m]]))))
  cell_rows <- lapply(split(seq_len(nrow(scored)), scored$case), function(idx) {
    s <- scored[idx, , drop = FALSE]
    o <- s$O[is.finite(s$O)]
    do.call(rbind, lapply(methods, function(m) cbind(
      s[1L, c("case", "family", "overlap", "affected", "true_overlap")], method = m,
      one_component = sum(s$components == 1L), O_defined = length(o),
      mean_O = if (length(o)) mean(o) else NA_real_, metrics(s$error, decisions[[m]][idx, ]))))
  })
  cells <- do.call(rbind, cell_rows)
  stopifnot(nrow(summary) == length(methods), nrow(cells) == 45L * length(methods),
    all(cells$attempted == 20L), all(cells$one_component + cells$O_defined == cells$attempted))
  write_csv(full, output("datasets"))
  write_csv(fits$components, output("components"))
  write_csv(summary, output("metrics"))
  write_csv(cells, output("cells"))
  if (phase == "confirmation") writeLines(selection$commit, file.path(results_dir, "selection-commit.txt"))
  writeLines(c(paste("Phase", phase), capture.output(sessionInfo())), file.path(results_dir, paste0(phase, "-runtime.txt")))
}

phase_section <- function(phase) {
  name <- function(x) paste0(phase, "-", x, ".csv")
  s <- read_csv(name("metrics"))
  cells <- read_csv(name("cells"))
  selected <- read_csv("selected-rule.csv")
  label <- function(m) ifelse(m == "selected", sprintf("Rule %s, t = %.2f", selected$rule, selected$threshold), "LRT alone")
  desc <- sprintf("All %d datasets passed the frozen-fit comparison. The largest absolute difference in either limit or the offset was %s. BIC selected one component in %d/%d datasets.",
                  s$attempted[1L], format(s$max_abs_parity_difference[1L], digits = 5), s$one_component[1L], s$attempted[1L])
  out <- c(paste("##", tools::toTitleCase(phase)), "", desc, "",
    report_table(c("Decision", "All declined", "Declined among E > 10%", "Declined among E > 25%", "Declined among E < 5%", "Technical withholdings"),
      data.frame(label(s$method), fraction(s$declined, s$attempted), fraction(s$declined_over10, s$n_over10),
        fraction(s$declined_over25, s$n_over25), fraction(s$declined_under5, s$n_under5), s$technical_withheld)),
    "Error among reported intervals, as a percentage of the true interval width", "",
    report_table(c("Decision", "Number reported", "Mean", "Median", "90th percentile", "Maximum", "E > 10%"),
      data.frame(label(s$method), s$reported, format_number(s$reported_mean, 2), format_number(s$reported_median, 2),
        format_number(s$reported_p90, 2), format_number(s$reported_max, 2), fraction(s$reported_over10, s$reported))))
  if (phase == "development") {
    old <- read_csv("development-archived-step-b-check.csv")
    out <- c(out, sprintf("The current LRT-alone decision agreed with archived Step B in %d/%d datasets. The archive used the earlier fitted limits.",
      old$decisions_agree, old$datasets), "")
  }
  families <- c(normal = "Normal", lognormal = "Lognormal", shifted = "Shifted lognormal")
  overlap_levels <- c("separated", "moderate", "heavy")
  affected <- c(.1, .3, .5, .7, .9)
  overlap_summary <- do.call(rbind, lapply(overlap_levels, function(ov) {
    x <- cells[cells$overlap == ov & cells$method == "lrt", ]
    keep <- x$O_defined > 0L
    estimate <- if (any(keep)) weighted.mean(x$mean_O[keep], x$O_defined[keep]) else NA_real_
    c(tools::toTitleCase(ov), sprintf("%.1f%%", 100 * x$true_overlap[1L]),
      if (is.finite(estimate)) sprintf("%.1f%%", 100 * estimate) else "undefined",
      sprintf("%d/%d", sum(x$O_defined), sum(x$attempted)))
  }))
  out <- c(out, "Fitted overlap by design level, pooled over all datasets with a defined O.", "",
    report_table(c("Overlap", "True", "Mean fitted O", "Defined / attempted"), overlap_summary))
  for (f in names(families)) {
    fc <- cells[cells$family == f, ]
    tr <- do.call(rbind, lapply(overlap_levels, function(ov) {
      x <- fc[fc$overlap == ov & fc$method == "lrt", ]
      x <- x[match(affected, x$affected), ]
      c(sprintf("%s (%.1f%% true)", tools::toTitleCase(ov), 100 * x$true_overlap[1L]),
        ifelse(is.na(x$mean_O), "undefined (0/20)", sprintf("%.1f%% (%d/20)", 100 * x$mean_O, x$O_defined)))
    }))
    out <- c(out, paste("###", families[[f]]), "",
      "Mean fitted O, with the number of defined estimates out of 20 in parentheses. The remaining fits selected one component.", "",
      report_table(c("Overlap", "10% affected", "30%", "50%", "70%", "90%"), tr))
    both <- "selected" %in% fc$method
    dec <- do.call(rbind, lapply(overlap_levels, function(ov) {
      x <- fc[fc$overlap == ov & fc$method == "lrt", ]
      x <- x[match(affected, x$affected), ]
      v <- fraction(x$declined, x$attempted)
      if (both) {
        y <- fc[fc$overlap == ov & fc$method == "selected", ]
        y <- y[match(affected, y$affected), ]
        v <- paste(fraction(y$declined, y$attempted), v, sep = " / ")
      }
      c(tools::toTitleCase(ov), v)
    }))
    out <- c(out, if (both) "Declines per cell, selected rule / LRT alone." else "Declines per cell, LRT alone.", "",
       report_table(c("Overlap", "10% affected", "30%", "50%", "70%", "90%"), dec))
  }
  out
}

write_report <- function(phase) {
  selected <- read_csv("selected-rule.csv")
  out <- c("# Routine-data separation study", "")
  if (phase == "confirmation") {
    s <- subset(read_csv("confirmation-metrics.csv"), method == "selected")
    verdict <- passes(s)
    text <- if (is.na(verdict)) "could not be evaluated" else if (verdict) "met both criteria" else "did not meet both criteria"
    out <- c(out, sprintf("**Confirmation.** Rule %s at t = %.2f %s. It declined %s of estimates with error over 10%% and %s of estimates with error under 5%%.",
      selected$rule, selected$threshold, text, fraction(s$declined_over10, s$n_over10), fraction(s$declined_under5, s$n_under5)),
      "", "The prespecified targets were at least two-thirds and at most 5%, respectively.", "", phase_section("confirmation"))
  } else if (selected$status == "selected") {
    out <- c(out, sprintf("Development selected rule %s at t = %.2f. Confirmation is pending.", selected$rule, selected$threshold), "")
  } else out <- c(out, "**Development stopped before confirmation.** None of the 14 candidates met the 5% cap.", "")
  out <- c(out, phase_section("development"), "### Development candidates", "")
  cand <- read_csv("development-candidates.csv")
  out <- c(out, report_table(c("Rule", "Threshold", "Declined among E > 10%", "Declined among E < 5%", "Eligible"),
    data.frame(cand$rule, sprintf("%.2f", cand$threshold), fraction(cand$declined_over10, cand$n_over10),
      fraction(cand$declined_under5, cand$n_under5), ifelse(cand$eligible, "Yes", "No"))),
    "The endpoint error E is the larger absolute error in the two limits divided by the true interval width, expressed as a percentage. The unit of analysis is one dataset. The 15 cases within each family and replicate share a screened sample.", "",
    "The dataset, component, cell and summary tables are in [results](results/). The [protocol](PROTOCOL.md) specifies the rules, withheld-interval messages and success criteria.",
    "", "### Reproducing the study", "",
    "From the anchoredRI package directory, run the commands below with base R. The worker count is optional. Per-dataset checkpoints are written beneath this study directory and can be resumed.", "",
    "```sh", "Rscript tools/disagreement-v2/test.R", "Rscript tools/disagreement-v2/run.R development 9",
    "Rscript tools/disagreement-v2/report.R development", "```", "",
    "If development selects a rule, commit the updated protocol, selected-rule.csv and study code before opening confirmation. On the published selection commit, run the following.", "",
    "```sh", "Rscript tools/disagreement-v2/run.R confirmation 9", "Rscript tools/disagreement-v2/report.R confirmation", "```")
  writeLines(out, file.path(study_dir, "REPORT.md"))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 1L || !args %in% c("development", "confirmation"))
    stop("Usage: Rscript report.R development|confirmation")
  summarize_phase(args[1L])
  write_report(args[1L])
}
