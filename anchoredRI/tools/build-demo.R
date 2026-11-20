# Run from the package source directory: Rscript tools/build-demo.R [workers]
# The same routine draw is reused for every method and both screening rules.
stopifnot(file.exists("DESCRIPTION"), requireNamespace("refineR", quietly = TRUE))
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(path)
args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L
stopifnot(length(workers) == 1L, is.finite(workers), workers >= 1L)
screened_sizes <- c(20L, 40L, 80L, 120L, 400L)
reference_sizes <- c(120L, 400L)
screening_rules <- c("clinical", "in_range")
grid <- expand.grid(scenario = 1:4, replicate = 1:12, n_routine = c(1000L, 4000L))

one_dataset <- function(i) {
  s <- grid$scenario[i]
  k <- grid$replicate[i]
  n <- grid$n_routine[i]
  clinic <- n %/% 4L
  d <- simulate_ri_study(s, k, clinic, n - clinic, 0L)
  rows <- compare_ri_methods(d$routine$value, d$routine$setting, include_stock = TRUE)
  rows$n_anchor <- 0L
  rows$selection <- "any"
  for (rule in screening_rules) {
    anchor <- simulate_ri_study(s, k, clinic, n - clinic, max(screened_sizes), rule)$anchor
    for (size in screened_sizes) {
      screened <- screened_comparisons(d$routine$value, d$routine$setting, head(anchor, size))
      screened$n_anchor <- size
      screened$selection <- rule
      rows <- rbind(rows, screened)
    }
    for (size in reference_sizes) {
      sample <- head(anchor, size)
      limits <- mean(sample) + qnorm(c(.025, .975)) * sd(sample)
      reference <- interval_rows(paste0("screened_", size), rbind(limits))
      reference$n_anchor <- size
      reference$selection <- rule
      rows <- rbind(rows, reference)
    }
  }
  rows$scenario <- s
  rows$replicate <- k
  rows$n_routine <- n
  if (any(!is.finite(rows$lower) | !is.finite(rows$upper) | rows$lower >= rows$upper)) {
    stop(sprintf("Incomplete fit at scenario %d, draw %d, n=%d", s, k, n))
  }
  rows
}
indices <- seq_len(nrow(grid))
runs <- if (.Platform$OS.type == "unix" && workers > 1L) {
  parallel::mclapply(indices, one_dataset, mc.cores = workers, mc.set.seed = FALSE)
} else {
  lapply(indices, one_dataset)
}
stopifnot(all(vapply(runs, is.data.frame, logical(1))))
results <- do.call(rbind, runs)
rows_per_dataset <- 4L + length(screening_rules) * (3L * length(screened_sizes) + length(reference_sizes))
stopifnot(nrow(results) == nrow(grid) * rows_per_dataset)
results$enrichment_gap <- NULL
for (column in c("lower", "upper")) results[[column]] <- round(results[[column]], 8)
write.table(results, "assets/fits.csv", row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = ",", na = "")
packages <- c("refineR", "ash", "future", "future.apply", "parallelly", "globals", "listenv", "digest")
versions <- vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))
writeLines(c(R.version.string, sprintf("%s %s", packages, versions),
             "96 synthetic datasets; 12 replicates per case mix and routine sample size.",
             "Screened sample sizes: 20, 40, 80, 120, 400. Fixed screened-only comparisons: 120 and 400.",
             "All intervals are point estimates on the identity input scale."),
           "assets/versions.txt")
