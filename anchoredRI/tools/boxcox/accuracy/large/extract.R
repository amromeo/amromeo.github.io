args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
root <- args[1]
archive <- args[2]
tables <- file.path(root, "tables")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
source("tools/boxcox/reliability-v3/common.R")
source("tools/boxcox/post-round3/comparator.R")
defs <- subset(loadTestsetDefinition(), N == 50000L)
keys <- paste0("r-", defs$Index)
paths <- list.files(file.path(root, "fits"), "^r-[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE)
ids <- sub("[.]rds$", "", basename(paths))
stopifnot(length(paths) == 1440L, !anyDuplicated(ids), setequal(ids, keys))
paths <- paths[match(keys, ids)]
template <- do.call(rbind, readRDS(file.path(archive, "anchoredRI-native-performance.rds"))$resRIs)
custom_files <- list.files(archive, "^custom-[0-9]+[.]csv$", full.names = TRUE)
custom <- subset(do.call(rbind, lapply(custom_files, read.csv)), algo == "anchoredRI")
stopifnot(nrow(template) == 1440L, nrow(custom) == 1440L, !anyDuplicated(custom$Index),
          setequal(custom$Index, defs$Index), setequal(template$Index, defs$Index))
template <- template[match(defs$Index, template$Index), ]
custom <- custom[match(defs$Index, custom$Index), ]
tables_by_fit <- lapply(seq_along(paths), function(i) {
    z <- readRDS(paths[i])
    x <- z$pair
    s <- defs[i, ]
    one <- s$Analyte == "CRP"
    probs <- if (one) .95 else c(.025, .975)
    stopifnot(z$key == keys[i], all(x$Index == s$Index), identical(as.character(x$method), c("anchoredRI", "armB")),
              all(x$routine_hash == z$routine), all(x$anchor_hash == z$anchor),
              all(x$ribench_hashcode == s$HashCode), all(x$anchor_seed == custom$anchor_seed[i]),
              all(x$n_anchor == 40L), all(is.na(x$reported_upper) == (x$declined | x$technical)))
    diag <- x[, c("dataset_id", "stage", "Index", "method", "lower", "upper", "lrt_p", "lower_flag",
                  "declined", "technical", "transform", "components", "shift", "at_bound")]
    old <- z$frozen
    g <- old$guard
    frozen <- diag[1L, ]
    frozen$method <- "frozenActions"
    frozen$lower <- if (one) NA_real_ else g$limits[1]
    frozen$upper <- tail(g$limits, 1)
    frozen[c("transform", "components", "shift", "at_bound")] <- old[c("transform", "components", "shift", "at_bound")]
    frozen[c("lrt_p", "lower_flag", "declined", "technical")] <- g[c("lrt_p", "lower_flag", "declined", "technical")]
    original <- diag[1L, ]
    original$method <- "archivedC2"
    original$lower <- if (one) NA_real_ else template$LRL[i]
    original$upper <- template$URL[i]
    original$transform <- custom$model[i]
    original$components <- custom$components[i]
    original$shift <- custom$shift[i]
    original$at_bound <- custom$at_bound[i]
    original$lrt_p <- custom$lrt_p[i]
    original$lower_flag <- !one && is.finite(original$lower) && original$lower < x$routine_min[1]
    original$technical <- !is.finite(original$lrt_p) || !is.finite(original$upper) || (!one && !is.finite(original$lower))
    original$declined <- !original$technical && (original$lrt_p < .01 || original$lower_flag)

    cand <- if (length(z$models)) cbind(dataset_id = z$key, v2_table(z$models, probs)) else NULL
    if (!is.null(cand)) {
        cand$selected_baseline <- !is.na(x$components[1]) & cand$transform == x$transform[1] & cand$k == x$components[1]
        cand$selected_armB <- !is.na(x$components[2]) & cand$transform == x$transform[2] & cand$k == x$components[2]
    }
    params <- exits <- list()
    for (j in seq_len(nrow(x))) {
        if (is.na(x$components[j])) next
        m <- z$models[[paste(x$transform[j], x$components[j])]]
        p <- m$p
        stopifnot(nrow(p$w) == m$k, ncol(p$w) == 1L)
        params[[j]] <- data.frame(dataset_id = z$key, method = x$method[j], transform = m$transform,
            lambda = p$lambda, shift = p$shift, component = seq_len(m$k), screened_component = seq_len(m$k) == 1L,
            setting = "pooled", transformed_mean = p$mu + as.integer(m$transform == "identity"),
            transformed_sd = p$sigma, routine_weight = p$w[, 1])
        win <- m$starts[[m$origin]]
        line <- grepl("line[ -]?search|lnsrch", win$message, ignore.case = TRUE) || m$code == 52L
        iter <- m$code == 1L || grepl("iterat.*limit|limit.*iterat", win$message, ignore.case = TRUE)
        if (line || iter) exits[[j]] <- data.frame(dataset_id = z$key, method = x$method[j], transform = m$transform,
            components = m$k, code = m$code, origin = win$label, message = win$message, iteration_limit = iter,
            line_search = line, initial_retained = m$initial_retained, initial_ll = win$initial,
            final_ll = win$final, retained_ll = m$ll, declined = x$declined[j], reported = x$report_status[j] == "reported")
    }
    timing <- data.frame(dataset_id = z$key, analyte = s$Analyte, N = s$N, unique_values = z$unique_values,
        run = z$run, attempt = z$attempt, shard = z$shard, image = z$image, began = z$began, ended = z$ended,
        fit_seconds = z$fit_seconds, fit_cpu_seconds = z$cpu, frozen_seconds = z$frozen_seconds)
    list(baseline = x[1L, ], armB = x[2L, ], diagnostics = rbind(diag, frozen, original), candidates = cand,
         parameters = do.call(rbind, params), exits = do.call(rbind, exits), timings = timing)
})
collect <- function(field) do.call(rbind, lapply(tables_by_fit, `[[`, field))
write_tab <- function(x, name) {
    for (field in names(x)[vapply(x, is.double, logical(1))])
        x[[field]] <- ifelse(is.na(x[[field]]), NA_character_, sprintf("%.17g", x[[field]]))
    write.csv(x, file.path(tables, paste0(name, ".csv")), row.names = FALSE, na = "")
}
for (field in c("diagnostics", "candidates", "timings")) write_tab(collect(field), field)
write_tab(collect("parameters"), "selected-component-parameters")
exits <- collect("exits")
if (is.null(exits)) exits <- data.frame(dataset_id = character(), method = character(),
                                       iteration_limit = logical(), line_search = logical())
write_tab(exits, "winning-optimizer-exits")

for (method in c("baseline", "armB")) {
    dest <- file.path(tables, paste0(method, "-fits.csv"))
    write_comparator(lapply(tables_by_fit, `[[`, method), dest)
    expected <- collect(method)
    check <- read.csv(dest, na.strings = "", stringsAsFactors = FALSE)
    for (field in c("lower", "upper", "loglik", "bic", "lrt_p", "lrt_statistic", "reported_lower", "reported_upper"))
        stopifnot(identical(as.double(check[[field]]), as.double(expected[[field]])))
}
saveRDS(list(large_archive = normalizePath(archive)), file.path(root, "large-config.rds"))
