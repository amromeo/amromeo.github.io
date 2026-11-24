args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
inputs <- readRDS(args[1])
out <- args[2]
workers <- as.integer(args[3])
dyn.load(args[4])
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("tools/boxcox/reliability-v3/common.R")
source("tools/boxcox/post-round3/comparator.R")
expected <- read.csv("tools/boxcox/reliability-v3/results/source-hashes.csv")
stopifnot(getRversion() == "4.5.0", packageVersion("RIbench") == "1.0.2",
          packageVersion("data.table") == "1.16.4", packageVersion("digest") == "0.6.39")
stopifnot(all(vapply(expected$file, function(f) digest::digest(file = f, algo = "sha256"),
                     character(1)) == expected$sha256))
frozen <- bc_frozen(TRUE)
adapter <- bc_capacity(3L)
template <- read.csv("tools/boxcox/anchorRefine-comparator/anchoredRI-comparator.csv",
                     na.strings = "", nrows = 1)[NA_integer_, ]
manifest <- c(r3_sources(), driver = digest::digest(file = "tools/boxcox/accuracy/large/fit.R", algo = "sha256"),
              native = digest::digest(file = "tools/boxcox/reliability-v3/native-recovery.cpp", algo = "sha256"),
              guard_export = digest::digest(file = "tools/boxcox/post-round3/comparator.R", algo = "sha256"))
fit_attempt <- function(expr) tryCatch(list(value = force(expr), error = NULL),
    error = function(e) list(value = NULL, error = conditionMessage(e)))

run_one <- function(pos) {
    input <- inputs[[pos]]
    s <- input$s
    key <- paste0("r-", s$Index)
    path <- file.path(out, paste0(key, ".rds"))
    x <- input$value
    a <- input$anchor
    probs <- input$probs
    stopifnot(s$N == 50000L, length(a) == 40L,
              digest::digest(x, algo = "sha256") == input$routine_hash,
              digest::digest(a, algo = "sha256") == input$anchor_hash,
              digest::digest(x, algo = "md5", serializeVersion = 2, ascii = TRUE) == s$HashCode)
    if (file.exists(path)) {
        z <- readRDS(path)
        stopifnot(z$key == key, identical(z$manifest, manifest), z$routine == input$routine_hash,
                  z$anchor == input$anchor_hash, z$run == Sys.getenv("GITHUB_RUN_ID"))
        return(TRUE)
    }
    began <- as.numeric(Sys.time())
    setting <- rep("pooled", length(x))
    d <- v2_data(x, setting, a)
    meta <- data.frame(schema_version = 1L, dataset_id = key, benchmark = "RIbench", stage = "50000",
        Index = s$Index, case = NA_integer_, replicate = NA_integer_, Analyte = s$Analyte, N = s$N,
        n_routine = length(x), n_anchor = length(a), distribution = s$Distribution,
        pathological_fraction = s$fractionPathol, scenario = s$Scenario, rounding_step = 10^(-s$decimals),
        routine_seed = s$startSeed, anchor_seed = 17000000L + s$Index, routine_hash = input$routine_hash,
        anchor_hash = input$anchor_hash, ribench_hashcode = s$HashCode)

    start <- proc.time()
    first <- fit_attempt(r3_baseline(d, s$Index, s$decimals, 3L))
    extended <- if (is.null(first$error)) {
        fit_attempt(r3_baseline(d, s$Index, s$decimals, 5L, first$value))
    } else first
    fit_time <- proc.time() - start
    models <- if (is.null(extended$error)) extended$value else first$value
    if (length(models)) stopifnot(max(vapply(models, function(m) as.integer(m$k), integer(1))) ==
                                     if (is.null(extended$error)) 5L else 3L)
    row <- function(attempt, cap, label) {
        if (is.null(attempt$error)) return(comparator_row(attempt$value, meta, x, a, probs, cap, label))
        r <- template
        r[names(meta)] <- meta
        r$method <- label
        r$fitter <- "multistart-v3"
        r$fitter_source_commit <- "ec3dd35197759424d32a6a341320173e5db08c89"
        r$max_components <- cap
        r$lower_probability <- if (length(probs) == 2L) probs[1] else NA_real_
        r$upper_probability <- tail(probs, 1)
        r$routine_min <- min(x)
        r$technical <- TRUE
        r$declined <- r$lower_flag <- FALSE
        r$report_status <- "technical failure"
        r$reason <- attempt$error
        r
    }
    pair <- rbind(row(first, 3L, "anchoredRI"), row(extended, 5L, "armB"))
    # Reproduce the original one-group estimator on this runner and verify its study adapter.
    clock_ref <- proc.time()
    old <- frozen$fit_anchored_ri(x, setting, a, transform = "auto")
    now <- adapter$fit_anchored_ri(x, setting, a, transform = "auto")
    for (field in c("limits", "shift")) {
        stopifnot(isTRUE(all.equal(old[[field]], now[[field]], tolerance = 1e-10, scale = 1)),
                  all(abs(old[[field]] - now[[field]]) <= 1e-10))
    }
    stopifnot(identical(old$transform, now$transform), identical(old$components, now$components),
              identical(is.finite(old$selection$bic_original), is.finite(now$selection$bic_original)))
    b <- old$selection$bic_original - now$selection$bic_original
    stopifnot(all(abs(b[is.finite(b)]) <= 1e-10))
    frozen_guard <- bc_safeguard(now, x, a, probs)
    ref_time <- proc.time() - clock_ref
    historic <- list(transform = old$transform, components = old$components, shift = old$shift,
                     at_bound = old$at_bound, guard = frozen_guard)
    completed <- as.numeric(Sys.time())
    z <- list(key = key, manifest = manifest, routine = input$routine_hash, anchor = input$anchor_hash,
        run = Sys.getenv("GITHUB_RUN_ID"), attempt = Sys.getenv("GITHUB_RUN_ATTEMPT"),
        shard = Sys.getenv("BATCH"), image = Sys.getenv("ImageVersion"), models = models, pair = pair,
        frozen = historic, settings = "pooled", unique_values = length(unique(c(x, a))),
        began = began, ended = completed, fit_seconds = unname(fit_time["elapsed"]),
        cpu = unname(fit_time["user.self"] + fit_time["sys.self"]), frozen_seconds = unname(ref_time["elapsed"]))
    r3_save(z, path)
    TRUE
}
r3_pool(seq_along(inputs), workers, run_one)
