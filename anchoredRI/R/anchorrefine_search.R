# SPDX-License-Identifier: GPL-3.0-only
# Ported from anchorRefine, commit 55b9d09ba7e17b13e341e67b0aa03855351bc25b.
# Source: R/stage2-search.R, R/stage2-selected.R and the peak-only branch of R/stage2-grid.R. See inst/NOTICE for attribution.

# Stage 2 is derived from refineR 2.0.0, GPL (>= 3), by Tatjana Ammer,
# Matthias Beck, Christopher M Rank and Andre Schuetzenmeister.
# The unchanged routine histogram and its cost remain in the upstream package.

# These internal functions belong to refineR 2.0.0. The package pins that
# version because the Stage 2 peak choice uses its search and histogram code.
.refineR_function <- function(name) utils::getFromNamespace(name, "refineR")

# Use the same two future seed calls as refineR::findRI(NBootstrap=0).
stage2_job <- function(seed, work) {
    seeds <- future.apply::future_lapply(1L, function(i) .Random.seed,
                  future.chunk.size = Inf, future.seed = seed)
    future.apply::future_lapply(1L, function(i) work(), future.seed = seeds,
                               future.packages = c("ash", "refineR"))[[1L]]
}

stage2_routine <- function(x) {
    stopifnot(is.numeric(x))
    x <- x[!is.na(x)]
    x <- x[x >= 0]
    stopifnot(length(x) > 10L)
    x
}

# Upstream defineSearchRegions jitters rounded data in the working copy.
# length(x) after duplication also sets the ASH smoother downstream.
stage2_working_data <- function(x, roundingBase) {
    if (length(x) > 100000) {
        if (!is.na(roundingBase)) {
            x <- x + stats::runif(length(x), -0.5 * roundingBase, 0.5 * roundingBase - 1e-20)
        }
        qData <- as.numeric(stats::quantile(x, seq(0, 1, length.out = 50000), type = 7))
    } else if (!is.na(roundingBase)) {
        x <- rep(x, times = floor(100000 / length(x)))
        qData <- x + stats::runif(length(x), -0.5 * roundingBase, 0.5 * roundingBase - 1e-20)
    } else {
        qData <- as.numeric(stats::quantile(x, seq(0, 1, length.out = 50000), type = 7))
    }
    list(qData = qData, n = length(x))
}

# Bounds of the *first* coarse ASH in upstream estimateAB, before trimming.
stage2_first_window <- function(x) {
    q <- as.numeric(stats::quantile(x, c(0, 0.01, 0.95, 0.99, 0.999), na.rm = TRUE))
    opt <- q[c(2L, 4L)]
    span <- diff(opt)
    opt[1L] <- max(opt[1L] - 0.1 * span, q[1L], 1e-20)
    opt[2L] <- min(opt[2L] + 0.1 * span, q[5L])
    det <- q[c(2L, 4L)]
    det[1L] <- max(det[1L], 1e-20)
    det[1L] <- det[1L] - 0.01 * diff(det)
    if (!all(is.finite(c(opt, det))) || opt[2L] <= opt[1L] || opt[2L] <= det[1L]) {
        stop("election: empty upstream first window", call. = FALSE)
    }
    ab <- c(max(det[1L], 1e-20), min(opt[2L], det[2L]))
    ab[1L] <- ab[1L] - 0.01 * diff(ab)
    ab
}

stage2_prepare <- function(x, anchors, h, seed) {
    restore <- preserve_rng()
    on.exit(restore())
    tryCatch(stage2_job(seed, function() {
        raw <- stage2_routine(x)
        data <- rep(raw, times = ceiling(250 / length(raw))) - 0
        work <- stage2_working_data(data, .refineR_function("findRoundingBase")(data))
        stage2_first_election(work$qData, stage2_first_window(work$qData), 24, anchors, h, 24L)
    }), error = function(e) stage2_empty_election(paste("prepare", conditionMessage(e)), 24L))
}

stage2_trace <- function() {
    data.frame(phase = character(), pass = integer(), lambda = double(),
               coarse = integer(), fine = integer(), stringsAsFactors = FALSE)
}

stage2_track <- function(x, ab, mStart, ctx, lambda = NULL, initial = FALSE) {
    peak <- stage2_guided_peak(x, ab, mStart, ctx$election, lambda = lambda, initial = initial)
    ctx$trace <- rbind(ctx$trace, data.frame(phase = if (is.null(lambda)) "raw" else "transformed",
        pass = ctx$pass, lambda = if (is.null(lambda)) NA_real_ else lambda,
        coarse = as.integer(peak$coarse), fine = as.integer(peak$fine), stringsAsFactors = FALSE))
    peak
}

# Only the peak picker is replaced in this private copy. estimateAB's own
# iterative expansion and clipping of both raw intervals remain unchanged.
stage2_estimate_ab <- function(x, ctx) {
    select_window <- .refineR_function("estimateAB")
    scope <- new.env(parent = environment(select_window))
    scope$findMainPeak <- function(x, ab, mStart, withHeight = FALSE, prevPeak = NULL) {
        initial <- !ctx$first_peak_seen
        ctx$first_peak_seen <- TRUE
        stage2_track(x, ab, mStart, ctx, initial = initial)
    }
    environment(select_window) <- scope
    select_window(x)
}

# The lambda is passed explicitly to the peak picker. The width equations
# below are those of Tatjana Ammer in refineR::defineSearchRegions.
stage2_regions <- function(x, lambdaVec, roundingBase, abEst, ctx) {
    work <- stage2_working_data(x, roundingBase)
    qData <- work$qData
    if (is.null(abEst)) abEst <- stage2_estimate_ab(qData, ctx)
    qData <- qData[qData >= abEst$abSearchReg[1L] & qData <= abEst$abSearchReg[2L]]
    muEst <- sigmaEst <- matrix(nrow = length(lambdaVec), ncol = 3L)
    densM <- min(48, max(24, ceiling(50 - (1/10) * work$n / 1000)))
    for (i in seq_along(lambdaVec)) {
        lambda <- lambdaVec[i]
        transformed <- suppressWarnings(refineR::BoxCox(qData, lambda))
        transformed <- transformed[is.finite(transformed)]
        abTr <- range(transformed)
        span <- diff(abTr)
        abTr <- c(abTr[1L] - 0.1 * span, abTr[2L] + 0.1 * span)
        main <- stage2_track(transformed, abTr, densM, ctx, lambda)
        dens <- main$Dens
        mode <- main$modEst
        steps <- seq(0.5, 0.95, by = 0.05)
        levels <- steps * dens$y[main$peakInd]
        fractions <- widthsL <- widthsR <- NULL
        for (j in seq_along(levels)) {
            left <- which(dens$x < mode & dens$y <= levels[j])
            right <- which(dens$x > mode & dens$y <= levels[j])
            if (length(left) && length(right)) {
                widthsL <- c(widthsL, mode - max(dens$x[left]))
                widthsR <- c(widthsR, min(dens$x[right]) - mode)
                fractions <- c(fractions, steps[j])
            }
        }
        if (!length(fractions))
            stop("regions: no usable routine peak widths", call. = FALSE)
        centers <- apply(cbind(mode - widthsL, mode + widthsR), 1L, mean)
        muRange <- range(centers, mode)
        muEst[i, ] <- c(lambda, muRange[1L] - 0.5 * diff(muRange),
                                  muRange[2L] + 0.5 * diff(muRange))
        scales <- pmin(widthsL, widthsR) / sqrt((-2) * log(fractions))
        sigRange <- range(scales)
        sigmaEst[i, ] <- c(lambda, max(0, sigRange[1L] - 2.5 * diff(sigRange)),
                                     sigRange[2L] + diff(sigRange))
    }
    list(abEst = abEst, muEst = muEst, sigmaEst = sigmaEst)
}

stage2_valid_fit <- function(fit) {
    if (!inherits(fit, "RWDRI")) return(FALSE)
    parameters <- unlist(fit[c("Mu", "Sigma", "Lambda", "Shift")], use.names = FALSE)
    length(parameters) == 4L && all(is.finite(parameters)) && fit$Sigma > 0 &&
        (is.null(fit$Status) || identical(as.character(fit$Status), "ResultFound"))
}

stage2_guided <- function(x, seed, election) {
    stage2_job(seed, function() {
        ctx <- new.env(parent = emptyenv())
        ctx$election <- election
        ctx$first_peak_seen <- FALSE
        ctx$trace <- stage2_trace()
        step <- "initialization"
        tryCatch({
            raw <- stage2_routine(x)
            data <- rep(raw, times = ceiling(250 / length(raw))) - 0
            rb <- .refineR_function("findRoundingBase")(data)
            lambda <- (1/12 * (0:12))^1.8170595
            best <- c(NA, NA, NA, NA, 1e9, NA)
            ab <- NULL
            for (pass in 1:2) {
                ctx$pass <- pass
                step <- paste0("pass", pass, "_regions")
                regions <- stage2_regions(data, lambda, rb, ab, ctx)
                if (pass == 1L) {
                    step <- "histogram"
                    ab <- regions$abEst
                    hist <- .refineR_function("generateHistData")(data, ab$abHist, rb)
                }
                step <- paste0("pass", pass, "_grid")
                # The frozen peak-only arm delegates the entire grid to refineR.
                best <- .refineR_function("testParam")(lambda, best, data, hist, regions,
                                                        7, 0.01, 0.1)
                if (pass == 1L) {
                    if (!is.finite(best[1L])) stop("no winning first-pass lambda")
                    idx <- min(length(lambda)-2, max(3, which.min(abs(best[1L]-lambda))))
                    lambda <- unlist(lapply((idx-2):(idx+1), function(j) {
                        seq(lambda[j], lambda[j+1], length.out = 5)[2:4]
                    }), use.names = FALSE)
                }
            }
            fit <- list(Data = raw, Method = "refineR", PkgVersion = utils::packageVersion("refineR"),
                        Model = "BoxCox", NBootstrap = 0, Lambda = best[1L], Mu = best[2L],
                        Sigma = best[3L], P = best[4L], Shift = 0, Cost = best[5L],
                        abOr = hist$abOr, roundingBase = rb, rf = best[6L])
            class(fit) <- "RWDRI"
            if (!stage2_valid_fit(fit)) stop("invalid guided model")
            list(fit = fit, reason = "guided", trace = ctx$trace)
        }, error = function(e) list(fit = NULL, reason = paste(step, conditionMessage(e)),
                                    trace = ctx$trace))
    })
}

stage2_baseline <- function(x, seed) {
    error <- NULL
    fit <- tryCatch(refineR::findRI(Data = x, model = "BoxCox", NBootstrap = 0, seed = seed),
        error = function(e) { error <<- conditionMessage(e); NULL })
    status <- if (stage2_valid_fit(fit)) "ResultFound" else "fit_failed"
    list(fit = fit, status = status, baseline_error = error)
}

# Stage 2's selected K = 24 peak-only arm. An unclear election or failed
# guided search runs the original zero-bootstrap refineR Box-Cox fit.
stage2_selected_fit <- function(x, anchors, h, seed) {
    election <- stage2_prepare(x, anchors, h, seed)
    restore <- preserve_rng()
    trace <- stage2_trace()
    reason <- election$reason
    if (isTRUE(election$ok)) {
        attempt <- tryCatch(stage2_guided(x, seed, election),
            error = function(e) list(fit = NULL, reason = paste("guided", conditionMessage(e)), trace = trace))
        trace <- attempt$trace
        reason <- attempt$reason
        if (stage2_valid_fit(attempt$fit)) {
            return(list(fit = attempt$fit, status = "ResultFound", fallback = FALSE, reason = reason,
                        baseline_status = NA_character_, election = election, peak_trace = trace))
        }
    }
    restore()
    baseline <- stage2_baseline(x, seed)
    baseline$fallback <- TRUE
    baseline$baseline_status <- baseline$status
    baseline$reason <- reason
    baseline$election <- election
    baseline$peak_trace <- trace
    baseline
}
