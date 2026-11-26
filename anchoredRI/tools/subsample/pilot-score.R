# The pilot's provisional paired scores use the same grouped routine as Task E.
# The full benchmark must score the entire 1,440-row method cohort again.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L, as.character(utils::packageVersion("RIbench")) == "1.0.2")
endpoints <- read.csv(args[1], stringsAsFactors = FALSE)
design <- read.csv(args[2], stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(design$Index), !anyDuplicated(endpoints[c("Index", "method")]))
score_group <- function(rows) {
    truth <- design[match(rows$Index, design$Index), ]
    analyte <- unique(truth$Analyte)
    stopifnot(length(analyte) == 1L, !anyNA(truth$Index))
    estimates <- data.frame(Index = rows$Index, GT_LRL = truth$GT_LRL, GT_URL = truth$GT_URL,
                            LRL = rows$lower, URL = rows$upper)
    # Task E presents a zero lower bound to RIbench for one-sided CRP.
    if (analyte == "CRP") estimates$LRL <- 0
    method <- unique(rows$method)
    result <- RIbench:::computePerfMeasAll(
        analyte, method, list(estimates), truth, cutoffZ = 5
    )$resRIs[[1L]]
    data.frame(Index = rows$Index, method = method,
               z_lower = if (analyte == "CRP") NA_real_ else result$zzDev_LRL,
               z_upper = result$zzDev_URL, d = result$zzDevAbs_Ov,
               native_cutoff = result$zzDevAbsCutoff_Ov)
}
groups <- split(endpoints, list(endpoints$method, design$Analyte[match(endpoints$Index, design$Index)]),
                drop = TRUE)
options(digits = 17)
write.csv(do.call(rbind, lapply(groups, score_group)), args[3], row.names = FALSE, na = "")
