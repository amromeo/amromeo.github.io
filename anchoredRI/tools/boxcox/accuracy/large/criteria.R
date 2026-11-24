before <- subset(summ, group == "all" & filter == "before")
threshold <- before$capped_mean[before$method == "refineR"]
stopifnot(length(threshold) == 1L, is.finite(threshold))
criteria <- lapply(c("anchoredRI", "armB"), function(m) {
    raw <- scores[[m]]
    attempt <- before[before$method == m, ]
    e <- raw$zzDevAbs_Ov
    failures <- !is.finite(e) | raw$technical
    implausible <- is.finite(e) & e > 5
    burden <- sum(failures | implausible)
    score_pass <- bc_error_criterion(attempt, threshold)
    stopifnot(nrow(raw) == 1440L, !anyDuplicated(raw$Index))
    data.frame(method = m, attempted = nrow(raw), native_capped = attempt$capped_mean,
               native_n = attempt$capped_n, refineR_threshold = threshold, missing = attempt$missing,
               all_attempt_capped = attempt$all_attempt_capped, technical = sum(raw$technical),
               implausible = sum(implausible), implausible_or_failed = burden, max_implausible = 14L,
               score_pass, count_pass = burden <= 14L, passed = score_pass && burden <= 14L)
})
save_tab(do.call(rbind, criteria), "large-success-criteria")
pieces <- list()
for (m in c("anchoredRI", "armB")) {
    raw <- scores[[m]]
    groups_k <- list("1" = 1L, "2" = 2L, "3" = 3L, "4" = 4L, "5" = 5L,
                     "4 or 5" = 4:5, "1 to 3" = 1:3, "unselected" = NA_integer_)
    for (bucket in names(groups_k)) {
        q <- raw[raw$components %in% groups_k[[bucket]], ]
        k <- nrow(q)
        for (filter in c("before", "after")) {
            retained <- if (filter == "after") q[q$reported, ] else q
            e <- retained$zzDevAbs_Ov
            cs <- bc_cutoff_scores(e, retained$zzDevAbsCutoff_Ov)
            pieces[[length(pieces) + 1L]] <- data.frame(method = m, components = bucket, filter,
                attempts = nrow(raw), selected = k, share = k/nrow(raw), reported = sum(q$reported),
                declined = sum(q$declined), decline_rate = if (k) mean(q$declined) else NA_real_,
                technical = sum(q$technical), scored = sum(is.finite(e)), native_n = cs$native_n,
                native_capped = cs$native, mean_abs_z = mean(e[is.finite(e)]),
                implausible = sum(e > 5, na.rm = TRUE))
        }
    }
}
save_tab(do.call(rbind, pieces), "component-counts-and-errors")
