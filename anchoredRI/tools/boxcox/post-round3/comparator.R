# Export local improved-baseline estimates and the original Task E safeguard.
comparator_row <- function(models, meta, value, anchor, probs) {
  eligible <- Filter(function(m) m$k<=3L && m$transform %in% r3_families,models)
  tab <- v2_table(eligible,probs)
  selected <- tab[tab$selected,]
  stopifnot(nrow(selected)==1L)
  m <- eligible[[paste(selected$transform,selected$k)]]
  # The reliability fitter uses the Box-Cox identity x-1. The safeguard takes x.
  parameter <- m$p
  if (m$transform=="identity") parameter$mu <- parameter$mu+1
  candidate <- vector("list",m$k)
  candidate[[m$k]] <- parameter
  guard_fit <- list(transform=m$transform,components=m$k,shift=parameter$shift,candidate_parameters=candidate)
  guard <- bc_safeguard(guard_fit,value,anchor,probs)
  both <- length(probs)==2L
  raw <- if (both) c(selected$lower,selected$upper) else selected$upper
  stopifnot(isTRUE(all.equal(unname(guard$limits),unname(raw),tolerance=1e-12)))
  cbind(meta,data.frame(method="anchoredRI",fitter="multistart-v3",max_components=3L,
        fitter_source_commit="ec3dd35197759424d32a6a341320173e5db08c89",
        transform=m$transform,lambda=m$p$lambda,shift=m$p$shift,components=m$k,
        loglik=m$ll,bic=m$bic,convergence_code=m$code,at_bound=m$at_bound,
        lower_probability=if (both) probs[1] else NA_real_,upper_probability=tail(probs,1),
        lower=selected$lower,upper=selected$upper,routine_min=min(value),lrt_statistic=guard$lrt_statistic,
        lrt_p=guard$lrt_p,lower_flag=guard$lower_flag,declined=guard$declined,technical=guard$technical,
        report_status=if (guard$technical) "technical failure" else if (guard$declined) "declined" else "reported",
        reason=guard$reason,reported_lower=if (both) guard$reported[1] else NA_real_,
        reported_upper=tail(guard$reported,1),stringsAsFactors=FALSE))
}

write_comparator <- function(rows,out) {
  tab <- do.call(rbind,rows)
  stopifnot(!anyDuplicated(tab$dataset_id),!anyNA(tab$dataset_id))
  # Round-trip all finite binary64 values, including endpoints and very small LRT p-values.
  for (key in names(tab)[vapply(tab,is.double,logical(1))]) {
    x <- tab[[key]]
    tab[[key]] <- ifelse(is.na(x),NA_character_,sprintf("%.17g",x))
  }
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  write.csv(tab,out,row.names=FALSE,na="")
}
