# RIbench sets its cutoff metric to NA when a valid dataset has mean absolute z > 5.
bc_cutoff_scores <- function(error,native_cutoff) {
  stopifnot(length(error)==length(native_cutoff))
  observed <- is.finite(error)
  list(native=mean(native_cutoff,na.rm=TRUE),native_n=sum(is.finite(native_cutoff)),
       all_attempt=mean(ifelse(observed,pmin(error,5),5)))
}

bc_error_criterion <- function(arm,threshold) {
  is.finite(arm$capped_mean) && arm$capped_mean<=threshold &&
    (arm$missing==0L || arm$all_attempt_capped<=threshold)
}
