# SPDX-License-Identifier: GPL-3.0-only
.validate_ri_reporting <- function(safeguard, probs) {
  if (!is.logical(safeguard) || length(safeguard)!=1L || is.na(safeguard))
    stop("safeguard must be TRUE or FALSE.")
  if (!is.numeric(probs) || !length(probs) %in% 1:2 ||
      any(!is.finite(probs)) || any(probs<=0 | probs>=1) ||
      (length(probs)==2L && probs[1]>=probs[2]))
    stop("probs must contain one probability or two increasing probabilities between 0 and 1.")
}
.ri_limit_names <- function(probs) {
  if(length(probs)==2L)c("lower","upper") else if(probs<.5)"lower" else "upper"
}
.ri_forward <- function(x, transform, shift) {
  switch(transform, identity=x, log=log(x), shifted_log=log(x-shift),
         stop("Unknown fitted transform."))
}
.ri_inverse <- function(x, transform, shift) {
  if(transform=="identity") x else shift+exp(x)
}
.ri_reference_parameters <- function(fit) {
  z <- .ri_forward(fit$limits, fit$transform, fit$shift)
  if (length(z) != 2L || any(!is.finite(z)) || z[2] <= z[1])
    stop("Invalid fitted reference parameters.")
  list(mu = mean(z), sigma = unname(diff(z)) / (2 * stats::qnorm(.975)))
}
.ri_requested_limits <- function(fit, probs) {
  if (identical(as.numeric(probs), c(.025, .975))) return(unname(fit$limits))
  parameters <- .ri_reference_parameters(fit)
  .ri_inverse(parameters$mu + stats::qnorm(probs) * parameters$sigma,
              fit$transform, fit$shift)
}
# Use the original central-95% diagnostic limits, even for a one-sided request.
# This reproduces Step B and the saved RIbench diagnostic calculation.
.ri_checks <- function(fit, value, anchor, probs=c(.025,.975), enabled=TRUE) {
  if (!enabled) {
    estimates <- .ri_requested_limits(fit, probs)
    return(list(limits=stats::setNames(estimates,.ri_limit_names(probs)),
                lrt_statistic=NA_real_, lrt_p=NA_real_, lrt_flag=NA,
                lower_flag=NA, lower_check_skipped=TRUE, routine_min=NA_real_,
                technical=FALSE, reasons=character(), error=""))
  }
  statistic <- NA_real_
  estimates <- rep(NA_real_,length(probs))
  error <- ""
  tryCatch({
    estimates <- .ri_requested_limits(fit, probs)
    parameters <- .ri_reference_parameters(fit)
    y <- .ri_forward(anchor,fit$transform,fit$shift)
    own_sd <- sqrt(mean((y-mean(y))^2))
    if(any(!is.finite(y)) || !is.finite(own_sd) || own_sd<=0)
      stop("Anchor likelihood is unavailable.")
    statistic <- max(0,2*(sum(stats::dnorm(y,mean(y),own_sd,log=TRUE))-
                          sum(stats::dnorm(y,parameters$mu,parameters$sigma,log=TRUE))))
  },error=function(e)error<<-conditionMessage(e))
  estimates <- stats::setNames(estimates, .ri_limit_names(probs))
  .ri_flags(estimates, statistic, min(value), error)
}
.ri_flags <- function(estimates, statistic, routine_min, error="") {
  p <- stats::pchisq(statistic, df=2, lower.tail=FALSE)
  one_sided <- length(estimates) == 1L
  lower <- if(one_sided)FALSE else
    if(is.finite(estimates[1]))unname(estimates[1]<routine_min) else NA
  lrt <- if(is.finite(p))p<.01 else NA
  technical <- !is.finite(p) || any(!is.finite(estimates)) || is.na(lower) ||
    (length(estimates) == 2L && estimates[1] >= estimates[2])
  reasons <- c(if(isTRUE(lrt))"routine data disagree with the screened sample",
               if(isTRUE(lower))"estimated lower limit is below every observed result")
  list(limits=estimates,lrt_statistic=statistic,lrt_p=p,lrt_flag=lrt,
       lower_flag=lower,lower_check_skipped=one_sided,
       routine_min=routine_min,technical=technical,reasons=reasons,error=error)
}
.apply_ri_safeguard <- function(fit,value,anchor,safeguard=TRUE,probs=c(.025,.975)) {
  .validate_ri_reporting(safeguard,probs)
  check <- .ri_checks(fit,value,anchor,probs,enabled=safeguard)
  anchors <- tryCatch(.anchors_auto(anchor),error=function(e)
    list(limits=c(NA_real_,NA_real_),error=conditionMessage(e)))
  if(is.null(anchors$error))anchors$error <- ""
  if(length(probs)!=2L || !identical(as.numeric(probs),c(.025,.975))) {
    anchors$limits <- if(!is.null(anchors$mu))
      .ri_inverse(anchors$mu+stats::qnorm(probs)*anchors$sd,anchors$transform,anchors$shift) else
      rep(NA_real_,length(probs))
  }
  anchors$limits <- stats::setNames(anchors$limits,.ri_limit_names(probs))
  decision <- .ri_reporting_decision(check,safeguard,length(value))
  status <- decision$status
  fit$diagnostics <- list(anchored=list(limits=check$limits,central_limits=fit$limits,
      transform=fit$transform,shift=fit$shift,components=fit$components),
      anchors_only=anchors)
  fit$limits <- if(status=="reported")check$limits else
    stats::setNames(rep(NA_real_,length(probs)),.ri_limit_names(probs))
  fit$probs <- probs
  fit$status <- status
  fit$declined <- decision$declined
  fit$technical_withholding <- decision$technical_withholding
  fit$reason <- decision$reason
  fit$recommendation <- decision$recommendation
  fit$safeguard <- c(list(enabled=safeguard),check[setdiff(names(check),"limits")])
  fit
}

.ri_reporting_decision <- function(check,safeguard,n_routine) {
  technical <- safeguard && check$technical
  declined <- safeguard && !technical && length(check$reasons)>0L
  status <- if(technical)"withheld" else if(declined)"declined" else "reported"
  recommendation <- if(status!="reported")paste(
    "Obtain at least 120 screened patients or verify a suitable published interval with 20 screened patients.",
    if(n_routine>=50000L)"Consider an established indirect method such as refineR." else "") else ""
  list(status=status,declined=declined,technical_withholding=technical,
       reason=if(safeguard)c(check$reasons,if(technical)"safeguard could not be evaluated")else character(),
       recommendation=trimws(recommendation))
}
