# SPDX-License-Identifier: GPL-3.0-only
bootstrap_anchored_ri <- function(value, setting, anchor,
                                  transform = c("auto", "identity", "log", "shifted_log"),
                                  n_boot = 200L, conf = .95, seed = 123L,
                                  workers = 1L, shift_bounds = NULL,
                                  safeguard = TRUE, probs = c(.025,.975),
                                  rounding_unit = NULL, fit_seed = 1L) {
  transform <- match.arg(transform)
  integer_scalar <- function(x, minimum) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x==floor(x) && x>=minimum && x<=.Machine$integer.max
  if (!integer_scalar(n_boot,2L) || !integer_scalar(workers,1L) ||
      !integer_scalar(seed,0L) || !is.numeric(conf) || length(conf)!=1L ||
      !is.finite(conf) || conf<=0 || conf>=1)
    stop("Use n_boot >= 2, workers >= 1, a nonnegative integer seed, and 0 < conf < 1.")
  point <- fit_anchored_ri(value, setting, anchor, transform, shift_bounds, safeguard, probs,
                           rounding_unit = rounding_unit, fit_seed = fit_seed)
  existed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (existed) previous_seed <- get(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  on.exit({
    if (existed) assign(".Random.seed",previous_seed,envir=.GlobalEnv)
    else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max,n_boot)
  groups <- split(seq_along(value),as.character(setting))
  one_draw <- function(b) {
    set.seed(seeds[b])
    rows <- unlist(lapply(groups,function(g)g[sample.int(length(g),length(g),replace=TRUE)]),use.names=FALSE)
    a <- sample(anchor,length(anchor),replace=TRUE)
    f <- tryCatch(fit_anchored_ri(value[rows],setting[rows],a,transform,shift_bounds,
                                  safeguard, probs, rounding_unit = point$rounding_unit,
                                  fit_seed = fit_seed), error=identity)
    fit_error <- inherits(f,"error")
    invalid <- !fit_error && f$status=="reported" &&
      (any(!is.finite(f$limits)) || (length(probs)==2L && f$limits[1]>=f$limits[2]))
    failed <- fit_error || invalid
    take <- function(x,name)if(name %in% names(x))unname(x[name])else NA_real_
    row <- data.frame(replicate=b,
      lower=if(failed)NA_real_ else take(f$limits,"lower"),
      upper=if(failed)NA_real_ else take(f$limits,"upper"),
      diagnostic_lower=if(failed)NA_real_ else take(f$diagnostics$anchored$limits,"lower"),
      diagnostic_upper=if(failed)NA_real_ else take(f$diagnostics$anchored$limits,"upper"),
      anchors_lower=if(failed)NA_real_ else take(f$diagnostics$anchors_only$limits,"lower"),
      anchors_upper=if(failed)NA_real_ else take(f$diagnostics$anchors_only$limits,"upper"),
      transform=if(failed)NA_character_ else f$transform,
      shift=if(failed)NA_real_ else f$shift,
      components=if(failed)NA_integer_ else f$components,
      at_bound=if(failed)NA else f$at_bound,
      lrt_statistic=if(failed)NA_real_ else f$safeguard$lrt_statistic,
      lrt_p=if(failed)NA_real_ else f$safeguard$lrt_p,
      lrt_flag=if(failed)NA else f$safeguard$lrt_flag,
      lower_flag=if(failed)NA else f$safeguard$lower_flag,
      lower_check_skipped=!safeguard || length(probs)==1L,
      status=if(failed)"failed" else if(f$status=="reported")"ok" else f$status,
      reason=if(failed)"" else paste(f$reason,collapse="; "),
      error=if(fit_error)conditionMessage(f) else if(invalid)"Non-finite or unordered limits" else f$safeguard$error)
    list(row=row,diagnostics=if(failed)NULL else f$diagnostics)
  }
  if (.Platform$OS.type=="windows" && workers>1L) {
    warning("Parallel bootstrap is available on Unix; using one worker on Windows.")
    workers <- 1L
  }
  draws <- if(workers==1L)lapply(seq_len(n_boot),one_draw)else
    parallel::mclapply(seq_len(n_boot),one_draw,mc.cores=workers,mc.set.seed=FALSE)
  if (!all(vapply(draws,function(x)is.list(x)&&is.data.frame(x$row),logical(1))))
    stop("A bootstrap worker failed outside the fitter; retry with fewer workers.")
  diagnostics <- lapply(draws,function(x)x$diagnostics)
  draws <- do.call(rbind,lapply(draws,function(x)x$row))
  ok <- draws$status=="ok"
  limit_names <- .ri_limit_names(probs)
  ci <- matrix(NA_real_,length(probs),2,dimnames=list(limit_names,c("low","high")))
  if (point$status=="reported" && sum(ok)>=max(2L,ceiling(.8*n_boot)))
    ci[,] <- t(vapply(draws[ok,limit_names,drop=FALSE],
       quantile,numeric(2),probs=c((1-conf)/2,1-(1-conf)/2),names=FALSE))
  else if(point$status=="reported")
    warning("Fewer than 80% of bootstrap draws were accepted; confidence limits are unavailable.")
  list(fit=point,limits=point$limits,probs=point$probs,ci=ci,confidence=conf,bootstrap=draws,
       status=point$status,declined=point$declined,reason=point$reason,
       recommendation=point$recommendation,safeguard=point$safeguard,
       diagnostics=list(point=point$diagnostics,bootstrap=diagnostics),
       ci_conditional_on_safeguard=safeguard,
       requested=n_boot,successful=sum(ok),failed=sum(draws$status=="failed"),
       declined_draws=sum(draws$status=="declined"),
       withheld_draws=sum(draws$status=="withheld"),
       boundary_fits=sum(draws$at_bound[ok],na.rm=TRUE),seed=seed)
}
