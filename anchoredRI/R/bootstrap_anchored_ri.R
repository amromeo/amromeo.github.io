# SPDX-License-Identifier: GPL-3.0-only
bootstrap_anchored_ri <- function(value, setting, anchor,
                                  transform = c("auto", "identity", "log", "shifted_log"),
                                  n_boot = 200L, conf = .95, seed = 123L,
                                  workers = 1L, shift_bounds = NULL) {
  transform <- match.arg(transform)
  integer_scalar <- function(x, minimum) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x==floor(x) && x>=minimum && x<=.Machine$integer.max
  if (!integer_scalar(n_boot,2L) || !integer_scalar(workers,1L) ||
      !integer_scalar(seed,0L) || !is.numeric(conf) || length(conf)!=1L ||
      !is.finite(conf) || conf<=0 || conf>=1)
    stop("Use n_boot >= 2, workers >= 1, a nonnegative integer seed, and 0 < conf < 1.")
  point <- fit_anchored_ri(value,setting,anchor,transform,shift_bounds)
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
    rows <- unlist(lapply(groups,function(g)sample(g,length(g),replace=TRUE)),use.names=FALSE)
    a <- sample(anchor,length(anchor),replace=TRUE)
    f <- tryCatch(fit_anchored_ri(value[rows],setting[rows],a,transform,shift_bounds),error=identity)
    ok <- !inherits(f,"error") && length(f$limits)==2L && all(is.finite(f$limits)) && f$limits[1]<f$limits[2]
    data.frame(replicate=b,lower=if(ok)unname(f$limits[1])else NA_real_,
      upper=if(ok)unname(f$limits[2])else NA_real_,
      transform=if(ok)f$transform else NA_character_,shift=if(ok)f$shift else NA_real_,
      components=if(ok)f$components else NA_integer_,at_bound=if(ok)f$at_bound else NA,
      status=if(ok)"ok"else"failed",
      error=if(inherits(f,"error"))conditionMessage(f)else if(!ok)"Non-finite or unordered limits"else"")
  }
  if (.Platform$OS.type=="windows" && workers>1L) {
    warning("Parallel bootstrap is available on Unix; using one worker on Windows.")
    workers <- 1L
  }
  draws <- if(workers==1L)lapply(seq_len(n_boot),one_draw)else
    parallel::mclapply(seq_len(n_boot),one_draw,mc.cores=workers,mc.set.seed=FALSE)
  if (!all(vapply(draws,is.data.frame,logical(1))))
    stop("A bootstrap worker failed outside the fitter; retry with fewer workers.")
  draws <- do.call(rbind,draws)
  ok <- draws$status=="ok"
  ci <- matrix(NA_real_,2,2,dimnames=list(c("lower","upper"),c("low","high")))
  if (sum(ok)>=max(2L,ceiling(.8*n_boot)))
    ci[,] <- t(vapply(draws[ok,c("lower","upper"),drop=FALSE],
       quantile,numeric(2),probs=c((1-conf)/2,1-(1-conf)/2),names=FALSE))
  else warning("Fewer than 80% of bootstrap fits succeeded; confidence limits are unavailable.")
  list(fit=point,limits=point$limits,ci=ci,confidence=conf,bootstrap=draws,
       requested=n_boot,successful=sum(ok),failed=sum(!ok),
       boundary_fits=sum(draws$at_bound[ok],na.rm=TRUE),seed=seed)
}
