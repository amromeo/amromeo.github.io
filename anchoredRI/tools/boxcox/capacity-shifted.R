# Study-local copy of C2 shifted fitter. Only K extension and diagnostic returns differ.
# SPDX-License-Identifier: GPL-3.0-only
# A common, estimated location shift followed by log-normal mixture components.
# The log-Jacobian is included for every observation, including each anchor.
.anchored_shift_objective <- function(x, anchor, counts, k) {
  G <- ncol(counts)
  n_group <- colSums(counts)
  total_counts <- rowSums(counts)
  n <- sum(counts) + length(anchor)
  last_par <- last_result <- NULL
  function(par, gradient = FALSE) {
    if (!identical(par, last_par)) {
      mu <- par[seq_len(k)]
      sigma <- exp(par[k + seq_len(k)])
      shift <- tail(par, 1L)
      eta <- if (k == 1L) matrix(0, 1, G) else
        rbind(matrix(par[2L*k + seq_len(G*(k-1L))], ncol = G), rep(0, G))
      w <- exp(sweep(eta, 2, apply(eta, 2, max), "-"))
      w <- sweep(w, 2, colSums(w), "/")
      y <- log(x - shift)
      a <- log(anchor - shift)
      density <- vapply(seq_len(k), function(j) dnorm(y, mu[j], sigma[j], log = TRUE),
                         numeric(length(x)))
      ll <- sum(dnorm(a, mu[1], sigma[1], log = TRUE)) - sum(a) - sum(total_counts*y)
      contributions <- matrix(0, length(x), k)
      weight_score <- matrix(0, k - 1L, G)
      for (h in seq_len(G)) {
        parts <- sweep(density, 2, log(w[, h]), "+")
        peak <- do.call(pmax, as.data.frame(parts))
        scaled <- exp(parts - peak)
        denom <- rowSums(scaled)
        ll <- ll + sum(counts[, h] * (peak + log(denom)))
        mass <- scaled * (counts[, h] / denom)
        contributions <- contributions + mass
        if (k > 1L) weight_score[, h] <- (colSums(mass) - n_group[h]*w[, h])[-k]
      }
      u <- sweep(outer(y, mu, "-"), 2, sigma, "/")
      mean_score <- colSums(contributions*u) / sigma
      scale_score <- colSums(contributions*(u^2 - 1))
      ua <- (a-mu[1])/sigma[1]
      mean_score[1] <- mean_score[1] + sum(ua)/sigma[1]
      scale_score[1] <- scale_score[1] + sum(ua^2 - 1)
      shift_score <- sum((total_counts + rowSums(sweep(contributions*u, 2, sigma, "/")))/(x-shift)) +
        sum((1 + ua/sigma[1])/(anchor-shift))
      last_par <<- par + 0
      last_result <<- list(value = -ll/n,
        gradient = -c(mean_score, scale_score, weight_score, shift_score)/n,
        parameters = list(mu = mu, sigma = sigma, weights = w, shift = shift))
    }
    if (identical(gradient, "all")) return(last_result)
    if (gradient) last_result$gradient else last_result$value
  }
}

.fit_anchored_shifted <- function(value, setting, anchor, shift_bounds = NULL) {
  all_values <- c(value, anchor)
  location <- median(value)
  spread <- mad(value)
  if (spread <= 0) spread <- sd(value)
  if (!is.finite(spread) || spread <= 0 || length(unique(value)) < 4L)
    stop("Insufficient variation.")
  # Finite, data-derived search region prevents the support parameter from
  # collapsing onto the smallest result. Supplied bounds use original units.
  if (is.null(shift_bounds))
    shift_bounds <- min(all_values) + spread*c(-5, -.001)
  if (!is.numeric(shift_bounds) || length(shift_bounds) != 2L ||
      any(!is.finite(shift_bounds)) || shift_bounds[1] >= shift_bounds[2] ||
      shift_bounds[2] >= min(all_values))
    stop("shift_bounds must be two increasing finite numbers below every result.")
  bounds <- (shift_bounds-location)/spread
  v <- (value-location)/spread
  a <- (anchor-location)/spread
  x <- sort(unique(v))
  groups <- unique(as.character(setting))
  G <- length(groups)
  g <- match(as.character(setting), groups)
  counts <- matrix(tabulate(match(v,x)+(g-1L)*length(x), G*length(x)), ncol=G)
  n <- length(value)+length(anchor)
  zero <- min(bounds[2], max(bounds[1], -location/spread))
  shifts <- unique(c(zero, bounds[1]+.5*diff(bounds), bounds[1]+.9*diff(bounds)))
  fits <- lapply(seq_len(STUDY_K), function(k) {
    calc <- .anchored_shift_objective(x,a,counts,k)
    slots <- if (k==1L) list(numeric()) else if (k==2L) list(.15,.5,.85) else
      if(k==3L) list(c(.12,.85),c(.60,.90),c(.10,.40)) else
      list(seq(.12,.85,length.out=k-1L),seq(.60,.90,length.out=k-1L),seq(.10,.40,length.out=k-1L))
    lo <- c(rep(log(min(c(x,a))-bounds[2]),k), rep(log(.005),k),
            rep(-7,G*(k-1L)), bounds[1])
    hi <- c(rep(log(max(c(x,a))-bounds[1]),k), rep(log(5),k),
            rep(7,G*(k-1L)), bounds[2])
    starts <- unlist(lapply(shifts, function(d) lapply(slots, function(probs) {
      y <- log(v-d); aa <- log(a-d)
      c(mean(aa), unname(quantile(y,probs)), log(min(4.9,max(.006,sd(aa)))),
        rep(log(min(4.9,max(.006,sd(y)))),k-1L), rep(0,G*(k-1L)),d)
    })), recursive=FALSE)
    candidates <- lapply(starts,function(init) {
      tryCatch(optim(init,calc,function(p)calc(p,TRUE),method="L-BFGS-B",lower=lo,upper=hi,
            control=list(maxit=700L,factr=1e6)), error=function(e)NULL)
    })
    start_audit <- candidates
    candidates <- Filter(function(f) !is.null(f) && f$convergence==0L && is.finite(f$value),candidates)
    if (!length(candidates)) return(NULL)
    best <- candidates[[which.min(vapply(candidates,function(f)f$value,numeric(1)))]]
    p <- calc(best$par,"all")$parameters
    limits <- location + spread*(p$shift + exp(p$mu[1]+qnorm(c(.025,.975))*p$sigma[1]))
    if (any(!is.finite(limits)) || limits[1]>=limits[2]) return(NULL)
    list(starts=start_audit, p=p, limits=limits, bic=2*n*(best$value+log(spread))+length(best$par)*log(n),
      at_bound=any(abs(best$par-lo)<1e-5 | abs(best$par-hi)<1e-5),
      shift_at_bound=abs(p$shift-bounds[1])<1e-5 || abs(p$shift-bounds[2])<1e-5)
  })
  scores <- vapply(fits,function(f)if(is.null(f)) Inf else f$bic,numeric(1))
  candidate_limits <- t(vapply(fits,function(f)if(is.null(f))c(NA_real_,NA_real_)else f$limits,numeric(2)))
  dimnames(candidate_limits) <- list(as.character(seq_len(STUDY_K)),c("lower","upper"))
  shifts_out <- vapply(fits,function(f)if(is.null(f))NA_real_ else location+spread*f$p$shift,numeric(1))
  bound_hits <- vapply(fits,function(f)if(is.null(f))NA else f$at_bound,logical(1))
  if (!any(is.finite(scores))) return(list(limits=c(lower=NA_real_,upper=NA_real_),
    selected_bic=Inf,bic_original=scores,candidate_shift=shifts_out,candidate_at_bound=bound_hits))
  k <- which.min(scores); best <- fits[[k]]
  colnames(best$p$weights) <- groups
  list(candidate_parameters=lapply(fits,function(f) if(is.null(f)) NULL else list(mu=f$p$mu+log(spread),sigma=f$p$sigma,w=f$p$weights)),
       candidate_fits=fits, limits=setNames(best$limits,c("lower","upper")),components=k,weights=best$p$weights,
       bic=setNames(scores,seq_len(STUDY_K)),bic_original=setNames(scores,seq_len(STUDY_K)),selected_bic=scores[k],
       at_bound=best$at_bound,shift_at_bound=best$shift_at_bound,
       candidate_limits=candidate_limits,candidate_at_bound=bound_hits,
       candidate_shift=setNames(shifts_out,seq_len(STUDY_K)),transform="shifted_log",
       shift=shifts_out[k],shift_bounds=shift_bounds)
}
