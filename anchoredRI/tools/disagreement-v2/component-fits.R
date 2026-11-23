# SPDX-License-Identifier: GPL-3.0-only
# Copy of the frozen helpers. The two additional return fields expose the selected mixture.
.fit_anchored_fixed <- function(value, setting, anchor, transform = c("identity", "log")) {
  transform <- match.arg(transform)
  if (!is.numeric(value) || !is.numeric(anchor) || length(anchor) < 2L ||
      length(value) != length(setting) || length(unique(setting)) != 2L ||
      anyNA(setting) || any(!is.finite(c(value, anchor)))) {
    stop("Supply numeric results from two settings and independent reference values.")
  }
  if (transform == "log" && any(c(value, anchor) <= 0)) stop("Log requires positive values.")
  forward <- if (transform == "log") log else identity
  inverse <- if (transform == "log") exp else identity
  y <- forward(value)
  center <- median(y)
  scale <- mad(y)
  if (scale <= 0) scale <- sd(y)
  if (!is.finite(scale) || scale <= 0 || length(unique(y)) < 4L) stop("Insufficient variation.")
  z <- (y - center) / scale
  a <- (forward(anchor) - center) / scale
  x <- sort(unique(z))
  groups <- unique(as.character(setting))
  g <- match(as.character(setting), groups)
  counts <- matrix(tabulate(match(z, x) + (g - 1L) * length(x), 2L * length(x)),
                   ncol = 2L, dimnames = list(NULL, groups))
  n_group <- colSums(counts)
  n <- length(z) + length(a)
  # Component means, log SDs, and non-baseline log odds for each setting.
  fit_k <- function(k) {
    unpack <- function(par) {
      eta <- rbind(matrix(tail(par, 2L * (k - 1L)), ncol = 2L), c(0, 0))
      weights <- exp(sweep(eta, 2, apply(eta, 2, max), "-"))
      list(mu = par[seq_len(k)], sigma = exp(par[k + seq_len(k)]),
           w = sweep(weights, 2, colSums(weights), "/"))
    }
    calc <- function(par, gradient = FALSE) {
      p <- unpack(par)
      log_density <- vapply(seq_len(k), function(j) {
        dnorm(x, p$mu[j], p$sigma[j], log = TRUE)
      }, numeric(length(x)))
      ll <- sum(dnorm(a, p$mu[1], p$sigma[1], log = TRUE))
      contributions <- matrix(0, length(x), k)
      weight_score <- matrix(0, k - 1L, 2L)
      for (h in 1:2) {
        parts <- sweep(log_density, 2, log(p$w[, h]), "+")
        peak <- do.call(pmax, as.data.frame(parts))
        scaled <- exp(parts - peak)
        denom <- rowSums(scaled)
        ll <- ll + sum(counts[, h] * (peak + log(denom)))
        if (gradient) {
          mass <- scaled * (counts[, h] / denom)
          contributions <- contributions + mass
          weight_score[, h] <- (colSums(mass) - n_group[h] * p$w[, h])[-k]
        }
      }
      if (!gradient) return(-ll / n)
      u <- sweep(outer(x, p$mu, "-"), 2, p$sigma, "/")
      mean_score <- colSums(contributions * u) / p$sigma
      scale_score <- colSums(contributions * (u^2 - 1))
      ua <- (a - p$mu[1]) / p$sigma[1]
      mean_score[1] <- mean_score[1] + sum(ua) / p$sigma[1]
      scale_score[1] <- scale_score[1] + sum(ua^2 - 1)
      -c(mean_score, scale_score, weight_score) / n
    }
    slots <- if (k == 2L) list(.15, .5, .85) else list(c(.12,.85), c(.60,.90), c(.10,.40))
    ref_sd <- min(19.9, max(.06, sd(a)))
    lo <- c(rep(min(c(z,a)), k), rep(log(.05), k), rep(-7, 2L * (k-1L)))
    hi <- c(rep(max(c(z,a)), k), rep(log(20), k), rep(7, 2L * (k-1L)))
    fits <- lapply(slots, function(probs) {
      init <- c(mean(a), unname(quantile(z, probs)), log(ref_sd),
                rep(log(1), k-1L), rep(0, 2L * (k-1L)))
      tryCatch(optim(init, calc, function(p) calc(p, TRUE), method = "L-BFGS-B",
            lower = lo, upper = hi, control = list(maxit = 450L, factr = 1e6)),
            error = function(e) NULL)
    })
    fits <- Filter(function(f) !is.null(f) && f$convergence == 0L && is.finite(f$value), fits)
    if (!length(fits)) return(NULL)
    fit <- fits[[which.min(vapply(fits, function(f) f$value, numeric(1)))]]
    list(p = unpack(fit$par), bic = 2*n*fit$value + length(fit$par)*log(n),
         at_bound = any(abs(fit$par-lo)<1e-5 | abs(fit$par-hi)<1e-5))
  }
  all_z <- c(z, a)
  one_mu <- mean(all_z)
  one_sd <- sqrt(mean((all_z - one_mu)^2))
  one <- list(p = list(mu = one_mu, sigma = one_sd, w = matrix(1, 1L, 2L)),
              bic = -2 * sum(dnorm(all_z, one_mu, one_sd, log = TRUE)) + 2*log(n),
              at_bound = FALSE)
  fits <- list(one, fit_k(2L), fit_k(3L))
  bics <- vapply(fits, function(f) if (is.null(f)) Inf else f$bic, numeric(1))
  k <- which.min(bics)
  p <- fits[[k]]$p
  colnames(p$w) <- groups
  candidate_limits <- t(vapply(fits, function(f) {
    if (is.null(f)) return(c(NA_real_, NA_real_))
    inverse(center + scale * (f$p$mu[1] + qnorm(c(.025, .975)) * f$p$sigma[1]))
  }, numeric(2)))
  dimnames(candidate_limits) <- list(as.character(1:3), c("lower", "upper"))
  limits <- inverse(center + scale * (p$mu[1] + qnorm(c(.025,.975))*p$sigma[1]))
  # Density conversion back to the measurement scale. This constant cancels
  # across component counts, but is essential when comparing transformations.
  correction <- 2*n*log(scale) + if (transform == "log") 2*sum(log(c(value, anchor))) else 0
  list(limits = setNames(limits, c("lower","upper")), components = k,
       component_mu = p$mu, component_sd = p$sigma,
       weights = p$w, bic = setNames(bics, 1:3), at_bound = fits[[k]]$at_bound,
       bic_original = setNames(bics + correction, 1:3),
       selected_bic = bics[k] + correction, transform = transform, shift = 0,
       candidate_shift = setNames(rep(0, 3), 1:3), shift_bounds = NULL,
       candidate_limits = candidate_limits,
       candidate_at_bound = vapply(fits, function(f) {
         if (is.null(f)) NA else f$at_bound
       }, logical(1)))
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
  g <- match(as.character(setting), groups)
  counts <- matrix(tabulate(match(v,x)+(g-1L)*length(x), 2L*length(x)), ncol=2L)
  n <- length(value)+length(anchor)
  zero <- min(bounds[2], max(bounds[1], -location/spread))
  shifts <- unique(c(zero, bounds[1]+.5*diff(bounds), bounds[1]+.9*diff(bounds)))
  fits <- lapply(1:3, function(k) {
    calc <- .anchored_shift_objective(x,a,counts,k)
    slots <- if (k==1L) list(numeric()) else if (k==2L) list(.15,.5,.85) else
      list(c(.12,.85),c(.60,.90),c(.10,.40))
    lo <- c(rep(log(min(c(x,a))-bounds[2]),k), rep(log(.005),k),
            rep(-7,2L*(k-1L)), bounds[1])
    hi <- c(rep(log(max(c(x,a))-bounds[1]),k), rep(log(5),k),
            rep(7,2L*(k-1L)), bounds[2])
    starts <- unlist(lapply(shifts, function(d) lapply(slots, function(probs) {
      y <- log(v-d); aa <- log(a-d)
      c(mean(aa), unname(quantile(y,probs)), log(min(4.9,max(.006,sd(aa)))),
        rep(log(min(4.9,max(.006,sd(y)))),k-1L), rep(0,2L*(k-1L)),d)
    })), recursive=FALSE)
    candidates <- lapply(starts,function(init) {
      tryCatch(optim(init,calc,function(p)calc(p,TRUE),method="L-BFGS-B",lower=lo,upper=hi,
            control=list(maxit=700L,factr=1e6)), error=function(e)NULL)
    })
    candidates <- Filter(function(f) !is.null(f) && f$convergence==0L && is.finite(f$value),candidates)
    if (!length(candidates)) return(NULL)
    best <- candidates[[which.min(vapply(candidates,function(f)f$value,numeric(1)))]]
    p <- calc(best$par,"all")$parameters
    limits <- location + spread*(p$shift + exp(p$mu[1]+qnorm(c(.025,.975))*p$sigma[1]))
    if (any(!is.finite(limits)) || limits[1]>=limits[2]) return(NULL)
    list(p=p, limits=limits, bic=2*n*(best$value+log(spread))+length(best$par)*log(n),
      at_bound=any(abs(best$par-lo)<1e-5 | abs(best$par-hi)<1e-5),
      shift_at_bound=abs(p$shift-bounds[1])<1e-5 || abs(p$shift-bounds[2])<1e-5)
  })
  scores <- vapply(fits,function(f)if(is.null(f)) Inf else f$bic,numeric(1))
  candidate_limits <- t(vapply(fits,function(f)if(is.null(f))c(NA_real_,NA_real_)else f$limits,numeric(2)))
  dimnames(candidate_limits) <- list(as.character(1:3),c("lower","upper"))
  shifts_out <- vapply(fits,function(f)if(is.null(f))NA_real_ else location+spread*f$p$shift,numeric(1))
  bound_hits <- vapply(fits,function(f)if(is.null(f))NA else f$at_bound,logical(1))
  if (!any(is.finite(scores))) return(list(limits=c(lower=NA_real_,upper=NA_real_),
    selected_bic=Inf,bic_original=scores,candidate_shift=shifts_out,candidate_at_bound=bound_hits))
  k <- which.min(scores); best <- fits[[k]]
  colnames(best$p$weights) <- groups
  list(limits=setNames(best$limits,c("lower","upper")),components=k,weights=best$p$weights,
       component_mu=best$p$mu, component_sd=best$p$sigma,
       bic=setNames(scores,1:3),bic_original=setNames(scores,1:3),selected_bic=scores[k],
       at_bound=best$at_bound,shift_at_bound=best$shift_at_bound,
       candidate_limits=candidate_limits,candidate_at_bound=bound_hits,
       candidate_shift=setNames(shifts_out,1:3),transform="shifted_log",
       shift=shifts_out[k],shift_bounds=shift_bounds)
}
