# SPDX-License-Identifier: GPL-3.0-only
# Task E Arm B. Its identity working coordinate is x - 1.
.routine_scale <- function(x) {
  spread <- mad(x)
  if (spread > 0) spread else sd(x)
}

.measurement_unit <- function(x, unit) {
  on_grid <- function(step) {
    scaled <- x / step
    if (any(!is.finite(scaled))) return(FALSE)
    nearest <- round(scaled)
    on_step <- x == nearest * step
    # Decimal grid points such as 0.05 are integer multiples of 1/20.
    # At large indices division can recover the stored value more accurately.
    divisor <- round(1 / step)
    if (step < 1 && is.finite(divisor) && abs(step * divisor - 1) <= .Machine$double.eps)
      on_step <- on_step | x == nearest / divisor
    # A tiny nonzero observation must not pass as zero on a coarser inferred grid.
    tolerance <- if (is.null(unit)) pmin(1e-6, 8 * .Machine$double.eps * abs(scaled)) else 1e-6
    all(on_step | abs(scaled - nearest) < tolerance)
  }
  if (!is.null(unit)) {
    if (!is.numeric(unit) || length(unit) != 1L || !is.finite(unit) || unit < 0 ||
        (unit > 0 && !on_grid(unit)))
      stop("rounding_unit must be zero or a positive reporting increment shared by all results.")
    return(as.double(unit))
  }
  steps <- as.vector(outer(c(5, 2.5, 2, 1), 10^(3:-8)))
  for (step in steps) if (on_grid(step)) return(step)
  0
}

.mixture_data <- function(value, setting, anchor, shift_bounds) {
  points <- sort(unique(c(value, anchor)))
  groups <- unique(as.character(setting))
  ix <- match(value, points)
  group_index <- match(as.character(setting), groups)
  counts <- tabulate(ix + (group_index - 1L) * length(points), length(points) * length(groups))
  spread <- .routine_scale(value)
  if (!is.finite(spread) || spread <= 0 || length(unique(value)) < 4L)
    stop("Insufficient variation.")
  bounds <- if (is.null(shift_bounds)) min(points) + spread * c(-5, -.001) else shift_bounds
  if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds)) ||
      bounds[1] >= bounds[2] || bounds[2] >= min(points))
    stop("shift_bounds must be two increasing finite numbers below every result.")
  list(value = value, anchor = anchor, x = as.double(points), groups = groups,
       counts = matrix(as.double(counts), ncol = length(groups)),
       anchors = as.double(tabulate(match(anchor, points), length(points))),
       g = length(groups), n = length(value) + length(anchor),
       location = median(value), spread = spread, bounds = bounds,
       shifts = unique(c(max(bounds[1], min(bounds[2], 0)), mean(bounds), bounds[1] + .9 * diff(bounds))),
       zero = bounds[1] <= 0 && bounds[2] >= 0, positive = points[1] > 0)
}

.mixture_coordinates <- function(d, family) {
  switch(family,
         identity = c(median(d$value) - 1, .routine_scale(d$value)),
         log = c(median(log(d$value)), .routine_scale(log(d$value))),
         shifted_log = c(log(d$spread), 1))
}

.mixture_bounds <- function(d, family, k, coords) {
  values <- range(d$x)
  bounds <- switch(family, identity = values - 1, log = log(values),
                   shifted_log = log(values - rev(d$bounds)))
  shifted <- family == "shifted_log"
  sd_range <- if (shifted) c(.005, 5) else c(.05, 20)
  b <- (bounds - coords[1]) / coords[2]
  offsets <- (d$bounds - d$location) / d$spread
  n_odds <- d$g * (k - 1L)
  list(lo = c(rep(b[1], k), rep(log(sd_range[1]), k), rep(-7, n_odds), if (shifted) offsets[1]),
       hi = c(rep(b[2], k), rep(log(sd_range[2]), k), rep(7, n_odds), if (shifted) offsets[2]))
}

.pack_mixture <- function(d, family, p, coords) {
  k <- length(p$mu)
  odds <- if (k > 1L) as.vector(log(sweep(p$w[-k, , drop = FALSE], 2, p$w[k, ], "/")))
  c((p$mu - coords[1]) / coords[2], log(p$sigma / coords[2]), odds,
    if (family == "shifted_log") (p$shift - d$location) / d$spread)
}

.unpack_mixture <- function(d, family, k, v, coords) {
  odds <- if (k == 1L) matrix(0, 1, d$g) else
    rbind(matrix(v[2*k + seq_len(d$g * (k-1L))], ncol = d$g), 0)
  w <- exp(odds)
  list(mu = coords[1] + coords[2] * v[seq_len(k)],
       sigma = coords[2] * exp(v[k + seq_len(k)]), w = sweep(w, 2, colSums(w), "/"),
       shift = if (family == "shifted_log") d$location + d$spread * tail(v, 1L) else 0)
}

.grid_starts <- function(d, family, k, coords) {
  shifted <- family == "shifted_log"
  offsets <- if (shifted) d$shifts else 0
  slots <- if (k == 1L) list(numeric()) else if (k == 2L) list(.15, .5, .85) else
    list(seq(.12, .85, length.out = k-1L), seq(.60, .90, length.out = k-1L),
         seq(.10, .40, length.out = k-1L))
  starts <- list()
  for (shift in offsets) {
    y <- if (family == "identity") d$value - 1 else log(d$value - shift)
    a <- if (family == "identity") d$anchor - 1 else log(d$anchor - shift)
    for (q in slots) {
      analytic <- k == 1L && !shifted
      if (analytic) {
        both <- c(y, a)
        mu <- mean(both)
        sigma <- sqrt(mean((both - mu)^2))
      } else {
        mu <- c(mean(a), unname(quantile(y, q)))
        sigma <- c(sd(a), rep(if (shifted) sd(y) else .routine_scale(y), k-1L))
        if (shifted) sigma <- coords[2] * pmin(4.9, pmax(.006, sigma / coords[2]))
        else sigma[1] <- coords[2] * min(19.9, max(.06, sigma[1] / coords[2]))
      }
      p <- list(mu = mu, sigma = sigma, w = matrix(1/k, k, d$g), shift = shift)
      starts[[length(starts)+1L]] <- list(p = p, analytic = analytic)
    }
  }
  starts
}

.split_component <- function(p, j, apart = FALSE) {
  k <- length(p$mu)
  p$mu <- c(p$mu, p$mu[j])
  p$sigma <- c(p$sigma, p$sigma[j])
  p$w <- rbind(p$w, p$w[j, , drop = FALSE] / 2)
  p$w[j, ] <- p$w[j, ] / 2
  if (apart) {
    if (j != 1L) p$mu[j] <- p$mu[j] - .3 * p$sigma[j]
    p$mu[k+1L] <- p$mu[k+1L] + .3 * p$sigma[j]
  }
  p
}

.nested_starts <- function(previous) {
  p <- previous$p
  widest <- which.max(p$sigma)
  k <- length(p$mu) + 1L
  split <- list(p = .split_component(p, widest, TRUE))
  # Relabel only nonreference components to keep log odds inside their bounds.
  for (j in unique(c(widest, seq_along(p$mu)))) {
    duplicate <- .split_component(p, j)
    for (last in rev(2:k)) {
      order <- c(setdiff(seq_len(k), last), last)
      weights <- duplicate$w[order, , drop = FALSE]
      if (any(abs(log(sweep(weights, 2, weights[k, ], "/"))) > 7 + 1e-11)) next
      duplicate$mu <- duplicate$mu[order]
      duplicate$sigma <- duplicate$sigma[order]
      duplicate$w <- weights
      return(list(list(p = duplicate), split))
    }
  }
  list(list(p = .split_component(p, widest)), split)
}
