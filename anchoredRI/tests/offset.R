library(anchoredRI)
# Check numerical results before applying the reporting safeguard.
fit_anchored_ri <- function(...) anchoredRI::fit_anchored_ri(..., safeguard=FALSE)
bootstrap_anchored_ri <- function(...) anchoredRI::bootstrap_anchored_ri(..., safeguard=FALSE)

# Check the analytic gradient, including the transformation Jacobian.
native <- getFromNamespace("C_anchored_fit_start", "anchoredRI")
calc_factory <- function(x, a, counts, k) {
  points <- sort(unique(c(x, a)))
  mass <- matrix(0, length(points), ncol(counts))
  mass[match(x, points), ] <- counts
  anchors <- as.double(tabulate(match(a, points), length(points)))
  n <- sum(counts) + length(a)
  spec <- as.double(c(k, ncol(counts), n, 0, 1, 0, 1, 0, NA_real_))
  function(par, gradient = FALSE) {
    z <- .Call(native, points, mass, anchors, spec, par, rep(-Inf, length(par)),
               rep(Inf, length(par)), numeric(), 1000L, 0L)
    if (gradient) z[[6]] else -z[[3]] / n
  }
}
x <- c(.2,.5,1,1.7,3);a <- c(.3,.7,1.2)
counts <- cbind(c(2,4,3,2,1),c(1,2,4,2,3))
for(k in 1:5) {
  par <- c(seq(-.1,.5,length.out=k),rep(log(.6),k),rep(.2,2*(k-1)),-.5)
  calc <- calc_factory(x,a,counts,k)
  numeric_gradient <- vapply(seq_along(par),function(i) {
    plus <- minus <- par;plus[i] <- plus[i]+1e-6;minus[i] <- minus[i]-1e-6
    (calc(plus)-calc(minus))/2e-6
  },numeric(1))
  stopifnot(max(abs(numeric_gradient-calc(par,TRUE)))<1e-5)
  if(k==1L) {
    exact <- sum(dlnorm(a-par[3],par[1],exp(par[2]),log=TRUE))+
      sum(rowSums(counts)*dlnorm(x-par[3],par[1],exp(par[2]),log=TRUE))
    stopifnot(abs(calc(par)+exact/(sum(counts)+length(a)))<1e-10)
  }
}
set.seed(43)
v <- 10+rlnorm(400,log(15),.5);a <- 10+rlnorm(30,log(15),.5)
g <- rep(c("a","b"),each=200)
f <- fit_anchored_ri(v,g,a,"shifted_log",shift_bounds=c(0,12))
scaled <- fit_anchored_ri(10*v,g,10*a,"shifted_log",shift_bounds=c(0,120))
stopifnot(max(abs(scaled$limits/10-f$limits))<.01,
  abs(scaled$shift/10-f$shift)<.01,
  abs(scaled$selected_bic-f$selected_bic-2*430*log(10))<.01)
auto <- fit_anchored_ri(v,g,a,"auto")
stopifnot(nrow(auto$selection)==15L,auto$selected_bic==min(auto$selection$bic_original),
  auto$requested_transform=="auto",auto$transform %in% c("identity","log","shifted_log"))
for(family in c("identity","log")) {
  fixed <- fit_anchored_ri(v,g,a,family)
  all_values <- if(family=="log")log(c(v,a))else c(v,a)
  exact_bic <- -2*sum(dnorm(all_values,mean(all_values),
    sqrt(mean((all_values-mean(all_values))^2)),log=TRUE))+2*log(length(all_values))
  if(family=="log")exact_bic <- exact_bic+2*sum(log(c(v,a)))
  stopifnot(abs(fixed$bic_original[1]-exact_bic)<1e-8)
}
# The public bootstrap preserves the caller RNG and independently refits shifts.
seed_before <- .Random.seed
b <- bootstrap_anchored_ri(v,g,a,"shifted_log",n_boot=3,seed=7,shift_bounds=c(0,12))
stopifnot(identical(seed_before,.Random.seed),b$successful==3L,
  length(unique(round(b$bootstrap$shift,5)))>1L,
  all(b$ci[,1]<=b$ci[,2]),all(b$bootstrap$transform=="shifted_log"))
b2 <- bootstrap_anchored_ri(v,g,a,"shifted_log",n_boot=3,seed=7,
  shift_bounds=c(0,12),workers=if(.Platform$OS.type=="unix")2L else 1L)
stopifnot(isTRUE(all.equal(b$bootstrap,b2$bootstrap)))
fails <- function(expr)inherits(try(force(expr),silent=TRUE),"try-error")
stopifnot(fails(fit_anchored_ri(v,g,a,"auto",shift_bounds=c(0,min(v)+1))),
  fails(fit_anchored_ri(v,g,a,"log",shift_bounds=c(0,12))),
  fails(bootstrap_anchored_ri(v,g,a,n_boot=1)),
  fails(bootstrap_anchored_ri(v,g,a,conf=1)))
