# SPDX-License-Identifier: GPL-3.0-only
# All study code is isolated from the installed/package estimator.
bc_root <- 'tools/boxcox'
bc_frozen <- function(one_group = FALSE) {
  e <- new.env(parent=globalenv())
  prefix <- if(one_group) 'c2-' else ''
  for(f in c('fit_anchored_ri.R','fit_anchored_shifted.R'))
    sys.source(file.path(bc_root,'frozen',paste0(prefix,f)),e)
  e
}
bc_capacity <- function(max_components=5L) {
  e <- new.env(parent=globalenv()); e$STUDY_K <- max_components
  sys.source(file.path(bc_root,'capacity-fixed.R'),e)
  sys.source(file.path(bc_root,'capacity-shifted.R'),e)
  e
}
bc_transform <- function(x, lambda) {
  if(lambda == 1) return(x-1)
  t <- log(x)
  if(abs(lambda)<1e-7) return(t+lambda*t^2/2+lambda^2*t^3/6)
  expm1(lambda*t)/lambda
}
bc_lambda_derivative <- function(x,lambda) {
  t <- log(x)
  if(abs(lambda)<1e-5) return(t^2/2+lambda*t^3/3+lambda^2*t^4/8)
  (exp(lambda*t)*(lambda*t-1)+1)/lambda^2
}
bc_inverse <- function(y,lambda) {
  if(abs(lambda)<1e-12) return(exp(y))
  z <- 1+lambda*y
  out <- rep(NA_real_,length(y)); ok <- is.finite(z)&z>0
  out[ok] <- exp(log1p(lambda*y[ok])/lambda)
  out
}
# Compressed routine values retain their multiplicities and setting labels.
bc_counts <- function(value,setting) {
  x <- sort(unique(value)); groups <- unique(as.character(setting))
  counts <- matrix(tabulate(match(value,x)+(match(as.character(setting),groups)-1L)*length(x),
                            length(x)*length(groups)),ncol=length(groups))
  list(x=x,counts=counts,groups=groups)
}
# The objective uses original x units. center/scale are fixed numerical
# coordinates within one optimization, never fitted parameters.
bc_objective <- function(x,a,counts,k,center=0,scale=1,fixed_lambda=NULL,fixed_shift=NULL) {
  G <- ncol(counts); total <- rowSums(counts); ng <- colSums(counts)
  n <- sum(total)+length(a); base <- 2*k+G*(k-1L)
  last <- answer <- NULL
  function(par,gradient=FALSE) {
    if(!identical(par,last)) {
      lambda <- if(is.null(fixed_lambda)) par[base+1L] else fixed_lambda
      d <- if(is.null(fixed_shift)) tail(par,1L) else fixed_shift
      mu <- par[seq_len(k)]; sigma <- exp(par[k+seq_len(k)])
      eta <- if(k==1L) matrix(0,1,G) else rbind(matrix(par[2*k+seq_len(G*(k-1L))],ncol=G),0)
      w <- exp(sweep(eta,2,apply(eta,2,max),'-')); w <- sweep(w,2,colSums(w),'/')
      t <- x-d; ta <- a-d
      identity_case <- !is.null(fixed_lambda)&&lambda==1&&!is.null(fixed_shift)
      if(!identity_case && any(c(t,ta)<=0)) stop('Box-Cox domain violation')
      y <- (bc_transform(t,lambda)-center)/scale
      ya <- (bc_transform(ta,lambda)-center)/scale
      jac <- if(identity_case) 0 else (lambda-1)*(sum(total*log(t))+sum(log(ta)))
      ll <- sum(dnorm(ya,mu[1],sigma[1],log=TRUE))+jac-n*log(scale)
      density <- vapply(seq_len(k),function(j)dnorm(y,mu[j],sigma[j],log=TRUE),numeric(length(y)))
      mass <- matrix(0,length(y),k); ws <- matrix(0,k-1L,G)
      for(h in seq_len(G)) {
        parts <- sweep(density,2,log(w[,h]),'+'); peak <- apply(parts,1,max)
        scaled <- exp(parts-peak); denom <- rowSums(scaled)
        ll <- ll+sum(counts[,h]*(peak+log(denom)))
        mh <- scaled*(counts[,h]/denom); mass <- mass+mh
        if(k>1) ws[,h] <- (colSums(mh)-ng[h]*w[,h])[-k]
      }
      u <- sweep(outer(y,mu,'-'),2,sigma,'/'); ua <- (ya-mu[1])/sigma[1]
      ms <- colSums(mass*u)/sigma; ss <- colSums(mass*(u^2-1))
      ms[1] <- ms[1]+sum(ua)/sigma[1]; ss[1] <- ss[1]+sum(ua^2-1)
      dy <- -rowSums(sweep(mass*u,2,sigma,'/'))/scale
      da <- -ua/sigma[1]/scale
      score <- c(ms,ss,ws)
      if(is.null(fixed_lambda)) score <- c(score,
        sum(dy*bc_lambda_derivative(t,lambda))+sum(da*bc_lambda_derivative(ta,lambda))+
          sum(total*log(t))+sum(log(ta)))
      if(is.null(fixed_shift)) score <- c(score,
        -sum(dy*t^(lambda-1))-sum(da*ta^(lambda-1))-
          (lambda-1)*(sum(total/t)+sum(1/ta)))
      last <<- par+0
      answer <<- list(value=-ll/n,gradient=-score/n,
        parameters=list(mu=center+scale*mu,sigma=scale*sigma,w=w,lambda=lambda,shift=d),
        loglik=ll)
    }
    if(identical(gradient,'all')) return(answer)
    if(gradient) answer$gradient else answer$value
  }
}

# Fixed-shape endpoints use the frozen optimization in affine coordinates.
# The independent general Box-Cox objective is checked at every fitted candidate.
bc_fixed <- function(value,setting,anchor,lambda) {
  stopifnot(lambda %in% c(0,1))
  e <- bc_capacity(3L)
  f <- e$.fit_anchored_fixed(value,setting,anchor,if(lambda==1)'identity' else 'log')
  count <- bc_counts(value,setting); n <- length(value)+length(anchor)
  discrepancies <- numeric(3)
  for(k in 1:3) {
    p <- f$candidate_parameters[[k]]
    if(is.null(p)) next
    mu <- p$mu-if(lambda==1)1 else 0
    eta <- if(k==1) numeric() else as.vector(log(sweep(p$w[-k,,drop=FALSE],2,p$w[k,],'/')))
    calc <- bc_objective(count$x,anchor,count$counts,k,fixed_lambda=lambda,fixed_shift=0)
    ll <- -n*calc(c(mu,log(p$sigma),eta))
    expected <- (length(c(mu,log(p$sigma),eta))*log(n)-f$bic_original[k])/2
    discrepancies[k] <- abs(ll-expected)
  }
  stopifnot(all(discrepancies <= n*1e-8))
  f$objective_discrepancy <- discrepancies; f
}

bc_fit <- function(value,setting,anchor,seed_fits=NULL) {
  count <- bc_counts(value,setting); n <- length(value)+length(anchor); G <- ncol(count$counts)
  spread <- mad(value); if(spread<=0) spread <- sd(value)
  stopifnot(is.finite(spread),spread>0)
  bounds <- min(c(value,anchor))+spread*c(-5,-.001)
  shifts <- unique(c(min(bounds[2],max(bounds[1],0)),mean(bounds),bounds[1]+.9*diff(bounds)))
  lambdas <- c(-1,-.5,0,.5,1,2)
  if(is.null(seed_fits)) {
    e <- bc_capacity(3L)
    seed_fits <- lapply(c('identity',if(all(c(value,anchor)>0))'log','shifted_log'),
      function(tr)e$fit_anchored_ri(value,setting,anchor,transform=tr))
  }
  candidates <- list()
  for(shifted in c(FALSE,TRUE)) {
    if(!shifted && any(c(value,anchor)<=0)) next
    for(k in 1:3) {
      slots <- if(k==1) list(numeric()) else if(k==2) list(.15,.5,.85) else
        list(c(.12,.85),c(.60,.90),c(.10,.40))
      starts <- list()
      for(d in if(shifted)shifts else 0) for(lambda in lambdas) {
        y <- bc_transform(value-d,lambda); a <- bc_transform(anchor-d,lambda)
        for(probs in slots) starts[[length(starts)+1L]] <- list(lambda=lambda,shift=d,
          mu=c(mean(a),unname(quantile(y,probs))),
          sigma=c(sd(a),rep(sd(y),k-1L)),w=matrix(1/k,k,G),source='grid')
      }
      for(f in seed_fits) {
        p <- f$candidate_parameters[[k]]; if(is.null(p)) next
        lambda <- if(f$transform=='identity')1 else 0
        d <- f$candidate_shift[k]
        if((!shifted && d!=0)||(shifted && (d<bounds[1]||d>bounds[2]))) next
        starts[[length(starts)+1L]] <- list(lambda=lambda,shift=d,
          mu=p$mu-if(lambda==1)1 else 0,sigma=p$sigma,w=p$w,source=f$transform)
      }
      audit <- lapply(starts,function(s) {
        y <- bc_transform(value-s$shift,s$lambda)
        center <- median(y); scale <- mad(y); if(scale<=0) scale <- sd(y)
        if(!is.finite(scale)||scale<=0) return(list(error='Invalid numerical scale',start=s))
        # Wide mean bounds cover every observed value across the shape domain.
        corners <- unlist(lapply(c(-2,0,2),function(l)
          bc_transform(range(c(value,anchor))-if(shifted)rep(rev(bounds),each=1) else 0,l)))
        lo <- c(rep((min(corners)-center)/scale,k),rep(log(.005),k),rep(-7,G*(k-1)),-2,
                if(shifted)bounds[1])
        hi <- c(rep((max(corners)-center)/scale,k),rep(log(5),k),rep(7,G*(k-1)),2,
                if(shifted)bounds[2])
        eta <- if(k==1)numeric() else as.vector(log(sweep(s$w[-k,,drop=FALSE],2,s$w[k,],'/')))
        init <- c((s$mu-center)/scale,log(pmin(4.9,pmax(.006,s$sigma/scale))),eta,s$lambda,
                  if(shifted)s$shift)
        calc <- bc_objective(count$x,anchor,count$counts,k,center,scale,
                             fixed_shift=if(shifted)NULL else 0)
        fit <- tryCatch(optim(init,calc,function(p)calc(p,TRUE),method='L-BFGS-B',lower=lo,upper=hi,
                             control=list(maxit=700L,factr=1e6)),error=function(e)list(error=conditionMessage(e)))
        fit$start <- s; fit$lower <- lo; fit$upper <- hi; fit$center <- center; fit$scale <- scale
        if(is.null(fit$error)) {
          fit$parameters <- calc(fit$par,'all')$parameters
          fit$loglik <- -n*fit$value
          fit$at_bound <- any(abs(fit$par-lo)<1e-5|abs(fit$par-hi)<1e-5)
        }
        fit
      })
      valid <- which(vapply(audit,function(f)is.null(f$error)&&identical(f$convergence,0L)&&
                              is.finite(f$value),logical(1)))
      candidate <- list(k=k,shifted=shifted,starts=audit,bic=Inf)
      if(length(valid)) {
        best <- audit[[valid[which.min(vapply(audit[valid],function(f)f$value,numeric(1)))]]]
        candidate$best <- best
        candidate$bic <- -2*best$loglik+(2*k+G*(k-1)+1+as.integer(shifted))*log(n)
      }
      candidates[[length(candidates)+1L]] <- candidate
    }
  }
  scores <- vapply(candidates,function(f)f$bic,numeric(1))
  if(!any(is.finite(scores))) stop('No Box-Cox candidate converged')
  best <- candidates[[which.min(scores)]]; p <- best$best$parameters
  list(transform='boxcox',lambda=p$lambda,shift=p$shift,components=best$k,
       parameters=p,weights=p$w,selected_bic=best$bic,at_bound=best$best$at_bound,
       limits=p$shift+bc_inverse(p$mu[1]+qnorm(c(.025,.975))*p$sigma[1],p$lambda),
       candidates=candidates)
}

bc_safeguard <- function(f,value,anchor,probs=c(.025,.975)) {
  if(f$transform=='boxcox') {
    p <- f$parameters; a <- bc_transform(anchor-p$shift,p$lambda)
    limits <- p$shift+bc_inverse(p$mu[1]+qnorm(probs)*p$sigma[1],p$lambda)
  } else {
    p <- f$candidate_parameters[[f$components]]
    a <- if(f$transform=='identity')anchor else log(anchor-f$shift)
    y <- p$mu[1]+qnorm(probs)*p$sigma[1]
    limits <- if(f$transform=='identity')y else f$shift+exp(y)
  }
  own_sd <- sqrt(mean((a-mean(a))^2))
  statistic <- max(0,2*(sum(dnorm(a,mean(a),own_sd,log=TRUE))-
                         sum(dnorm(a,p$mu[1],p$sigma[1],log=TRUE))))
  pv <- pchisq(statistic,2,lower.tail=FALSE)
  lower <- length(probs)==2L && is.finite(limits[1]) && limits[1]<min(value)
  technical <- !is.finite(pv)||any(!is.finite(limits))
  decline <- !technical&&(pv<.01||lower)
  list(limits=limits,reported=if(technical||decline)rep(NA_real_,length(probs))else limits,
       lrt_statistic=statistic,lrt_p=pv,lower_flag=lower,declined=decline,technical=technical,
       reason=if(technical)'check or interval unavailable' else paste(c(
         if(pv<.01)'routine data disagree with the screened sample',
         if(lower)'estimated lower limit is below every observed result'),collapse='; '))
}
