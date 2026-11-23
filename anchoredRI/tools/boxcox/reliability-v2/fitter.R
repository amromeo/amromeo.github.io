# SPDX-License-Identifier: GPL-3.0-only
# Multistart v2 works in Box-Cox coordinates, including x-1 for identity.
source('tools/boxcox/models.R')
v2_families <- c('identity','log','shifted_log','boxcox_unshifted','boxcox_shifted')
v2_shape <- function(f) startsWith(f,'boxcox')
v2_shift <- function(f) f %in% c('shifted_log','boxcox_shifted')
v2_scale <- function(x) { s <- mad(x); if(s>0) s else sd(x) }
v2_data <- function(value,setting,anchor) {
 stopifnot(length(value)==length(setting),length(anchor)>1L,all(is.finite(c(value,anchor))))
 points <- sort(unique(c(value,anchor))); groups <- unique(as.character(setting))
 ix <- match(value,points); g <- length(groups)
 counts <- matrix(as.double(tabulate(ix+(match(setting,groups)-1L)*length(points),length(points)*g)),ncol=g)
 anchors <- as.double(tabulate(match(anchor,points),length(points)))
 spread <- v2_scale(value); stopifnot(is.finite(spread),spread>0,g %in% 1:2)
 b <- min(points)+spread*c(-5,-.001)
 list(value=value,anchor=anchor,x=as.double(points),counts=counts,anchors=anchors,g=g,n=sum(counts)+sum(anchors),
      location=median(value),spread=spread,bounds=b,shifts=unique(c(max(b[1],min(b[2],0)),mean(b),b[1]+.9*diff(b))),
      zero=b[1]<=0&&b[2]>=0,positive=points[1]>0)
}
v2_seed <- function(index,f,k) as.integer(index+1000000L*match(f,v2_families)+100000L*k)
v2_coords <- function(d,f,p) {
 if(f=='identity')return(c(median(d$value)-1,v2_scale(d$value)))
 if(f=='log')return(c(median(log(d$value)),v2_scale(log(d$value))))
 if(f=='shifted_log')return(c(log(d$spread),1))
 y <- bc_transform(d$value-p$shift,p$lambda); c(median(y),v2_scale(y))
}
v2_box <- function(d,f,k,cs) {
 fitted_shape <- v2_shape(f); shifted <- v2_shift(f)
 rx <- range(d$x)
 if(fitted_shape) {
  span <- if(shifted) rx-rev(d$bounds) else rx
  bounds <- range(unlist(lapply(c(-2,0,2),function(l)bc_transform(span,l))))
 } else if(f=='shifted_log') bounds <- log(rx-rev(d$bounds))
 else bounds <- if(f=='identity')rx-1 else log(rx)
 sr <- if(f %in% c('identity','log'))c(.05,20) else c(.005,5)
 b <- (bounds-cs[1])/cs[2]
 wb <- (d$bounds-d$location)/d$spread
 list(lo=c(rep(b[1],k),rep(log(sr[1]),k),rep(-7,d$g*(k-1)),if(fitted_shape)-2,if(shifted)wb[1]),
      hi=c(rep(b[2],k),rep(log(sr[2]),k),rep(7,d$g*(k-1)),if(fitted_shape)2,if(shifted)wb[2]))
}
v2_native <- function(d,f,k,cs,par,b,fit=TRUE) {
 spec <- as.double(c(k,d$g,d$n,cs,d$location,d$spread,
            if(v2_shape(f))NA_real_ else if(f=='identity')1 else 0,if(v2_shift(f))NA_real_ else 0))
 .Call('reliability_fit',d$x,d$counts,d$anchors,spec,as.double(par),as.double(b$lo),as.double(b$hi),fit)
}
v2_pack <- function(d,f,p,cs) {
 k <- length(p$mu)
 eta <- if(k>1L)as.vector(log(sweep(p$w[-k,,drop=FALSE],2,p$w[k,],'/'))) else numeric()
 c((p$mu-cs[1])/cs[2],log(p$sigma/cs[2]),eta,if(v2_shape(f))p$lambda,if(v2_shift(f))(p$shift-d$location)/d$spread)
}
v2_unpack <- function(d,f,k,v,cs) {
 last <- 2*k+d$g*(k-1L)
 eta <- if(k>1L)rbind(matrix(v[2*k+seq_len(d$g*(k-1))],ncol=d$g),0) else matrix(0,1,d$g)
 w <- exp(eta); w <- sweep(w,2,colSums(w),'/')
 list(mu=cs[1]+cs[2]*v[seq_len(k)],sigma=cs[2]*exp(v[k+seq_len(k)]),w=w,
      lambda=if(v2_shape(f))v[last+1L] else if(f=='identity')1 else 0,
      shift=if(v2_shift(f))d$location+d$spread*tail(v,1L) else 0)
}
v2_grid <- function(d,f,k) {
 locations <- if(v2_shift(f))d$shifts else 0
 shapes <- if(v2_shape(f))c(-1,-.5,0,.5,1,2) else if(f=='identity')1 else 0
 slots <- if(k==1)list(numeric()) else if(k==2)list(.15,.5,.85) else
  list(seq(.12,.85,length.out=k-1),seq(.60,.90,length.out=k-1),seq(.10,.40,length.out=k-1))
 starts <- list()
 for(shift in locations) for(lambda in shapes) {
  y <- bc_transform(d$value-shift,lambda); a <- bc_transform(d$anchor-shift,lambda)
  fixed <- f %in% c('identity','log')
  for(q in slots) {
   sy <- if(fixed)v2_scale(y) else sd(y)
   mu <- c(mean(a),unname(quantile(y,q)))
   sigma <- c(sd(a),rep(sy,k-1L))
   if(k==1L && fixed) {
    all <- c(y,a); mu <- mean(all); sigma <- sqrt(mean((all-mu)^2))
   }
   p <- list(mu=mu,sigma=sigma,w=matrix(1/k,k,d$g),lambda=lambda,shift=shift)
   starts[[length(starts)+1L]] <- list(p=p,label=if(k==1L&&fixed)'analytic MLE' else 'grid')
  }
 }
 starts
}
v2_split <- function(p,j,apart=FALSE) {
 k <- length(p$mu)
 q <- p
 q$mu <- c(p$mu,p$mu[j]); q$sigma <- c(p$sigma,p$sigma[j])
 q$w <- rbind(p$w,p$w[j,,drop=FALSE]/2); q$w[j,] <- q$w[j,]/2
 if(apart) {
  if(j!=1L) q$mu[j] <- q$mu[j]-.3*q$sigma[j]
  q$mu[k+1L] <- q$mu[k+1L]+.3*q$sigma[j]
 }
 q
}
# Keep the anchored component first. Relabel only disease components to recover
# an exact embedding if a change of reference log odds would cross a box bound.
v2_duplicate <- function(prev) {
 p <- prev$p; k <- length(p$mu)+1L; widest <- which.max(p$sigma)
 jobs <- list(list(p=v2_split(p,widest,TRUE),label='widest split',cs=prev$cs))
 for(j in unique(c(widest,seq_along(p$mu)))) {
  q <- v2_split(p,j)
  for(last in rev(2:k)) {
   ord <- c(setdiff(seq_len(k),last),last); w <- q$w[ord,,drop=FALSE]
   if(any(abs(log(sweep(w,2,w[k,],'/')))>7+1e-11))next
   q$mu <- q$mu[ord]; q$sigma <- q$sigma[ord]; q$w <- w
   return(c(list(list(p=q,label=if(j==widest)'exact widest duplicate' else 'exact feasible duplicate',cs=prev$cs)),jobs))
  }
 }
 c(list(list(p=v2_split(p,widest),label='widest duplicate outside odds box',cs=prev$cs)),jobs)
}
v2_candidate <- function(d,f,k,index,models) {
 starts <- v2_grid(d,f,k)
 get_fit <- function(family,count=k) models[[paste(family,count)]]
 if(k>1L) starts <- c(starts,v2_duplicate(get_fit(f,k-1L)))
 use <- function(parent,label,translate=FALSE) {
  if(is.null(parent))return()
  p <- parent$p
  if(translate) {
   shift <- if(d$zero)0 else d$shifts[1]
   p$mu <- p$mu+p$shift-shift; p$shift <- shift
  }
  # The free-shape family has the original start-dependent numerical scale.
  cs <- if(label=='unshifted Box-Cox')parent$cs else NULL
  starts[[length(starts)+1L]] <<- list(p=p,label=label,cs=cs)
 }
 if(f=='shifted_log' && d$zero) use(get_fit('log'),'unshifted log')
 if(v2_shape(f)) {
  use(get_fit('identity'),'fixed identity',translate=v2_shift(f))
  if(!v2_shift(f)||d$zero)use(get_fit('log'),'fixed log')
  if(v2_shift(f)) {
   use(get_fit('shifted_log'),'fixed shifted log')
   if(d$zero)use(get_fit('boxcox_unshifted'),'unshifted Box-Cox')
  }
 }
 audit <- list(); best <- NULL
 run <- function(cs,b,init,label) {
  proj <- pmin(b$hi,pmax(b$lo,init)); projected <- any(abs(proj-init)>1e-10)
  ans <- v2_native(d,f,k,cs,proj,b)
  ok0 <- is.finite(ans[[2]]) && abs(ans[[2]])<1e90
  ok1 <- is.finite(ans[[3]]) && abs(ans[[3]])<1e90 && all(is.finite(ans[[1]])) &&
          all(ans[[1]]>=b$lo-1e-12 & ans[[1]]<=b$hi+1e-12)
  use_start <- ok0 && (!ok1||ans[[2]]>ans[[3]])
  ll <- if(use_start)ans[[2]] else if(ok1)ans[[3]] else -Inf
  par <- if(use_start)proj else ans[[1]]
  hit <- any(abs(par-b$lo)<1e-5|abs(par-b$hi)<1e-5)
  a <- length(audit)+1L
  audit[[a]] <<- list(label=label,initial=ans[[2]],final=ans[[3]],retained=ll,
      code=ans[[4]],message=ans[[5]],projected=projected,initial_retained=use_start,at_bound=hit,
      initial_par=proj,final_par=ans[[1]],center=cs[1],scale=cs[2],lo=b$lo,hi=b$hi)
  if(is.finite(ll)&&(is.null(best)||ll>best$ll)) {
   best <<- list(ll=ll,p=v2_unpack(d,f,k,par,cs),cs=cs,par=par,box=b,origin=a,
                code=ans[[4]],initial_retained=use_start,at_bound=hit)
  }
 }
 for(s in starts) {
  cs <- if(is.null(s$cs))v2_coords(d,f,s$p) else s$cs
  b <- v2_box(d,f,k,cs)
  par <- v2_pack(d,f,s$p,cs)
  if(s$label=='analytic MLE') {
   # The frozen one-component normal has an unconstrained closed-form MLE.
   b$lo <- pmin(b$lo,par); b$hi <- pmax(b$hi,par)
  }
  run(cs,b,par,s$label)
 }
 if(is.null(best))stop('No feasible deterministic start for ',f,' K=',k)
 set.seed(v2_seed(index,f,k),kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
 start_best <- best
 # Perturb means in units of each component SD. Odds and log SDs are dimensionless.
 sd_par <- c(.45*exp(start_best$par[k+seq_len(k)]),rep(.35,k),rep(.65,d$g*(k-1)),
             if(v2_shape(f)).30,if(v2_shift(f)).12*diff(d$bounds)/d$spread)
 for(i in 1:20)run(start_best$cs,start_best$box,start_best$par+rnorm(length(sd_par))*sd_par,
                    sprintf('seeded perturbation %02d',i))
 best$k <- k; best$transform <- f; best$seed <- v2_seed(index,f,k); best$starts <- audit
 best$bic <- -2*best$ll+(2*k+d$g*(k-1)+v2_shape(f)+v2_shift(f))*log(d$n)
 best
}
v2_fit <- function(value,setting,anchor,index,k_baseline=3L,k_boxcox=3L) {
 d <- v2_data(value,setting,anchor); models <- list()
 for(k in seq_len(max(k_baseline,k_boxcox))) for(f in v2_families) {
  if((v2_shape(f)&&k>k_boxcox)||(!v2_shape(f)&&k>k_baseline))next
  if(!d$positive && f %in% c('log','boxcox_unshifted'))next
  models[[paste(f,k)]] <- v2_candidate(d,f,k,index,models)
 }
 models
}
v2_table <- function(models,probs) {
 rows <- lapply(models,function(m) {
  p <- m$p; y <- p$mu[1]+qnorm(probs)*p$sigma[1]
  limits <- if(m$transform=='identity')y+1 else p$shift+bc_inverse(y,p$lambda)
  data.frame(method=if(v2_shape(m$transform))'armA' else 'baseline',transform=m$transform,k=m$k,
      eligible=TRUE,bic=m$bic,loglik=m$ll,lambda=p$lambda,shift=p$shift,
      lower=if(length(probs)>1)limits[1] else NA_real_,upper=tail(limits,1),
      origin=m$starts[[m$origin]]$label,code=m$code,initial_retained=m$initial_retained,at_bound=m$at_bound,
      starts=length(m$starts),nonzero_codes=sum(vapply(m$starts,function(s)s$code!=0L,logical(1))),
      perturbation_seed=m$seed)
 })
 x <- do.call(rbind,rows); rownames(x)<-NULL
 x <- x[order(match(x$transform,v2_families),x$k),]; x$selected <- FALSE
 for(m in unique(x$method)) { ix <- which(x$method==m); x$selected[ix[which.min(x$bic[ix])]]<-TRUE }
 x
}
