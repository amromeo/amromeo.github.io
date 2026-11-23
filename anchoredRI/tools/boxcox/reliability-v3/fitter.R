# Study-local multistart v3. The frozen Box-Cox likelihood and v2.1 starts remain in use.
source('tools/boxcox/reliability-v2/fitter.R')
r3_families <- c('identity','log','shifted_log')
r3_C <- function(s) sqrt(exp(s^2)*expm1(s^2))
r3_rule <- function(d,decimals) {
 delta <- 10^(-decimals)
 stopifnot(length(delta)==1L,is.finite(delta),delta>0,
           max(abs(c(d$value,d$anchor)/delta-round(c(d$value,d$anchor)/delta)))<1e-6)
 blog <- if(d$positive) max(.005,.05*v2_scale(log(d$value))) else .005
 c(delta=delta,identity_floor=.05*d$spread,log_multiplier=r3_C(blog))
}
r3_eval <- function(d,f,k,cs,par,b,rules=numeric(),mode=1L,maxit=1000L,exact_zero=FALSE) {
 spec <- as.double(c(k,d$g,d$n,cs,if(exact_zero)0 else d$location,if(exact_zero)1 else d$spread,
                if(v2_shape(f))NA_real_ else if(f=='identity')1 else 0,if(v2_shift(f))NA_real_ else 0))
 ans<-.Call('round3_native',d$x,d$counts,d$anchors,spec,as.double(par),as.double(b$lo),as.double(b$hi),
            as.double(rules),as.integer(maxit),as.integer(mode))
 names(ans)<-c('par','initial','final','code','message','gradient','initial_par','internal','initial_internal','evaluations')
 ans
}
r3_original_sd <- function(f,p) {
 if(f=='identity')p$sigma else exp(p$mu)*r3_C(p$sigma)
}
r3_required_sd <- function(f,p,rules) {
 med<-if(f=='identity')p$mu+1 else p$shift+exp(p$mu)
 common<-pmax(rules[1],rules[2],pmax(med,0)*rules[3])
 if(f=='shifted_log')pmax(common,exp(p$mu)*r3_C(.005)) else common
}
r3_candidate <- function(d,f,k,index,models,rules) {
 starts<-v2_grid(d,f,k)
 if(k>1L)starts<-c(starts,v2_duplicate(models[[paste(f,k-1L)]]))
 if(f=='shifted_log'&&d$zero)starts[[length(starts)+1L]]<-list(p=models[[paste('log',k)]]$p,label='unshifted log')
 best<-NULL;audit<-list()
 run<-function(cs,b,init,label) {
  z<-r3_eval(d,f,k,cs,init,b,rules); use_start<-is.finite(z$initial)&&(!is.finite(z$final)||z$initial>z$final)
  ll<-if(use_start)z$initial else z$final
  par<-if(use_start)z$initial_par else z$par
  p<-v2_unpack(d,f,k,par,cs)
  critical<-r3_required_sd(f,p,rules);dispersion<-r3_original_sd(f,p)
  feasible<-all(is.finite(c(par,critical,dispersion)))&&all(dispersion>=critical*(1-1e-11))
  if(!feasible)ll<--Inf
  bound<-any(abs(par-b$lo)<1e-5|abs(par-b$hi)<1e-5)||any(abs(dispersion-critical)<=1e-8*pmax(1,critical))
  a<-length(audit)+1L
  audit[[a]]<<-list(label=label,initial=z$initial,final=z$final,retained=ll,code=z$code,message=z$message,
    projected=any(abs(z$initial_par-init)>1e-10),initial_retained=use_start,at_bound=bound)
  if(is.finite(ll)&&(is.null(best)||ll>best$ll))best<<-list(ll=ll,p=p,par=par,cs=cs,box=b,
          origin=a,code=z$code,initial_retained=use_start,at_bound=bound)
 }
 for(s in starts) {
  cs<-if(is.null(s$cs))v2_coords(d,f,s$p) else s$cs
  b<-v2_box(d,f,k,cs); par<-v2_pack(d,f,s$p,cs)
  # The historical K=1 fit has no artificial upper bound below its analytic SD.
  if(s$label=='analytic MLE')b$hi<-pmax(b$hi,par)
  run(cs,b,par,s$label)
 }
 if(is.null(best))stop('No feasible deterministic start for ',f,' K=',k,' Index ',index)
 set.seed(v2_seed(index,f,k),kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
 base<-best;scale<-c(.45*exp(base$par[k+seq_len(k)]),rep(.35,k),rep(.65,d$g*(k-1L)),
                      if(v2_shift(f)).12*diff(d$bounds)/d$spread)
 for(j in 1:20)run(base$cs,base$box,base$par+rnorm(length(scale))*scale,sprintf('seeded perturbation %02d',j))
 best$k<-k;best$transform<-f;best$seed<-v2_seed(index,f,k);best$starts<-audit
 best$bic<--2*best$ll+(2*k+d$g*(k-1L)+v2_shift(f))*log(d$n)
 best
}
r3_baseline <- function(d,index,decimals,kmax=3L,models=list()) {
 rules<-r3_rule(d,decimals)
 for(k in seq_len(kmax))for(f in r3_families) {
  name<-paste(f,k)
  if((f=='log'&&!d$positive)||!is.null(models[[name]]))next
  models[[name]]<-r3_candidate(d,f,k,index,models,rules)
 }
 models
}
r3_zero <- function(d,models,rules) {
 checks<-list()
 for(m in models) {
  if(m$transform!='log')next
  k<-m$k;p<-m$p;cs<-c(0,1)
  par<-v2_pack(d,'log',p,cs)
  # At an excluded zero shift this is an objective identity, separate from fit eligibility.
  lo<-rep(-1e6,length(par)+1L);hi<--lo
  fixed<-r3_eval(d,'log',k,cs,par,list(lo=head(lo,-1),hi=head(hi,-1)),mode=0L,exact_zero=TRUE)
  free<-r3_eval(d,'shifted_log',k,cs,c(par,0),list(lo=lo,hi=hi),mode=0L,exact_zero=TRUE)
  sd0<-r3_required_sd('log',p,rules);med<-exp(p$mu)
  floor_gap<-max(abs(sd0-r3_required_sd('shifted_log',p,rules)),
             abs(sd0-r3_required_sd('identity',list(mu=med-1),rules)))
  err<-abs(fixed$final-free$final);tol<-128*.Machine$double.eps*max(1,d$n,abs(fixed$final),abs(free$final))
  stopifnot(err<=tol,floor_gap<=1e-11*max(1,sd0))
  checks[[length(checks)+1L]]<-data.frame(k=k,zero_eligible=d$zero,fixed_loglik=fixed$final,
       zero_shift_loglik=free$final,absolute_difference=err,roundoff=tol,physical_floor_difference=floor_gap)
 }
 if(length(checks))do.call(rbind,checks) else NULL
}
r3_bt <- function(models,d,rules,probs=c(.025,.975)) {
 tab<-v2_table(models,probs)
 tab$minimum_original_sd<-tab$minimum_sd_to_floor<-NA_real_
 tab$rounding_unit<-unname(rules[1])
 for(i in seq_len(nrow(tab))) {
  m<-models[[paste(tab$transform[i],tab$k[i])]]
  sd<-r3_original_sd(m$transform,m$p);floor<-r3_required_sd(m$transform,m$p,rules)
  tab$minimum_original_sd[i]<-min(sd);tab$minimum_sd_to_floor[i]<-min(sd/floor)
 }
 stopifnot(all(tab$minimum_sd_to_floor>=1-1e-11))
 tab
}

# A seed is expressed in its own saved affine coordinates and box, including the scale-dependent
# shape-model SD limits. Parameter transport does not re-standardize or project a valid donor.
r3_seed <- function(m) {
 list(f=m$transform,k=m$k,par=unname(m$par),cs=unname(m$cs),box=m$box,ll=unname(m$ll),
      code=m$code,source=m$starts[[m$origin]]$label)
}
r3_endpoint <- function(d,s) {
 p<-v2_unpack(d,s$f,s$k,s$par,s$cs)
 p$shift+bc_inverse(p$mu[1]+qnorm(c(.025,.975))*p$sigma[1],p$lambda)
}
r3_check_seed <- function(d,s) {
 z<-r3_eval(d,s$f,s$k,s$cs,s$par,s$box,mode=0L)
 tol<-128*.Machine$double.eps*max(1,d$n,abs(s$ll))
 stopifnot(is.finite(z$final),abs(z$final-s$ll)<=tol,all(abs(z$initial_par-s$par)<=1e-11*pmax(1,abs(s$par))))
 z$final-s$ll
}
r3_continue <- function(d,s,force=FALSE) {
 start_delta<-r3_check_seed(d,s);log<-list()
 if(!force&&s$code==0L)return(list(seed=s,blocks=log,transport_error=start_delta))
 point<-s$par; best<-s$ll;best_par<-point
 tol<-128*.Machine$double.eps*max(1,d$n,abs(best))
 for(i in seq_len(1000L)) {
  if(i==8L||i%%30L==0L){cat('continuation block',i,s$f,s$k,'best log likelihood',format(best,digits=12),'\n');flush.console()}
  z<-r3_eval(d,s$f,s$k,s$cs,point,s$box,maxit=10000L)
  log[[length(log)+1L]]<-list(engine='L-BFGS-B',block=i,code=z$code,initial=z$initial,final=z$final,
                             message=z$message,evaluations=z$evaluations)
  if(is.finite(z$final)&&z$final>best){best<-z$final;best_par<-z$par}
  if(z$code==0L&&is.finite(z$final)&&z$final>=best-tol) {
   s$par<-z$par;s$ll<-z$final;s$code<-z$code
   return(list(seed=s,blocks=log,transport_error=start_delta))
  }
  point<-if(is.finite(z$final)&&z$final>=best-tol)z$par else best_par
  # Try a fresh limited-memory solve first. A repeated failure gets a bounded
  # analytic-gradient bridge; only a genuine final L-BFGS-B code zero is accepted.
  if(i%%2L==0L) {
   evaluate<-function(p)r3_eval(d,s$f,s$k,s$cs,p,s$box,mode=0L)
   bridge<-tryCatch(nlminb(point,function(p)-evaluate(p)$final/d$n,function(p)evaluate(p)$gradient,
                     lower=s$box$lo,upper=s$box$hi,control=list(iter.max=10000L,eval.max=30000L)),
                    error=function(e)list(par=point,convergence=99L,message=conditionMessage(e)))
   w<-evaluate(bridge$par)
   log[[length(log)+1L]]<-list(engine='nlminb bridge',block=i,code=bridge$convergence,initial=z$final,
                       final=w$final,message=bridge$message,evaluations=NA_integer_)
   if(is.finite(w$final)&&w$final>=best-tol){point<-w$par;if(w$final>best){best<-w$final;best_par<-point}}
  }
 }
 stop('Continuation unresolved after 1,000 blocks for ',s$f,' K=',s$k)
}
