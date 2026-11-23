args <- commandArgs(TRUE); stopifnot(length(args)==1L);dyn.load(args[1])
source('tools/boxcox/reliability-v2/fitter.R')
set.seed(4731)
x <- round(rlnorm(150,.2,.45),2); a <- round(rlnorm(32,.1,.3),2)
for(g in 1:2) {
 setting <- rep(seq_len(g),length.out=length(x)); d <- v2_data(x,setting,a);c <- bc_counts(x,setting)
 for(f in v2_families) for(k in 1:3) {
  lambda <- if(v2_shape(f)).31 else if(f=='identity')1 else 0
  shift <- if(v2_shift(f))d$bounds[1]+.7*diff(d$bounds) else 0
  p <- list(mu=seq(-.2,.4,length.out=k),sigma=seq(.6,.9,length.out=k),w=matrix(1/k,k,g),lambda=lambda,shift=shift)
  cs <- v2_coords(d,f,p); b <- v2_box(d,f,k,cs); init <- v2_pack(d,f,p,cs)
  result <- v2_native(d,f,k,cs,init,b,FALSE)
  ref <- bc_objective(c$x,a,c$counts,k,cs[1],cs[2],
    fixed_lambda=if(v2_shape(f))NULL else lambda,fixed_shift=if(v2_shift(f))NULL else 0)
  refpar <- init; if(v2_shift(f))refpar[length(refpar)] <- shift
  expected <- ref(refpar,'all'); grad <- expected$gradient
  if(v2_shift(f))grad[length(grad)] <- grad[length(grad)]*d$spread
  num <- vapply(seq_along(init),function(i) {
   hi <- lo <- init; hi[i]<-hi[i]+1e-6;lo[i]<-lo[i]-1e-6
   -(v2_native(d,f,k,cs,hi,b,FALSE)[[2]]-v2_native(d,f,k,cs,lo,b,FALSE)[[2]])/(2e-6*d$n)
  },numeric(1))
  stopifnot(abs(result[[2]]-expected$loglik)<1e-8,
          max(abs(result[[6]]-grad))<1e-9,max(abs(result[[6]]-num))<2e-6)
 }
}
models <- v2_fit(x,rep('p',length(x)),a,1L,2L,2L)
rows <- v2_table(models,c(.025,.975))
stopifnot(nrow(rows)==10L,sum(rows$selected)==2L)
for(m in models) {
 stopifnot(sum(startsWith(vapply(m$starts,`[[`,'','label'),'seeded perturbation'))==20L)
 witness <- vapply(m$starts,`[[`,0.0,'initial')
 stopifnot(m$ll+1e-9>=max(witness))
}
one <- models[['identity 1']]
stopifnot(abs(one$ll-sum(dnorm(c(x,a),mean(c(x,a)),sqrt(mean((c(x,a)-mean(c(x,a)))^2)),log=TRUE)))<1e-9)
xneg <- x-4; dn <- v2_data(xneg,rep('p',length(x)),a-4)
id <- v2_candidate(dn,'identity',1,1,list())
stopifnot(!dn$positive,is.finite(id$ll))
cat('Passed native/R likelihood, gradient, analytic identity, negative values, multi-setting, nesting-start and 20-perturbation checks.\n')

# Regression for the original interior SD starts and original fixed-fit seeds.
d <- v2_data(x,rep('p',length(x)),a)
p <- list(mu=c(0,1),sigma=c(1e-12,1e10),w=matrix(c(.5,.5),2L),lambda=.5,shift=0)
for(f in c('boxcox_shifted','boxcox_unshifted','shifted_log')) {
 cs <- v2_coords(d,f,p)
 stopifnot(max(abs(v2_original_scales(p,d,f)$sigma/cs[2]-c(.006,4.9)))<1e-12)
}
for(f in c('identity','log'))stopifnot(abs(v2_original_scales(p,d,f)$sigma[1]/v2_coords(d,f,p)[2]-.06)<1e-12)
for(k in 1:2) {
 tags <- vapply(models[[paste('boxcox_unshifted',k)]]$starts,`[[`,'','label')
 stopifnot(all(c('original seed identity','original seed log','fixed identity','fixed log') %in% tags))
}
cat('Passed legacy deterministic start and frozen seed checks.\n')
