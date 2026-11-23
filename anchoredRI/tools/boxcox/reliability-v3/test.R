args<-commandArgs(TRUE);stopifnot(length(args)==1L);dyn.load(args[1])
source('tools/boxcox/reliability-v3/fitter.R')
set.seed(301);value<-round(rlnorm(180,2.6,.16),1);anchor<-round(rlnorm(40,2.6,.10),1)
for(g in 1:2) {
 d<-v2_data(value,rep(seq_len(g),length.out=length(value)),anchor);rules<-r3_rule(d,1)
 for(f in r3_families)for(k in 1:3) {
  s<-v2_grid(d,f,k)[[1]]$p; cs<-v2_coords(d,f,s); b<-v2_box(d,f,k,cs)
  q<-v2_pack(d,f,s,cs); q[k+seq_len(k)]<-seq(.18,.62,length.out=k)
  if(v2_shift(f))q[length(q)]<-mean(tail(c(b$lo[length(q)],b$hi[length(q)]),2))
  w<-r3_eval(d,f,k,cs,q,b,rules,mode=2L)
  numerical<-vapply(seq_along(q),function(i) {
   a<-h<-q;a[i]<-a[i]-1e-6;h[i]<-h[i]+1e-6
   -(r3_eval(d,f,k,cs,h,b,rules,mode=2L)$final-r3_eval(d,f,k,cs,a,b,rules,mode=2L)$final)/(2e-6*d$n)
  },numeric(1))
  stopifnot(max(abs(w$gradient-numerical)/pmax(1,abs(numerical)))<3e-5)
  ref<-v2_native(d,f,k,cs,w$par,b,FALSE)
  stopifnot(abs(ref[[2]]-w$final)<1e-8)
  p<-v2_unpack(d,f,k,w$par,cs);stopifnot(all(r3_original_sd(f,p)>=r3_required_sd(f,p,rules)*(1-1e-11)))
 }
}
d<-v2_data(value,rep('pooled',length(value)),anchor); rules<-r3_rule(d,1)
models<-r3_baseline(d,1L,1L,3L);tab<-r3_bt(models,d,rules);zero<-r3_zero(d,models,rules)
stopifnot(nrow(tab)==9,nrow(zero)==3,all(zero$absolute_difference==0))
for(k in 1:3) {
 p<-models[[paste('log',k)]]$p;stopifnot(all(p$sigma>=.005-1e-12))
 for(f in r3_families) {
  m<-models[[paste(f,k)]];stopifnot(sum(startsWith(vapply(m$starts,`[[`,'','label'),'seeded perturbation'))==20L)
 }
}
cat('Passed original-scale floor, objective, analytic-gradient, zero-shift and retained-start tests.\n')
