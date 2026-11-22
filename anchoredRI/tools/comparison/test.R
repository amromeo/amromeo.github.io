source('tools/stress-design.R')
for(p in sort(list.files('R',full.names=TRUE,pattern='[.]R$')))source(p)
source('tools/comparison/mixture.R')
for(G in 1:2)for(A in c(0L,5L))for(k in 1:3) {
 x<-seq(.1,3,length.out=30);counts<-matrix(rep(1:3,length.out=30*G),ncol=G)
 a<-if(A)seq(.3,1,length.out=A)else numeric()
 fn<-.comparison_shift_objective(x,a,counts,k)
 par<-c(seq(.1,.5,length.out=k),rep(log(.7),k),rep(.1,G*(k-1)),-1)
 num<-vapply(seq_along(par),function(j){lo<-hi<-par;lo[j]<-lo[j]-1e-6;hi[j]<-hi[j]+1e-6;(fn(hi)-fn(lo))/2e-6},numeric(1))
 stopifnot(max(abs(fn(par,TRUE)-num))<1e-6)
}
for(case in 1:3) {
 d<-stress_data(case,2,n_routine=500)
 original<-fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')
 generalized<-fit_comparison_mixture(d$value,d$setting,d$anchor,transform='auto')
 stopifnot(isTRUE(all.equal(original$limits,generalized$limits,tolerance=1e-10)),
   isTRUE(all.equal(original$selection,generalized$selection,tolerance=1e-10)))
 pooled<-fit_comparison_mixture(d$value,rep('pooled',length(d$value)),d$anchor,transform='auto')
 unanchored<-fit_comparison_mixture(d$value,rep('pooled',length(d$value)),numeric(),transform='auto')
 stopifnot(stress_valid(pooled$limits),stress_valid(comparison_identify(unanchored)$limits),
   ncol(pooled$weights)==1L,ncol(unanchored$weights)==1L)
 stopifnot(isTRUE(all.equal(pooled$limits,fit_comparison_mixture(d$value,rep('renamed',length(d$value)),d$anchor,transform='auto')$limits)))
 cat('Passed frozen equivalence and pooled fits for case',case,'\n')
}
cat('All comparison tests passed.\n')
