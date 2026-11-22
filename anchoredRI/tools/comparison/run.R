# Rscript tools/comparison/run.R CASE [workers=4] [reps=20]
source('tools/stress-design.R')
for(p in sort(list.files('R',pattern='[.]R$',full.names=TRUE)))source(p)
source('tools/comparison/mixture.R')
a<-commandArgs(trailingOnly=TRUE);case<-as.integer(a[1]);workers<-if(length(a)>1)as.integer(a[2])else 4L
reps<-if(length(a)>2)as.integer(a[3])else 20L
stopifnot(case%in%1:45,reps%in%1:20,workers>=1,
 as.character(packageVersion('refineR'))=='2.0.0')
outdir<-'assets/comparison-study';dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
paths<-c(sort(list.files('R',full.names=TRUE)), 'tools/stress-design.R',sort(list.files('tools/comparison',pattern='[.]R$',full.names=TRUE)))
hash<-tools::md5sum(paths)
write.csv(data.frame(path=names(hash),md5=unname(hash)),file.path(outdir,sprintf('hashes-%02d.csv',case)),row.names=FALSE)
run_one<-function(rep) {
 d<-stress_data(case,rep);pool<-rep('pooled',length(d$value));rows<-list()
 make<-function(method,f,elapsed,error='') {
  ok<-!inherits(f,'error')&&stress_valid(f$limits)
  row<-stress_row(d,'anchored40',if(ok)f$limits else c(NA_real_,NA_real_),
   model=if(ok)f$transform else '',shift=if(ok&&!is.null(f$shift))f$shift else NA_real_,
   components=if(ok&&!is.null(f$components))f$components else NA_integer_,
   at_bound=if(ok&&!is.null(f$at_bound))f$at_bound else NA,
   error=if(inherits(f,'error'))conditionMessage(f)else if(!ok)'Invalid limits'else error,elapsed=elapsed)
  row$method<-method;row$n_screened<-if(method%in%c('anchored40','pooled40'))40L else 0L
  row$identified_component<-if(ok&&!is.null(f$identified_component))f$identified_component else NA_integer_
  row
 }
 for(method in c('refineR','gmm_largest','pooled40','anchored40')) {
  path<-file.path(outdir,sprintf('fit-%02d-%02d-%s.rds',case,rep,method))
  if(file.exists(path)) {cached<-readRDS(path);stopifnot(identical(cached$hash,hash));rows[[method]]<-cached$rows;next}
  t0<-proc.time()[['elapsed']]
  f<-tryCatch(switch(method,
   refineR={fit<-refineR::findRI(d$value,model=stress_design$refine_model[case],NBootstrap=0,seed=1200000L+case*100L+rep,NCores=1)
    list(limits=refineR::getRI(fit,pointEst='fullDataEst',UCMargins=FALSE)$PointEst,transform=stress_design$refine_model[case])},
   gmm_largest=fit_comparison_mixture(d$value,pool,numeric(),transform='auto'),
   pooled40=fit_comparison_mixture(d$value,pool,d$anchor,transform='auto'),
   anchored40=fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')),error=identity)
  elapsed<-proc.time()[['elapsed']]-t0
  if(method=='gmm_largest'&&!inherits(f,'error')) {
   rr<-rbind(make('gmm_largest',comparison_identify(f,'largest'),elapsed),make('gmm_lowest',comparison_identify(f,'lowest'),elapsed))
  } else {
   rr<-make(method,f,elapsed)
   if(method=='gmm_largest')rr<-rbind(rr,transform(rr,method='gmm_lowest'))
  }
  stopifnot(identical(hash,tools::md5sum(paths)))
  saveRDS(list(hash=hash,rows=rr,fit=f),path)
  write.csv(rr,sub('[.]rds$','.csv',path),row.names=FALSE,na='')
  rows[[method]]<-rr
  message(sprintf('case %02d replicate %02d %s complete (%.1f sec)',case,rep,method,elapsed))
 }
 for(n in c(40L,120L,400L)) {
  row<-stress_row(d,paste0('direct',n),quantile(head(d$screened,n),c(.025,.975),type=7,names=FALSE),model='Empirical quantiles (type 7)')
  row$identified_component<-NA_integer_;rows[[paste0('direct',n)]]<-row
 }
 do.call(rbind,rows)
}
results<-parallel::mclapply(seq_len(reps),run_one,mc.cores=workers,mc.set.seed=FALSE)
stopifnot(all(vapply(results,is.data.frame,logical(1))))
write.csv(do.call(rbind,results),file.path(outdir,sprintf('results-%02d.csv',case)),row.names=FALSE,na='')
writeLines(c(capture.output(sessionInfo()),paste('refineR',packageVersion('refineR'))),file.path(outdir,sprintf('versions-%02d.txt',case)))
