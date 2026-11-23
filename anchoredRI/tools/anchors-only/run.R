source('tools/stress-design.R');source('tools/anchors-only/estimators.R')
out<-'assets/anchors-only';dir.create(out,recursive=TRUE,showWarnings=FALSE)
paths<-c(sort(list.files('R',full.names=TRUE)),'tools/stress-design.R','tools/anchors-only/PROTOCOL.md','tools/anchors-only/estimators.R','tools/anchors-only/run.R')
hash<-tools::md5sum(paths)
archive<-'assets/comparison-study/summary/comparison-results.zip'
z<-unz(archive,'all-results.csv');prior<-read.csv(z)
prior<-prior[prior$method=='anchored40',];stopifnot(nrow(prior)==900L,!anyDuplicated(prior[c('case','replicate')]))
rows<-list();diagnostics<-list();cache<-new.env();i<-0L
for(case in 1:45)for(rep in 1:20) {
 d<-stress_data(case,rep);old<-prior[prior$case==case&prior$replicate==rep,]
 stopifnot(nrow(old)==1L,d$seed==old$seed,max(abs(d$truth-c(old$true_lower,old$true_upper)))<1e-10,d$actual_affected==old$actual_affected)
 key<-paste(d$family,rep,sep='-')
 if(!exists(key,cache,inherits=FALSE)) {
  begin<-proc.time()[['elapsed']];oracle<-tryCatch(ao_oracle(d$anchor,d$family),error=identity)
  t1<-proc.time()[['elapsed']]-begin;begin<-proc.time()[['elapsed']];auto<-tryCatch(ao_auto(d$anchor),error=identity)
  assign(key,list(anchor=d$anchor,oracle=oracle,auto=auto,elapsed=c(oracle=t1,auto=proc.time()[['elapsed']]-begin)),cache)
 }
 fits<-get(key,cache);stopifnot(identical(fits$anchor,d$anchor))
 for(method in c('oracle','auto')) {
  f<-fits[[method]];ok<-!inherits(f,'error')&&stress_valid(f$limits)
  r<-stress_row(d,'direct40',if(ok)f$limits else c(NA_real_,NA_real_),
    model=if(ok)f$transform else '',shift=if(ok)f$shift else NA_real_,components=1L,
    at_bound=if(ok)f$at_bound else NA,elapsed=fits$elapsed[[method]],
    error=if(inherits(f,'error'))conditionMessage(f)else if(!ok)f$error else '')
  r$method<-paste0('anchors_',method);i<-i+1L;rows[[i]]<-r
  if(method=='auto'&&!inherits(f,'error'))diagnostics[[length(diagnostics)+1L]]<-cbind(case=case,replicate=rep,family=d$family,anchor_sample=key,f$diagnostics)
 }
}
stopifnot(length(ls(cache))==60L,identical(hash,tools::md5sum(paths)))
write.csv(do.call(rbind,rows),file.path(out,'estimates.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,diagnostics),file.path(out,'candidates.csv'),row.names=FALSE,na='')
write.csv(prior,file.path(out,'anchoredRI-matched.csv'),row.names=FALSE,na='')
write.csv(data.frame(path=names(hash),md5=unname(hash)),file.path(out,'source-hashes.csv'),row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),paste('archive MD5:',tools::md5sum(archive))),file.path(out,'versions.txt'))
cat('Completed 900 paired evaluations; 60 unique anchor samples; 1800 estimator rows.\n')
