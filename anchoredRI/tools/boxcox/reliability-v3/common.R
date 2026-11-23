library(RIbench)
data.table::setDTthreads(1L)
source('tools/boxcox/inputs.R')
source('tools/boxcox/reliability-v3/fitter.R')
r3_sources <- function() {
 files<-c('tools/boxcox/inputs.R','tools/boxcox/models.R','tools/boxcox/capacity-fixed.R','tools/boxcox/capacity-shifted.R',
   file.path('tools/boxcox/reliability-v2',c('fitter.R','likelihood.cpp')),
   file.path('tools/boxcox/reliability-v3',c('native.cpp','fitter.R','common.R','build.R','replay.R')))
 vapply(files,function(f)digest::digest(file=f,algo='sha256'),character(1))
}
r3_pedigree <- function(runtime,archive,out) {
 session<-capture.output(sessionInfo()); matched<-NA
 if(runtime=='Actions') {
  old<-readLines(file.path(archive,'versions-1.txt')); libs<-function(x)grep('^(BLAS|LAPACK):',x,value=TRUE)
  matched<-Sys.getenv('ImageVersion')=='20261116.202.1'&&identical(old[1],session[1])&&identical(libs(old),libs(session))&&
   as.character(packageVersion('RIbench'))=='1.0.2'&&as.character(packageVersion('data.table'))=='1.16.4'&&
   as.character(packageVersion('digest'))=='0.6.39'
 }
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 writeLines(c(paste('runtime',runtime,'C2 image match',matched),paste('ImageVersion',Sys.getenv('ImageVersion')),session),
            file.path(out,'runtime.txt'))
 if(runtime=='Actions'&&!matched)stop('Archived C2 runtime mismatch.')
 hashes<-r3_sources()
 write.csv(data.frame(file=names(hashes),sha256=unname(hashes)),file.path(out,'sources.csv'),row.names=FALSE)
 hashes
}
r3_inputs <- function(s,data_dir,hashes) {
 d<-bc_ribench_input(s,data_dir,hashes$sha256[match(s$Index,hashes$Index)])
 data<-v2_data(d$value,d$setting,d$anchor)
 list(d=data,routine=digest::digest(d$value,algo='sha256'),anchor=digest::digest(d$anchor,algo='sha256'),probs=d$probs)
}
r3_save <- function(result,file) {
 temp<-paste0(file,'.tmp-',Sys.getpid());saveRDS(result,temp);stopifnot(file.rename(temp,file))
}
r3_nesting <- function(models,d) {
 rows<-list()
 for(m in models) {
  add<-function(nested,kind,eligible=TRUE) {
   ll<-if(is.null(nested))NA_real_ else nested$ll;tol<-128*.Machine$double.eps*max(1,d$n,abs(m$ll),abs(ll),na.rm=TRUE)
   status<-if(!eligible)'structurally excluded' else if(!is.finite(ll)||!is.finite(m$ll))'unavailable' else 'checked'
   rows[[length(rows)+1L]]<<-data.frame(transform=m$transform,k=m$k,relationship=kind,loglik=m$ll,
      nested_loglik=ll,loss=ll-m$ll,roundoff=tol,status=status,violation=status=='checked'&&ll-m$ll>tol)
  }
  if(m$k>1L)add(models[[paste(m$transform,m$k-1L)]],'adjacent K')
  if(m$transform=='shifted_log')add(models[[paste('log',m$k)]],'shifted vs zero',d$positive&&d$zero)
 }
 do.call(rbind,rows)
}
r3_npar <- function(s,g) 2*s$k+g*(s$k-1L)+as.integer(v2_shape(s$f))+as.integer(v2_shift(s$f))
r3_winner <- function(seeds,d) {
 bic<-vapply(seeds,function(s)-2*s$ll+r3_npar(s,d$g)*log(d$n),numeric(1))
 ordered<-order(bic,vapply(seeds,function(s)match(s$f,v2_families),integer(1)),
                     vapply(seeds,function(s)as.integer(s$k),integer(1)))
 names(seeds)[ordered[1]]
}
r3_row <- function(seed,d) {
 lim<-r3_endpoint(d,seed)
 data.frame(transform=seed$f,k=seed$k,ll=seed$ll,bic=-2*seed$ll+r3_npar(seed,d$g)*log(d$n),
            lower=lim[1],upper=lim[2],code=seed$code)
}
r3_traces <- function(trace,tag) {
 if(!length(trace))return(NULL)
 x<-do.call(rbind,lapply(trace,as.data.frame,stringsAsFactors=FALSE));x$tag<-tag;x
}
r3_pool <- function(ids,workers,fun) {
 res<-parallel::mclapply(ids,fun,mc.cores=min(length(ids),workers),mc.preschedule=FALSE,mc.set.seed=FALSE)
 if(length(res)!=length(ids)||!all(vapply(res,isTRUE,logical(1))))stop('Unfinished jobs. See the worker log; existing checkpoints are reusable.')
}
# Parameter and likelihood transport use round-trip binary64 text, including the original numerical box.
r3_exchange <- function(result) {
 rows<-lapply(names(result$seeds),function(name) {
  s<-result$seeds[[name]];x<-data.frame(Index=result$Index,runtime=result$runtime,transform=s$f,k=s$k,
    selected=name==result$winner,code=s$code,n=length(s$par),routine_hash=result$routine,anchor_hash=result$anchor,
    center=sprintf('%.17g',s$cs[1]),scale=sprintf('%.17g',s$cs[2]),loglik=sprintf('%.17g',s$ll),
    lower=sprintf('%.17g',result$arm_tab[name,'lower']),upper=sprintf('%.17g',result$arm_tab[name,'upper']))
  for(field in c('par','lo','hi')) {
   v<-if(field=='par')s$par else s$box[[field]]
   for(j in 1:10)x[[sprintf('%s%02d',field,j)]]<-if(j<=length(v))sprintf('%.17g',v[j])else ''
  }
  x
 })
 do.call(rbind,rows)
}
r3_read_exchange <- function(path) {
 files<-sort(list.files(path,'^boxcox_(shifted|unshifted)-K[1-3][.]csv$',full.names=TRUE))
 stopifnot(length(files)==6L)
 x<-do.call(rbind,lapply(files,read.csv,colClasses='character',check.names=FALSE));x$Index<-as.integer(x$Index)
 x$k<-as.integer(x$k);x$selected<-as.logical(x$selected);x$lower<-as.double(x$lower);x$upper<-as.double(x$upper)
 stopifnot(nrow(x)==936L,all(table(x$Index)==6),identical(sort(unique(x$Index)),1:156),sum(x$selected)==156L,
           !anyDuplicated(x[c('Index','transform','k')]))
 x
}
r3_from_row <- function(r) {
 stopifnot(nrow(r)==1L)
 take<-function(f)as.double(unlist(r[sprintf('%s%02d',f,seq_len(as.integer(r$n)))],use.names=FALSE))
 list(f=r$transform,k=r$k,par=take('par'),cs=c(as.double(r$center),as.double(r$scale)),
      box=list(lo=take('lo'),hi=take('hi')),ll=as.double(r$loglik),code=as.integer(r$code),source='transferred fit')
}
