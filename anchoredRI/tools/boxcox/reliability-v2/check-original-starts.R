# Compare the original Arm A deterministic starts with the saved v2.1 starts.
# Fitted endpoints, selected solutions and reference limits do not enter this check.
arg <- commandArgs(TRUE);stopifnot(length(arg)==5L,arg[1] %in% c('local','Actions'))
runtime <- arg[1]; original <- arg[2]; improved <- arg[3];config <- readRDS(arg[4]);dest <- arg[5]
source('tools/boxcox/inputs.R')
defs <- RIbench::loadTestsetDefinition(); hashes <- read.csv(file.path(config$output,'anchor-hashes.csv'))
compare <- function(a,b) {
 stopifnot(length(a)==length(b));max(abs(a-b)/pmax(1,abs(a),abs(b)))
}
out <- list()
for(i in 1:156) {
 if(runtime=='local')old <- readRDS(file.path(original,sprintf('r-%d.rds',i)))$armA$fit
 else {
  fs<-list.files(original,sprintf('^armA-%03d[.]rds$',i),full.names=TRUE,recursive=TRUE)
  stopifnot(length(fs)==1L);old<-readRDS(fs)
 }
 new<-readRDS(file.path(improved,sprintf('fit-%03d.rds',i)))$models
 s<-defs[defs$Index==i,];d<-bc_ribench_input(s,config$small_data,hashes$sha256[match(i,hashes$Index)])
 spread<-mad(d$value);if(spread<=0)spread<-sd(d$value);location<-median(d$value)
 for(c in old$candidates) {
  family<-if(c$shifted)'boxcox_shifted' else 'boxcox_unshifted';k<-c$k
  incoming<-new[[paste(family,k)]]$starts
  newtags<-vapply(incoming,`[[`,'','label')
  oldtags<-vapply(c$starts,function(a)if(a$start$source=='grid')'grid' else paste('original seed',a$start$source),'')
  for(tag in unique(c('grid',oldtags,grep('^original seed ',newtags,value=TRUE)))) {
   aa<-c$starts[oldtags==tag];bb<-incoming[newtags==tag]; matched<-min(length(aa),length(bb))
   errors<-vapply(seq_len(matched),function(j) {
    x<-aa[[j]]; y<-bb[[j]];p<-x$start; eta<-if(k==1L)numeric() else as.vector(log(sweep(p$w[-k,,drop=FALSE],2,p$w[k,],'/')))
    par<-c((p$mu-x$center)/x$scale,log(pmin(4.9,pmax(.006,p$sigma/x$scale))),eta,p$lambda,if(c$shifted)p$shift)
    par<-pmax(x$lower,pmin(x$upper,par))
    actual<-y$initial_par
    if(c$shifted)actual[length(actual)]<-location+spread*tail(actual,1)
    max(compare(par,actual),compare(c(x$center,x$scale),c(y$center,y$scale)))
   },numeric(1))
   out[[length(out)+1L]]<-data.frame(Index=i,runtime=runtime,transform=family,k=k,source=tag,
       original=length(aa),improved=length(bb),paired=matched,within_1e_10=sum(errors<=1e-10),
       max_scaled_difference=if(matched)max(errors) else NA_real_)
  }
 }
}
write.csv(do.call(rbind,out),dest,row.names=FALSE)
