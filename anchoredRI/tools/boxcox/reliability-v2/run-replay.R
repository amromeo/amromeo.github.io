# Only the 156-dataset reliability experiment can enter through this driver.
args <- commandArgs(TRUE)
stopifnot(length(args)==8L,args[1] %in% c('local','Actions'))
runtime <- args[1]; data_dir <- args[2]; archive <- args[3]; out <- args[4]
shard <- as.integer(args[5]); shards <- as.integer(args[6]); workers <- as.integer(args[7])
stopifnot(shards>=1,shards<=156,shard>=1,shard<=shards,workers>=1)
library(RIbench)
data.table::setDTthreads(1L)
source('tools/boxcox/inputs.R');source('tools/boxcox/reliability-v2/fitter.R');dyn.load(args[8])
dir.create(out,recursive=TRUE,showWarnings=FALSE)
versions <- capture.output(sessionInfo())
if(runtime=='Actions') {
 old <- readLines(file.path(archive,'versions-1.txt'))
 libs <- function(s)grep('^(BLAS|LAPACK):',s,value=TRUE)
 required <- Sys.getenv('ImageVersion')=='20261116.202.1' && identical(old[1],versions[1]) &&
   identical(libs(old),libs(versions)) && as.character(packageVersion('RIbench'))=='1.0.2' &&
   as.character(packageVersion('data.table'))=='1.16.4' && as.character(packageVersion('digest'))=='0.6.39'
 writeLines(c(paste('C2 runtime matched',required),paste('ImageVersion',Sys.getenv('ImageVersion')),versions),
            file.path(out,'runtime.txt'))
 if(!required)stop('C2 runtime mismatch before fitting')
} else writeLines(versions,file.path(out,'runtime.txt'))
sources <- c('tools/boxcox/models.R','tools/boxcox/inputs.R',
 file.path('tools/boxcox/reliability-v2',c('fitter.R','likelihood.cpp','build.R','run-replay.R')))
source_hash <- vapply(sources,function(f)digest::digest(file=f,algo='sha256'),character(1))
write.csv(data.frame(file=names(source_hash),sha256=unname(source_hash)),file.path(out,'sources.csv'),row.names=FALSE)
hashes <- read.csv(file.path(archive,'anchor-sample-hashes.csv'))
defs <- loadTestsetDefinition()
indices <- (1:156)[((1:156)-1L)%%shards==shard-1L]
fit_one <- function(index) {
 file <- file.path(out,sprintf('fit-%03d.rds',index)); s <- defs[defs$Index==index,]
 stopifnot(nrow(s)==1L,s$Analyte=='Hb')
 d <- bc_ribench_input(s,data_dir,hashes$sha256[match(index,hashes$Index)])
 routine_hash <- digest::digest(d$value,algo='sha256'); anchor_hash <- digest::digest(d$anchor,algo='sha256')
 if(file.exists(file)) {
  result <- readRDS(file)
  stopifnot(identical(result$sources,source_hash),identical(result$runtime,runtime),result$Index==index,
            result$routine_hash==routine_hash,result$anchor_hash==anchor_hash)
 } else {
  when <- proc.time()[3]
  fit <- v2_fit(d$value,d$setting,d$anchor,index); tab <- v2_table(fit,d$probs)
  details <- v2_data(d$value,d$setting,d$anchor)
  result <- list(Index=index,runtime=runtime,routine_hash=routine_hash,anchor_hash=anchor_hash,
    sources=source_hash,table=tab,models=fit,
    meta=data.frame(key=d$key,Index=index,group='RIbench small',subtype='Hb',routine_n=length(d$value),
      anchor_n=length(d$anchor),n=details$n,groups=details$g,positive=details$positive,
      shift_lo=details$bounds[1],shift_hi=details$bounds[2],zero_shift=details$zero),
    seconds=unname(proc.time()[3]-when))
  tmp <- paste0(file,'.tmp-',Sys.getpid());saveRDS(result,tmp);stopifnot(file.rename(tmp,file))
 }
 cat('saved',runtime,index,'\n');flush.console()
 result$models <- NULL;result
}
all <- parallel::mclapply(indices,fit_one,mc.cores=min(workers,length(indices)),mc.preschedule=FALSE,mc.set.seed=FALSE)
if(!all(vapply(all,function(r)is.list(r)&&!is.null(r$Index),logical(1))))stop('One or more dataset fits failed')
tab <- do.call(rbind,lapply(all,function(r) {
 data.frame(key=paste0('r-',r$Index),Index=r$Index,runtime=runtime,source='multistart v2',r$table,
   routine_hash=r$routine_hash,anchor_hash=r$anchor_hash,seconds=r$seconds)
}))
write.csv(tab,file.path(out,'candidates.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(all,`[[`,'meta')),file.path(out,'datasets.csv'),row.names=FALSE)
writeLines(sprintf('%s completed %d distinct datasets',runtime,length(unique(tab$Index))),file.path(out,'DONE.txt'))
