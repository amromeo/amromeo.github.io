# Run in the C2 GitHub Actions environment. No reference-truth scoring is used here.
# Inputs: mode, original data artifact, scored C2 artifact, output, shard, number of shards.
library(RIbench)
data.table::setDTthreads(1L)
source('tools/boxcox/models.R'); source('tools/boxcox/inputs.R'); source('tools/boxcox/candidates.R')
args <- commandArgs(TRUE)
stopifnot(length(args)==6L,args[1] %in% c('preflight','armA'))
mode <- args[1]; data_dir <- args[2]; archive <- args[3]; out <- args[4]
shard <- as.integer(args[5]); shards <- as.integer(args[6]); stopifnot(shard>=1,shard<=shards)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
versions <- capture.output(sessionInfo()); original <- readLines(file.path(archive,'versions-1.txt'))
lib_lines <- function(lines) grep('^(BLAS|LAPACK):',lines,value=TRUE)
image <- Sys.getenv('ImageVersion'); runtime_ok <- image=='20261116.202.1' &&
  identical(original[1],versions[1]) && identical(lib_lines(original),lib_lines(versions)) &&
  as.character(packageVersion('RIbench'))=='1.0.2' && as.character(packageVersion('data.table'))=='1.16.4' &&
  as.character(packageVersion('digest'))=='0.6.39'
writeLines(c(paste('C2 image 20261116.202.1; observed',image),paste('Runtime matched:',runtime_ok),
             paste('Git commit',Sys.getenv('GITHUB_SHA')),versions), file.path(out,'runtime.txt'))
if(!runtime_ok) stop('C2 runner image, R, package versions or numerical libraries differ. See runtime.txt.')
t <- loadTestsetDefinition(); hashes <- read.csv(file.path(archive,'anchor-sample-hashes.csv'))
all_indices <- if(mode=='preflight')c(100L,148L) else 1:156
indices <- all_indices[(seq_along(all_indices)-1L) %% shards == shard-1L]
model_sources <- file.path('tools/boxcox',c('capacity-fixed.R','capacity-shifted.R','models.R','inputs.R','candidates.R'))
md5 <- tools::md5sum(model_sources);write.csv(data.frame(file=names(md5),md5=unname(md5)),file.path(out,'model-sources.csv'),row.names=FALSE)
candidates <- summaries <- original_matches <- list()
for(index in indices) {
  s <- t[t$Index==index,]; stopifnot(nrow(s)==1,s$Analyte=='Hb')
  d <- bc_ribench_input(s,data_dir,hashes$sha256[match(index,hashes$Index)])
  input_hash <- digest::digest(d$anchor,algo='sha256')
  method <- if(mode=='preflight') 'anchoredRI' else 'armA'
  if(mode=='preflight') {
    fit <- bc_capacity(3L)$fit_anchored_ri(d$value,d$setting,d$anchor,'auto')
    frozen <- bc_frozen(TRUE)$fit_anchored_ri(d$value,d$setting,d$anchor,'auto')
    stopifnot(identical(fit$selection,frozen$selection),identical(fit$limits,frozen$limits))
    saved <- new.env();file <- sprintf('%d_%s_seed_%s_anchoredRI.Rdata',s$Index,s$Analyte,s$startSeed)
    load(file.path(archive,'Results','anchoredRI',s$Analyte,file),envir=saved)
    original_d <- do.call(rbind,lapply(list.files(archive,pattern='^custom-[0-9]+[.]csv$',full.names=TRUE),read.csv))
    original_d <- original_d[original_d$Index==index & original_d$algo=='anchoredRI',];stopifnot(nrow(original_d)==1)
    original_matches[[length(original_matches)+1L]] <- data.frame(Index=index,original_transform=original_d$model,
      original_k=original_d$components,replay_transform=frozen$transform,replay_k=frozen$components,
      original_lower=saved$fit$PointEst[1],original_upper=saved$fit$PointEst[2],
      replay_lower=unname(frozen$limits[1]),replay_upper=unname(frozen$limits[2]),
      max_difference=max(abs(frozen$limits-saved$fit$PointEst)))
  } else {
    # Checkpoints from this separate repeat never enter the main scoring run.
    fit <- tryCatch(bc_fit(d$value,d$setting,d$anchor),error=function(e)list(error=conditionMessage(e)))
  }
  tab <- bc_collect_candidates(fit,method,d$probs,all(c(d$value,d$anchor)>0),d$key)
  gap <- bc_runner_up(tab); best <- match(1L,tab$rank)
  candidates[[length(candidates)+1L]] <- tab
  summaries[[length(summaries)+1L]] <- data.frame(Index=index,analyte=s$Analyte,N=s$N,
    transform=tab$transform[best],components=tab$components[best],lambda=tab$lambda[best],shift=tab$shift[best],
    lower=tab$lower[best],upper=tab$upper[best],bic=gap$first_bic,bic_gap=gap$gap,
    anchor_seed=17000000L+index,anchor_sha256=input_hash,routine_hash=s$HashCode)
  if(mode=='armA')saveRDS(fit,file.path(out,sprintf('armA-%03d.rds',index)))
  cat('Completed',mode,index,'\n');flush.console()
}
write.csv(do.call(rbind,candidates),file.path(out,'candidates.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(out,'selected.csv'),row.names=FALSE)
if(length(original_matches))write.csv(do.call(rbind,original_matches),file.path(out,'replay-vs-original.csv'),row.names=FALSE)
