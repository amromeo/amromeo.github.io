# Read checkpoints and evaluate saved parameters. This script makes no optimizer calls.
args <- commandArgs(TRUE)
stopifnot(length(args)==4)
cfg <- readRDS(args[1]); replay <- args[2]; preflight <- args[3]; out <- args[4]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('tools/boxcox/models.R')
source('tools/boxcox/inputs.R')
source('tools/stress-design.R')
data.table::setDTthreads(1L)
defs <- RIbench::loadTestsetDefinition()
defs <- defs[defs$N<=50000,]
hash <- read.csv(file.path(cfg$output,'anchor-hashes.csv'))
keys <- c(paste0('r-',defs$Index),sprintf('s-%02d-%02d',rep(1:45,each=20),rep(1:20,45)))
replay_sets <- list.files(replay,'^armA-[0-9]+[.]rds$',recursive=TRUE,full.names=TRUE)
replay_index <- as.integer(sub('armA-([0-9]+)[.]rds','\\1',basename(replay_sets)))
stopifnot(length(keys)==5220,identical(sort(replay_index),1:156))
oldfits <- list.files(file.path(cfg$output,'fits'),pattern='[.]rds$')
constants <- list(); pieces <- list(); pos <- 0L
pack <- function(tab, key, method, runtime, source) {
 data.frame(key=key,method=method,runtime=runtime,source=source,transform=tab$transform,
            k=tab$components,eligible=tab$eligible,bic=tab$bic,stringsAsFactors=FALSE)
}
for(i in seq_along(keys)) {
 key <- keys[i]
 if(i<=nrow(defs)) {
  s <- defs[i,]; idx <- s$Index
  d <- bc_ribench_input(s,if(s$N<=5000)cfg$small_data else cfg$large_data,
                        hash$sha256[match(idx,hash$Index)])
  group <- if(s$N<=5000)'RIbench small' else 'RIbench 50000'
  subtype <- s$Analyte
 } else {
  v <- as.integer(strsplit(key,'-',fixed=TRUE)[[1]][-1])
  d <- stress_data(v[1],v[2]); idx <- NA_integer_; group <- 'Simulation'
  subtype <- as.character(stress_design$family[v[1]])
 }
 a <- d$anchor; x <- d$value; g <- length(unique(d$setting)); n <- length(x)+length(a)
 spread <- mad(x); if(spread<=0) spread <- sd(x)
 bounds <- min(c(x,a))+spread*c(-5,-.001)
 constants[[i]] <- data.frame(key=key,Index=idx,group=group,subtype=subtype,
   n=n,groups=g,routine_n=length(x),anchor_n=length(a),positive=min(c(x,a))>0,
   shift_lo=bounds[1],shift_hi=bounds[2],zero_shift=bounds[1]<=0 && bounds[2]>=0)
 old <- readRDS(file.path(cfg$output,'validation',paste0(key,'.rds')))
 pos <- pos+1L; pieces[[pos]] <- pack(old$candidate_table,key,'baseline','local','validated')
 if(paste0(key,'.rds') %in% oldfits) {
  part <- readRDS(file.path(cfg$output,'fits',paste0(key,'.rds')))
  for(method in c('armA','armB')) {
   pos <- pos+1L; pieces[[pos]] <- pack(part[[method]]$candidate_table,key,method,'local','checkpoint')
  }
 }
 if(!is.na(idx) && idx %in% replay_index) {
  f <- readRDS(replay_sets[match(idx,replay_index)])
  count <- bc_counts(x,d$setting); seeds <- list()
  for(candidate in f$candidates) {
   k <- candidate$k
   tr <- if(candidate$shifted)'boxcox_shifted' else 'boxcox_unshifted'
   pos <- pos+1L
   pieces[[pos]] <- data.frame(key=key,method='armA',runtime='Actions',source='replay',
       transform=tr,k=k,eligible=TRUE,bic=candidate$bic)
   for(start in candidate$starts) {
    p <- start$start; tr <- p$source
    if(!tr %in% c('identity','log','shifted_log')) next
    tag <- paste(tr,k); if(!is.null(seeds[[tag]])) next
    objective <- bc_objective(count$x,a,count$counts,k,fixed_lambda=p$lambda,fixed_shift=p$shift)
    eta <- if(k==1L)numeric() else as.vector(log(sweep(p$w[-k,,drop=FALSE],2,p$w[k,],'/')))
    ll <- -n*objective(c(p$mu,log(p$sigma),eta))
    seeds[[tag]] <- data.frame(key=key,method='baseline',runtime='Actions',source='saved Arm A seed',
      transform=tr,k=k,eligible=TRUE,bic=-2*ll+(2*k+g*(k-1L)+as.integer(tr=='shifted_log'))*log(n))
   }
  }
  for(tr in c('identity','log','shifted_log')) for(k in 1:3) {
   tag <- paste(tr,k)
   if(is.null(seeds[[tag]])) seeds[[tag]] <- data.frame(key=key,method='baseline',runtime='Actions',
     source='seed not retained',transform=tr,k=k,eligible=TRUE,bic=NA_real_)
  }
  pos <- pos+1L; pieces[[pos]] <- do.call(rbind,seeds)
 }
}
meta <- do.call(rbind,constants)
pre <- read.csv(file.path(preflight,'candidates.csv'))
pos <- pos+1L
pieces[[pos]] <- data.frame(key=pre$key,method='baseline',runtime='Actions preflight',
  source='separate Hb replay',transform=pre$transform,k=pre$components,eligible=pre$eligible,bic=pre$bic)
all <- do.call(rbind,pieces); stopifnot(!anyDuplicated(all[c('key','method','runtime','transform','k')]))
write.csv(meta,file.path(out,'datasets.csv'),row.names=FALSE)
write.csv(all,file.path(out,'saved-bics.csv'),row.names=FALSE)
