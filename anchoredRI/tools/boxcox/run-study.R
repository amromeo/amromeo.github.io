# Usage from anchoredRI/: Rscript tools/boxcox/run-study.R config.rds [validate|fit]
library(RIbench)
data.table::setDTthreads(1)
source('tools/boxcox/models.R')
source('tools/boxcox/candidates.R')
source('tools/boxcox/inputs.R')
source('tools/stress-design.R')
args <- commandArgs(TRUE); cfg <- readRDS(args[1]); phase <- args[2]
stopifnot(phase %in% c('validate','fit'),as.character(packageVersion('RIbench'))=='1.0.2')
dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
root <- normalizePath(cfg$output)
halt <- file.path(root,'HALT.txt')
if(file.exists(halt)) stop('Resolve recorded validation failure before resuming')
source_files <- c(file.path('tools/boxcox',c('models.R','inputs.R','capacity-fixed.R',
                           'capacity-shifted.R','candidates.R','run-study.R')),
                  list.files('tools/boxcox/frozen',full.names=TRUE),'tools/stress-design.R')
hashes <- tools::md5sum(source_files)
manifest <- file.path(root,'source-hashes.rds')
if(file.exists(manifest)) stopifnot(identical(readRDS(manifest),hashes)) else saveRDS(hashes,manifest)
writeLines(c(capture.output(sessionInfo()),paste('baseline',cfg$baseline)),file.path(root,'versions.txt'))
t <- loadTestsetDefinition(); t <- t[t$N<=50000,]; stopifnot(nrow(t)==4320)
archives <- list(small=cfg$small_archive,large=cfg$large_archive)
archived <- lapply(archives,function(p)read.csv(file.path(p,'anchoredRI-all-performance.csv')))
custom <- lapply(archives,function(p)do.call(rbind,lapply(list.files(p,pattern='^custom-[0-9]+[.]csv$',full.names=TRUE),read.csv)))
anchor_hashes <- do.call(rbind,lapply(archives,function(p)read.csv(file.path(p,'anchor-sample-hashes.csv'))))
# Verify the full independent-anchor set before any fitting, not just a shard.
actual <- vapply(seq_len(nrow(t)),function(i)digest::digest(bc_anchors(t[i,]),algo='sha256'),character(1))
stopifnot(!anyDuplicated(actual),identical(actual,anchor_hashes$sha256[match(t$Index,anchor_hashes$Index)]))
write.csv(data.frame(Index=t$Index,seed=t$Index+17000000L,sha256=actual),file.path(root,'anchor-hashes.csv'),row.names=FALSE)
keys <- c(paste0('r-',t$Index),sprintf('s-%02d-%02d',rep(1:45,each=20),rep(1:20,45)))
stopifnot(length(keys)==5220,!anyDuplicated(keys))
load_data <- function(key) {
  parts <- strsplit(key,'-',fixed=TRUE)[[1]]
  if(parts[1]=='s') {
    d <- stress_data(as.integer(parts[2]),as.integer(parts[3]))
    d$probs <- c(.025,.975); d$key <- key; return(d)
  }
  s <- t[t$Index==as.integer(parts[2]),]; stage <- if(s$N<=5000)'small' else 'large'
  bc_ribench_input(s, if(stage=='small')cfg$small_data else cfg$large_data,
                   anchor_hashes$sha256[match(s$Index,anchor_hashes$Index)])
}
near <- function(a,b) {
  if (!identical(dim(a),dim(b)) || length(a)!=length(b) || !identical(is.na(a),is.na(b)) ||
      !identical(is.finite(a),is.finite(b))) return(FALSE)
  both <- is.finite(a); infinite <- is.infinite(a)
  all(abs(a[both]-b[both])<=1e-10) && all(a[infinite]==b[infinite])
}
validate <- function(d) {
  frozen <- bc_frozen(TRUE); cap <- bc_capacity(3)
  baseline_seconds <- system.time(old <- frozen$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto'))[['elapsed']]
  now <- cap$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')
  stopifnot(near(now$limits,old$limits),identical(now$components,old$components),
            identical(now$transform,old$transform),near(now$shift,old$shift),
            near(now$selection$bic_original,old$selection$bic_original))
  delta <- 0
  for(lambda in c(1,if(all(c(d$value,d$anchor)>0))0)) {
    a <- bc_fixed(d$value,d$setting,d$anchor,lambda)
    b <- frozen$fit_anchored_ri(d$value,d$setting,d$anchor,transform=if(lambda==1)'identity' else 'log')
    stopifnot(near(a$candidate_limits,b$candidate_limits),near(a$bic_original,b$bic_original))
    delta <- max(delta,a$objective_discrepancy)
  }
  guard <- bc_safeguard(now,d$value,d$anchor,d$probs)
  cross <- NULL
  if(!is.null(d$s)) {
    row <- archived[[d$stage]]; row <- row[row$Index==d$s$Index,]
    expected <- if(length(d$probs)==1)row$URL else c(row$LRL,row$URL)
    stopifnot(nrow(row)==1)
    width <- row$URL-row$LRL
    stopifnot(is.finite(width),width>0)
    differences <- abs(guard$limits-expected)
    c <- custom[[d$stage]]; c <- c[c$Index==d$s$Index&c$algo=='anchoredRI',]
    stopifnot(nrow(c)==1)
    cross <- data.frame(Index=d$s$Index,Analyte=d$s$Analyte,N=d$s$N,
      local_transform=old$transform,archived_transform=c$model,
      local_components=old$components,archived_components=c$components,
      model_changed=old$transform!=c$model || old$components!=c$components,
      archived_width=width,one_sided_span=length(d$probs)==1L,
      lower_local=if(length(expected)==2L) guard$limits[1] else NA_real_,
      lower_archived=if(length(expected)==2L) expected[1] else NA_real_,
      upper_local=tail(guard$limits,1),upper_archived=tail(expected,1),
      max_absolute_difference=max(differences),
      max_fraction_of_width=max(differences)/width,
      outside_tolerance=any(differences>.001*width))
  }
  # Candidate parameters are diagnostic additions from the exactly matching copy.
  old$candidate_parameters <- now$candidate_parameters
  candidates <- bc_collect_candidates(now,'anchoredRI',d$probs,all(c(d$value,d$anchor)>0),d$key)
  gap <- bc_runner_up(candidates)
  list(key=d$key,baseline=old,safeguard=guard,baseline_seconds=baseline_seconds,
       candidate_table=candidates,bic_gap=gap,
       cross_machine=cross,max_loglik_discrepancy=delta)
}
validation_dir <- file.path(root,'validation'); fit_dir <- file.path(root,'fits')
dir.create(validation_dir,showWarnings=FALSE); dir.create(fit_dir,showWarnings=FALSE)
if(phase=='fit') {
  stopifnot(file.exists(file.path(root,'VALIDATED.txt')),
    all(file.exists(file.path(validation_dir,paste0(keys,'.rds')))))
}
worker <- function(key) {
  if(file.exists(halt)) return(FALSE)
  out <- file.path(if(phase=='validate')validation_dir else fit_dir,paste0(key,'.rds'))
  if(file.exists(out)) return(TRUE)
  tryCatch({
    d <- load_data(key)
    if(phase=='validate') result <- validate(d) else {
      old <- readRDS(file.path(validation_dir,paste0(key,'.rds')))
      e <- bc_capacity(3)
      seeds <- lapply(c('identity',if(all(c(d$value,d$anchor)>0))'log','shifted_log'),
        function(tr)e$fit_anchored_ri(d$value,d$setting,d$anchor,transform=tr))
      b <- bc_capacity(5)
      methods <- list(armA=function()bc_fit(d$value,d$setting,d$anchor,seeds),
                      armB=function()b$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto'))
      result <- list(key=key,s=d$s,case=d$case,replicate=d$replicate,routine_min=min(d$value),
                     anchor_hash=digest::digest(d$anchor,algo='sha256'))
      for(method in names(methods)) {
        elapsed <- system.time(f <- tryCatch(methods[[method]](),error=function(e)list(error=conditionMessage(e))))[['elapsed']]
        guard <- if(is.null(f$error))bc_safeguard(f,d$value,d$anchor,d$probs) else NULL
        candidates <- bc_collect_candidates(f,method,d$probs,all(c(d$value,d$anchor)>0),key)
        result[[method]] <- list(fit=f,safeguard=guard,seconds=elapsed,
                                 candidate_table=candidates,bic_gap=bc_runner_up(candidates))
      }
    }
    saveRDS(result,paste0(out,'.tmp')); stopifnot(file.rename(paste0(out,'.tmp'),out))
    cat(phase,key,'complete\n');flush.console(); TRUE
  },error=function(e) {
    writeLines(paste(key,conditionMessage(e)),halt)
    cat('HALT',key,conditionMessage(e),'\n');flush.console(); FALSE
  })
}
start <- Sys.time()
# Batches bound the amount of work that can continue after a failure.
for(begin in seq(1,length(keys),by=cfg$workers)) {
  batch <- keys[begin:min(length(keys),begin+cfg$workers-1)]
  outcome <- parallel::mclapply(batch,worker,mc.cores=cfg$workers,mc.set.seed=FALSE)
  if(file.exists(halt)||!all(vapply(outcome,isTRUE,logical(1)))) stop('Study halted; inspect HALT.txt')
}
if(phase=='validate') {
  cross <- do.call(rbind,lapply(keys,function(key)
    readRDS(file.path(validation_dir,paste0(key,'.rds')))$cross_machine))
  stopifnot(nrow(cross)==4320)
  source('tools/boxcox/export-candidates.R')
  bc_export_tables(root)
  write.csv(cross,file.path(root,'cross-machine-checks.csv'),row.names=FALSE)
  write.csv(cross[cross$model_changed|cross$outside_tolerance,],
            file.path(root,'cross-machine-exceptions.csv'),row.names=FALSE)
  summary <- data.frame(datasets=nrow(cross),model_changes=sum(cross$model_changed),
    model_change_fraction=mean(cross$model_changed),outside=sum(cross$outside_tolerance),outside_fraction=mean(cross$outside_tolerance),
    max_absolute_difference=max(cross$max_absolute_difference),
    max_fraction_of_width=max(cross$max_fraction_of_width))
  write.csv(summary,file.path(root,'cross-machine-summary.csv'),row.names=FALSE)
  if(sum(cross$model_changed)>0.05*nrow(cross)) {
    writeLines('Cross-machine review gate triggered. See cross-machine-exceptions.csv.',halt)
    stop('More than 5% of the completed RIbench audit selected a different model')
  }
}
write.csv(data.frame(phase=phase,start=as.character(start),finish=as.character(Sys.time()),
                    elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs'))),
          file.path(root,paste0(phase,'-runtime.csv')),row.names=FALSE)
writeLines(paste(length(keys),'datasets',phase,'completed'),
           file.path(root,if(phase=='validate')'VALIDATED.txt' else 'FITTED.txt'))
