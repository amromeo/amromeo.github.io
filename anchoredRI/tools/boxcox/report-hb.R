# Use local frozen C2 BICs for all candidates. The adapter adds candidate endpoints.
source('tools/boxcox/models.R'); source('tools/boxcox/candidates.R'); source('tools/boxcox/inputs.R')
cfg <- readRDS(commandArgs(TRUE)[1]); dest <- 'tools/boxcox/validation'
design <- RIbench::loadTestsetDefinition(); data.table::setDTthreads(1)
hashes <- read.csv(file.path(cfg$small_archive,'anchor-sample-hashes.csv'))
tables <- lapply(c(100L,148L),function(index) {
  s <- design[design$Index==index,]; stopifnot(s$Analyte=='Hb')
  d <- bc_ribench_input(s,cfg$small_data,hashes$sha256[match(index,hashes$Index)])
  f <- bc_capacity(3L)$fit_anchored_ri(d$value,d$setting,d$anchor,'auto')
  frozen <- bc_frozen(TRUE)$fit_anchored_ri(d$value,d$setting,d$anchor,'auto')
  stopifnot(identical(f$selection,frozen$selection),identical(f$limits,frozen$limits))
  bc_collect_candidates(f,'anchoredRI',d$probs,all(c(d$value,d$anchor)>0),d$key)
})
write.csv(do.call(rbind,tables),file.path(dest,'hb-100-148-local-candidates.csv'),row.names=FALSE)
write.csv(do.call(rbind,lapply(tables,bc_runner_up)),file.path(dest,'hb-100-148-local-gaps.csv'),row.names=FALSE)
