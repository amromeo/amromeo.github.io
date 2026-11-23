# The scoring step has no raw CSVs. The fitting shards retained their minima.
library(RIbench)
working_dir <- commandArgs(TRUE)[1]
sets <- loadTestsetDefinition()
sets <- sets[sets$N == as.integer(Sys.getenv('RIBENCH_N')), ]
files <- list.files(working_dir, pattern='^custom-[0-9]+[.]csv$', full.names=TRUE)
custom <- do.call(rbind, lapply(files, read.csv))
minima <- unique(custom[c('Index','routine_min')])
stopifnot(nrow(minima)==nrow(sets), !anyDuplicated(minima$Index), setequal(minima$Index, sets$Index))
rows <- lapply(c('anchoredRI','anchorsOnly','refineR'), function(method) {
  native <- read.csv(file.path(working_dir, paste0(method,'-all-performance.csv')))
  p <- merge(native,sets[c('Index','N','Distribution','fractionPathol')],by='Index')
  p <- merge(p,minima,by='Index')
  if (method == 'refineR') {
    p$transform <- 'BoxCox'; p$shift <- p$Est_Shift; p$components <- 1L
    p$lrt_p <- p$lrt_statistic <- NA_real_; p$at_bound <- NA; p$anchor_seed <- NA_integer_
  } else {
    d <- custom[custom$algo == method, c('Index','model','shift','components','lrt_p','lrt_statistic','at_bound','anchor_seed')]
    p <- merge(p,d,by='Index'); p$transform <- p$model
  }
  p$method <- method; p$nonp_lambda <- p$GT_Lambda
  p$lower_flag <- p$Analyte != 'CRP' & is.finite(p$LRL) & p$LRL < p$routine_min
  stopifnot(nrow(p)==nrow(sets), !anyDuplicated(p$Index))
  p[c('Analyte','Index','N','Distribution','fractionPathol','nonp_lambda','method','transform',
      'Est_Lambda','shift','components','at_bound','anchor_seed','lrt_statistic','lrt_p',
      'zzDev_LRL','zzDev_URL','zzDevAbs_Ov','LRL','URL','routine_min','lower_flag')]
})
write.csv(do.call(rbind,rows), file.path(working_dir,'per-dataset-diagnostics.csv'), row.names=FALSE)
