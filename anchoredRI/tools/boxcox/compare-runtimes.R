# Compare before-safeguard local Arm A fits with the 156 GitHub Actions repeats.
# Usage: Rscript tools/boxcox/compare-runtimes.R config.rds downloaded-shards/
args <- commandArgs(TRUE);cfg <- readRDS(args[1]);archive <- args[2]
source('tools/boxcox/candidates.R')
files <- list.files(archive,pattern='^selected[.]csv$',recursive=TRUE,full.names=TRUE)
stopifnot(length(files)>0)
source_hashes <- readRDS(file.path(cfg$output,'source-hashes.rds'))
for(p in dirname(files)) {
  s <- read.csv(file.path(p,'model-sources.csv')); runtime <- readLines(file.path(p,'runtime.txt'))
  stopifnot(nrow(s)==5L,!anyDuplicated(s$file),identical(s$md5,unname(source_hashes[s$file])),
            'Runtime matched: TRUE' %in% runtime)
}
other <- do.call(rbind,lapply(files,read.csv));stopifnot(nrow(other)==156L,!anyDuplicated(other$Index),setequal(other$Index,1:156))
local <- do.call(rbind,lapply(1:156,function(i) {
  obj <- readRDS(file.path(cfg$output,'fits',paste0('r-',i,'.rds')));tab <- obj$armA$candidate_table
  best <- match(1L,tab$rank)
  data.frame(Index=i,transform=tab$transform[best],components=tab$components[best],lambda=tab$lambda[best],
    shift=tab$shift[best],lower=tab$lower[best],upper=tab$upper[best],bic=tab$bic[best],
    anchor_sha256=obj$anchor_hash,routine_hash=obj$s$HashCode)
}))
x <- merge(local,other,by='Index',suffixes=c('_local','_actions'))
stopifnot(all(x$anchor_sha256_local==x$anchor_sha256_actions),all(x$routine_hash_local==x$routine_hash_actions))
model_pairs <- !is.na(x$components_local)&!is.na(x$components_actions)
x$model_changed <- ifelse(model_pairs,x$transform_local!=x$transform_actions|x$components_local!=x$components_actions,NA)
x$lambda_difference <- x$lambda_local-x$lambda_actions;x$shift_difference <- x$shift_local-x$shift_actions
x$width_actions <- x$upper_actions-x$lower_actions
x$lower_difference <- x$lower_local-x$lower_actions; x$upper_difference <- x$upper_local-x$upper_actions
x$max_limit_difference <- pmax(abs(x$lower_difference),abs(x$upper_difference))
x$fraction_of_width <- ifelse(x$width_actions>0,x$max_limit_difference/x$width_actions,NA_real_)
finite <- is.finite(x$fraction_of_width); x$outside_tolerance <- ifelse(finite,x$fraction_of_width>.001,NA)
maximum <- function(z,valid)if(any(valid))max(z[valid]) else NA_real_
summary <- data.frame(prespecified=156,model_pairs=sum(model_pairs),model_changes=sum(x$model_changed,na.rm=TRUE),
  limit_pairs=sum(finite),outside_tolerance=sum(x$outside_tolerance,na.rm=TRUE),
  largest_absolute=maximum(x$max_limit_difference,finite),largest_fraction=maximum(x$fraction_of_width,finite))
out <- file.path(cfg$output,'reproducibility');dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(x,file.path(out,'armA-reproducibility.csv'),row.names=FALSE)
write.csv(x[which(!model_pairs|!finite|x$model_changed|x$outside_tolerance),],file.path(out,'exceptions.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
section <- c('## Arm A reproducibility across runtimes','',
  sprintf('Among %d/%d datasets with a selected model in both runtimes, %d selected a different component count or shift specification.',summary$model_pairs,summary$prespecified,summary$model_changes),
  sprintf('For %d pairs with finite limits, %d exceeded 0.1%% of the Actions interval width. The largest difference was %.6g Hb units (largest relative difference %.3f%%).',summary$limit_pairs,summary$outside_tolerance,summary$largest_absolute,100*summary$largest_fraction))
writeLines(section,file.path(out,'REPORT.md'))
# Either the local scoring or the Actions repeat can finish first.
main <- file.path(cfg$output,'REPORT.md')
if(file.exists(main)) {
  lines <- readLines(main)
  marker <- which(lines %in% c(section[1],'The paired 156-dataset Arm A reproducibility results are pending.'))
  stopifnot(length(marker)==1L)
  temp <- tempfile('repro-report-',tmpdir=dirname(main))
  writeLines(c(lines[seq_len(marker-1L)],section),temp)
  stopifnot(file.rename(temp,main))
}
