library(RIbench)
source('tools/boxcox/models.R')
source('tools/boxcox/candidates.R')
source('tools/boxcox/export-candidates.R')
source('tools/stress-design.R')
cfg <- readRDS(commandArgs(TRUE)[1]); root <- cfg$output
stopifnot(file.exists(file.path(root,'FITTED.txt')),!file.exists(file.path(root,'HALT.txt')))
bc_export_tables(root,include_arms=TRUE)
t <- loadTestsetDefinition(); t <- t[t$N<=50000,]
methods <- c('armA','armB','anchoredRI','anchorsOnly','refineR')
diagnostics <- list(); summaries <- list(); pairs <- list(); flags <- list()
for(stage in c('small','large')) {
  sets <- t[if(stage=='small')t$N<=5000 else t$N==50000,]
  archive <- if(stage=='small')cfg$small_archive else cfg$large_archive
  w <- file.path(root,stage); dir.create(w,recursive=TRUE,showWarnings=FALSE)
  for(algo in c('anchorsOnly','refineR')) {
    dest <- file.path(w,'Results');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(file.path(archive,'Results',algo),dest,recursive=TRUE))
  }
  diag_stage <- list()
  for(i in seq_len(nrow(sets))) {
    s <- sets[i,]; key <- paste0('r-',s$Index)
    r <- readRDS(file.path(root,'fits',paste0(key,'.rds')))
    baseline <- readRDS(file.path(root,'validation',paste0(key,'.rds')))
    r$anchoredRI <- list(fit=baseline$baseline,safeguard=baseline$safeguard,seconds=baseline$baseline_seconds)
    for(method in c('armA','armB','anchoredRI')) {
      f <- r[[method]]$fit; g <- r[[method]]$safeguard
      probs <- if(s$Analyte=='CRP').95 else c(.025,.975)
      limits <- if(is.null(g))rep(NA_real_,length(probs)) else g$limits
      { # Every anchored baseline is the same-runtime frozen-source refit.
        fit <- data.frame(Percentile=probs,PointEst=limits,Runtime=c(r[[method]]$seconds,rep(NA_real_,length(probs)-1)))
        dest <- file.path(w,'Results',method,s$Analyte);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
        save(fit,file=file.path(dest,paste0(s$Index,'_',s$Analyte,'_seed_',s$startSeed,'_',method,'.Rdata')))
      }
      val <- function(x,default=NA)if(is.null(x))default else x
      diag_stage[[length(diag_stage)+1]] <- data.frame(Index=s$Index,Analyte=s$Analyte,N=s$N,
        method=method,transform=val(f$transform,''),lambda=val(f$lambda),shift=val(f$shift),
        components=val(f$components),at_bound=val(f$at_bound),anchor_seed=17000000L+s$Index,
        anchor_hash=r$anchor_hash,routine_min=r$routine_min,lrt_statistic=val(g$lrt_statistic),
        lrt_p=val(g$lrt_p),lower_flag=val(g$lower_flag,FALSE),declined=val(g$declined,FALSE),
        technical=val(g$technical,TRUE),reason=val(g$reason,val(f$error,'')),
        seconds=r[[method]]$seconds)
    }
  }
  diag_stage <- do.call(rbind,diag_stage)
  all <- list(); scored <- list()
  for(method in methods) {
    markers <- unique(sets$Analyte)
    res <- RIbench:::readResultFilesAll(analytes=markers,algo=method,baseDir=w,tableTCs=sets)
    ri <- if(inherits(res[[1]][[1]],'RWDRI'))
      RIbench:::computeRIsAll(analytes=markers,algo=method,resIn=res,tableTCs=sets,truncNormal=FALSE) else
      RIbench:::getRIsAllwithoutModel(analytes=markers,algo=method,resIn=res,tableTCs=sets)
    perf <- RIbench:::computePerfMeasAll(analytes=markers,algo=method,risIn=ri,tableTCs=sets,cutoffZ=5)
    all[[method]] <- perf; p <- do.call(rbind,perf$resRIs)
    p <- merge(p,sets[,c('Index','Distribution','fractionPathol')],by='Index')
    p$reported <- TRUE
    if(method %in% c('armA','armB','anchoredRI')) {
      d <- diag_stage[diag_stage$method==method,]
      p <- merge(p,d[,setdiff(names(d),c('Analyte','N'))],by='Index')
      p$reported <- !p$declined & !p$technical
    }
    p$stratum <- paste(p$Distribution,ifelse(p$fractionPathol<=.2,'low <=20%','high >20%'))
    scored[[method]] <- p
    write.csv(p,file.path(w,paste0(method,'-performance.csv')),row.names=FALSE)
    if(method %in% c('armA','armB','anchoredRI')) {
      error <- p$zzDevAbs_Ov
      groups <- list(all=rep(TRUE,nrow(p)),below_0.25=is.finite(error)&error<.25,
                     above_5=is.finite(error)&error>5,unscored=!is.finite(error))
      for(group in names(groups)) {
        q <- p[groups[[group]],]
        flags[[length(flags)+1L]] <- data.frame(stage=stage,method=method,error_group=group,
          attempts=nrow(q),lrt_flags=sum(q$lrt_p<.01,na.rm=TRUE),lower_flags=sum(q$lower_flag),
          declined=sum(q$declined),technically_withheld=sum(q$technical),reported=sum(q$reported))
      }
    }
    for(group in c('all','lambda=0.432',unique(p$Distribution),unique(p$stratum))) for(filter in c('before','after')) {
      q <- if(group=='all')p else if(group=='lambda=0.432')p[abs(p$GT_Lambda-.432)<1e-12,] else p[p$Distribution==group | p$stratum==group,]
      total <- nrow(q); reported <- sum(q$reported)
      declines <- if('declined' %in% names(q))sum(q$declined) else 0L
      technical <- if('technical' %in% names(q))sum(q$technical) else 0L
      if(filter=='after')q <- q[q$reported,]
      err <- q$zzDevAbs_Ov; missing <- !is.finite(err)
      summaries[[length(summaries)+1]] <- data.frame(stage=stage,group=group,method=method,filter=filter,
        attempts=total,reported=reported,declined=declines,technical_withheld=technical,scored=sum(!missing),missing=sum(missing),
        decline_or_withhold_rate=1-reported/total,mean_abs_z=mean(err,na.rm=TRUE),
        capped_mean=mean(q$zzDevAbsCutoff_Ov,na.rm=TRUE),
        all_attempt_capped=mean(ifelse(missing,5,q$zzDevAbsCutoff_Ov)),
        implausible=sum(err>5,na.rm=TRUE),implausible_or_missing=sum(missing|err>5),
        fit_seconds=sum(q$Runtime,na.rm=TRUE))
    }
  }
  saveRDS(all,file.path(w,'native-performance.rds'))
  category <- RIbench:::getBenchmarkResults(errorList=all,nameVec=methods,tableTCs=sets,cutoffZ=5,perfCombination='mean')
  saveRDS(category,file.path(w,'category-score.rds'))
  capture.output(print(category),file=file.path(w,'category-score.txt'))
  for(arm in c('armA','armB')) for(other in c('anchoredRI','anchorsOnly','refineR')) {
    x <- merge(scored[[arm]][,c('Index','zzDevAbs_Ov','reported','Distribution','fractionPathol','GT_Lambda')],
               scored[[other]][,c('Index','zzDevAbs_Ov','reported')],by='Index',suffixes=c('_arm','_other'))
    x$stratum <- paste(x$Distribution,ifelse(x$fractionPathol<=.2,'low <=20%','high >20%'))
    for(group in c('all','lambda=0.432',unique(x$stratum))) for(filter in c('before','after')) {
      q <- if(group=='all')x else if(group=='lambda=0.432') x[abs(x$GT_Lambda-.432)<1e-12,] else x[x$stratum==group,]
      attempted <- nrow(q)
      if(filter=='after')q <- q[q$reported_arm&q$reported_other,]
      delta <- q$zzDevAbs_Ov_arm-q$zzDevAbs_Ov_other; ok <- is.finite(delta)
      pairs[[length(pairs)+1]] <- data.frame(stage=stage,arm=arm,comparator=other,stratum=group,
        filter=filter,attempts=attempted,retained=nrow(q),complete=sum(ok),
        mean_difference=mean(delta[ok]),wins=sum(delta[ok]<0),ties=sum(delta[ok]==0),losses=sum(delta[ok]>0))
    }
  }
  # Preserve the complete native tables. Audit endpoints separately from the
  # primary mean-absolute-z > 5 definition.
  for(method in methods) {
    p <- scored[[method]]
    bad <- abs(p$zzDev_LRL)>5|abs(p$zzDev_URL)>5
    bad[p$Analyte=='CRP'] <- abs(p$zzDev_URL[p$Analyte=='CRP'])>5
    write.csv(p[which(bad),],file.path(w,paste0(method,'-extreme-errors.csv')),row.names=FALSE)
  }
  diagnostics[[stage]] <- do.call(rbind,lapply(c('armA','armB','anchoredRI'),function(method) {
    p <- scored[[method]]
    p[,c('Index','Analyte','method','transform','lambda','shift','components','at_bound',
          'anchor_seed','anchor_hash','routine_min','lrt_statistic','lrt_p','lower_flag',
          'declined','technical','reason','LRL','URL','zzDev_LRL','zzDev_URL','zzDevAbs_Ov',
          'Distribution','fractionPathol')]
  }))
}
summary <- do.call(rbind,summaries)
write.csv(summary,file.path(root,'summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,flags),file.path(root,'safeguard-by-error.csv'),row.names=FALSE)
write.csv(do.call(rbind,pairs),file.path(root,'paired-comparisons.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(root,'per-dataset-diagnostics.csv'),row.names=FALSE)
verdicts <- lapply(c('armA','armB'),function(arm) {
  select <- function(stage,method)summary[summary$stage==stage&summary$group=='all'&
    summary$method==method&summary$filter=='before',]
  small <- select('small',arm); large <- select('large',arm)
  data.frame(arm=arm,large_error_pass=large$all_attempt_capped<=select('large','refineR')$capped_mean,
    large_implausible_pass=large$implausible_or_missing<=14,
    small_error_pass=small$all_attempt_capped<=select('small','anchoredRI')$capped_mean)
})
verdicts <- do.call(rbind,verdicts)
write.csv(verdicts,file.path(root,'success-criteria.csv'),row.names=FALSE)
# Simulation comparisons retain frozen stress_metrics and paired 20-replicate MCSE.
sim <- list()
for(case in 1:45) for(replicate in 1:20) {
  key <- sprintf('s-%02d-%02d',case,replicate); d <- stress_data(case,replicate)
  r <- readRDS(file.path(root,'fits',paste0(key,'.rds')))
  old <- readRDS(file.path(root,'validation',paste0(key,'.rds')))
  r$anchoredRI <- list(fit=old$baseline,safeguard=old$safeguard,seconds=NA_real_)
  for(method in c('armA','armB','anchoredRI')) {
    f <- r[[method]]$fit; g <- r[[method]]$safeguard
    for(filter in c('before','after')) {
      limits <- if(is.null(g))c(NA_real_,NA_real_) else if(filter=='before')g$limits else g$reported
      row <- suppressWarnings(stress_row(d,method,limits,elapsed=r[[method]]$seconds))
      row$n_screened <- length(d$anchor); row$filter <- filter
      sim[[length(sim)+1]] <- row
    }
  }
}
sim <- do.call(rbind,sim)
read_comparator <- function(zip,member,method,label) {
  conn <- unz(zip,member,open='r'); on.exit(close(conn)); all <- read.csv(conn)
  rows <- all[all$method==method & all$replicate %in% 1:20,]; rows$method <- label
  stopifnot(nrow(rows)==900L,!anyDuplicated(rows[,c('case','replicate')]),
            all(rows$seed==8300000L+10000L*rows$case+rows$replicate))
  rows
}
refine <- read_comparator('assets/comparison-study/summary/comparison-results.zip','all-results.csv','refineR','refineR')
auto <- read_comparator('assets/anchors-only/anchors-only-results.zip','estimates.csv','anchors_auto','anchorsOnly')
legacy <- rbind(refine[,names(auto)],auto)
for(filter in c('before','after')) {
  legacy$filter <- filter
  sim <- rbind(sim,legacy[,names(sim)])
}
write.csv(sim,file.path(root,'simulation-results.csv'),row.names=FALSE)
regression <- list()
for(arm in c('armA','armB')) for(other in c('anchoredRI','anchorsOnly','refineR')) for(filter in c('before','after')) {
  x <- merge(sim[sim$method==arm&sim$filter==filter,],sim[sim$method==other&sim$filter==filter,],
    by=c('case','replicate'),suffixes=c('_arm','_other'))
  for(case in 1:45) {
    q <- x[x$case==case,]; delta <- q$width_error_arm-q$width_error_other; delta <- delta[is.finite(delta)]
    mcse <- if(length(delta)==20)sd(delta)/sqrt(20)else NA_real_
    regression[[length(regression)+1]] <- cbind(stress_design[case,c('family','overlap','affected')],
      data.frame(arm=arm,comparator=other,filter=filter,complete=length(delta),
        mean_difference=mean(delta),MCSE=mcse,wins=sum(delta<0),regression=mean(delta)>2*mcse))
  }
}
regression <- do.call(rbind,regression)
write.csv(regression,file.path(root,'simulation-paired.csv'),row.names=FALSE)

text <- c('# Task E results', '',
  'The baseline is anchoredRI at commit 55eee86b3f4111197b0043cc062ef0b781dc2057 with the frozen C2 adapter. The local refit used one routine group for RIbench and two groups for the simulation grid.', '',
  'All 5,220 same-machine input and equivalence checks passed, and the amended cross-machine review gate passed, before either experimental arm ran. Native RIbench scoring uses the same 2,880 small-N and 1,440 N=50,000 datasets and independent anchor seeds as C2.', '',
  '| Stage | Subset | Method | Safeguard | Scored | Mean absolute z | Capped mean | Reported / attempted |',
  '|---|---|---|---|---:|---:|---:|---:|')
for(i in which(summary$group %in% c('all','lambda=0.432'))) {
  r <- summary[i,]
  text <- c(text,sprintf('| %s | %s | %s | %s | %d | %.4f | %.4f | %d / %d |',
    r$stage,r$group,r$method,r$filter,r$scored,r$mean_abs_z,r$capped_mean,r$reported,r$attempts))
}
text <- c(text,'','## Prespecified success criteria','',
  'Success is judged before the safeguard. Missing fits receive a capped error of 5 and count toward the 1% implausible-result threshold. Each arm must meet all three criteria.','',
  '| Arm | N = 50,000 capped error at or below refineR | N = 50,000 implausible or missing (maximum 14/1,440) | Small-N capped error at or below local anchoredRI |',
  '|---|---|---|---|')
for(i in seq_len(nrow(verdicts))) {
  v <- verdicts[i,]; large <- summary[summary$method==v$arm & summary$stage=='large' & summary$group=='all' & summary$filter=='before',]
  pass <- function(x)if(isTRUE(x))'Pass' else if(isFALSE(x))'Fail' else 'Unavailable'
  text <- c(text,sprintf('| %s | %s | %s (%d/1,440) | %s |',v$arm,pass(v$large_error_pass),
    pass(v$large_implausible_pass),large$implausible_or_missing,pass(v$small_error_pass)))
}
text <- c(text,'',
  'After-safeguard errors describe the intervals still reported. Paired comparisons after the safeguard use only pairs retained by both methods. Declined estimates are never replaced by an anchors-only interval.', '',
  'The fitted Box-Cox likelihood retains the original estimator’s untruncated transformed-normal convention. Nonzero shapes have a restricted inverse domain. Invalid endpoints are withheld without clipping.', '',
  'Heavily mixed populations with many mildly abnormal patients, including ICU-only populations, can still produce wrong intervals that pass both checks.', '',
  'Simulation replicates 1–20 were used for development and regression checks. A paired increase exceeding twice its Monte Carlo standard error flags a regression.', '',
  'The frozen anchoredRI baseline was refitted in the same runtime as both arms. Anchors-only and refineR were retained from the archives and can have small cross-platform numerical differences. Those differences were not measured by refitting. Validation and fit wall times are saved separately. Input generation reuses the hash-verified archived data, so no new generation time is attributed to this study.')
ties <- read.csv(file.path(root,'baseline-near-tie-summary.csv'))
text <- c(text,'','## Baseline model selection','',
          '| Dataset group | Within 2 BIC points / attempted | Percent of all fits | With a finite runner-up | Percent of finite pairs within 2 | Median largest limit change (% of interval width) | Maximum (%) |',
          '|---|---:|---:|---:|---:|---:|---:|')
for(i in which(ties$group %in% c('all','small','50000','simulation'))) {
  r <- ties[i,]; text <- c(text,sprintf('| %s | %d / %d | %.2f | %d | %.2f | %.3f | %.3f |',r$group,r$within_2,
      r$attempts,100*r$share_all,r$finite_gaps,100*r$share_finite,100*r$median_fraction_of_width,100*r$max_fraction_of_width))
}
text <- c(text,'','The near-tie CSV reports both requested limits and their signed and absolute differences in measurement units for each fit. The width is the selected interval width. For CRP it is the zero-to-upper span. The per-method candidate and BIC-gap tables include all eligible models.','')
cross <- read.csv(file.path(root,'cross-machine-summary.csv'))
text <- c(text,sprintf('The frozen baseline selected a different transform or component count on %d/%d RIbench datasets (%.2f%%) across machines. %d datasets exceeded the 0.1%% limit tolerance. The largest relative endpoint change was %.3f%% of the archived width.',
   cross$model_changes,cross$datasets,100*cross$model_change_fraction,cross$outside,100*cross$max_fraction_of_width),
   'For Hb 100 the repeat on the C2 image selected the original identity model by 1.60 BIC points. For Hb 148 both new fits selected identity by about 37 points. C2 archived only the winners. [Hb candidate tables and original endpoints](https://github.com/amromeo/amromeo.github.io/blob/study/boxcox/anchoredRI/tools/boxcox/validation/HB-BIC.md)','',
   '## Simulation regressions','',
   '| Arm vs. local anchoredRI | Safeguard | Cells with 20 pairs | Cells flagged |',
   '|---|---|---:|---:|')
for(arm in c('armA','armB')) for(filter in c('before','after')) {
  q <- regression[regression$arm==arm & regression$comparator=='anchoredRI' & regression$filter==filter,]
  text <- c(text,sprintf('| %s | %s | %d / 45 | %d |',arm,filter,sum(q$complete==20L),sum(q$regression,na.rm=TRUE)))
}
text <- c(text,'','The paired-comparisons CSV gives wins, ties, losses and mean differences for each distribution and pathological fraction. The safeguard-by-error CSV gives both flag counts, technical failures and intervals reported.','')
repro <- file.path(root,'reproducibility','REPORT.md')
if(file.exists(repro)) text <- c(text,readLines(repro)) else text <- c(text,'The paired 156-dataset Arm A reproducibility results are pending.')
writeLines(text,file.path(root,'REPORT.md'))

writeLines('Native RIbench scoring and simulation summaries complete.',file.path(root,'SCORED.txt'))
