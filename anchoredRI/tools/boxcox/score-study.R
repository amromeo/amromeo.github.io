library(RIbench)
source('tools/boxcox/models.R')
source('tools/stress-design.R')
source('tools/anchors-only/estimators.R')
cfg <- readRDS(commandArgs(TRUE)[1]); root <- cfg$output
stopifnot(file.exists(file.path(root,'FITTED.txt')),!file.exists(file.path(root,'HALT.txt')))
t <- loadTestsetDefinition(); t <- t[t$N<=50000,]
methods <- c('armA','armB','anchoredRI','anchorsOnly','refineR')
diagnostics <- list(); summaries <- list(); pairs <- list()
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
    scored[[method]] <- p
    write.csv(p,file.path(w,paste0(method,'-performance.csv')),row.names=FALSE)
    for(group in c('all','lambda=0.432',unique(p$Distribution))) for(filter in c('before','after')) {
      q <- if(group=='all')p else if(group=='lambda=0.432')p[abs(p$GT_Lambda-.432)<1e-12,] else p[p$Distribution==group,]
      total <- nrow(q); reported <- sum(q$reported)
      if(filter=='after')q <- q[q$reported,]
      err <- q$zzDevAbs_Ov; missing <- !is.finite(err)
      summaries[[length(summaries)+1]] <- data.frame(stage=stage,group=group,method=method,filter=filter,
        attempts=total,reported=reported,scored=sum(!missing),missing=sum(missing),
        decline_or_withhold_rate=1-reported/total,mean_abs_z=mean(err,na.rm=TRUE),
        capped_mean=mean(q$zzDevAbsCutoff_Ov,na.rm=TRUE),
        all_attempt_capped=mean(ifelse(missing,5,pmin(err,5))),
        implausible=sum(err>5,na.rm=TRUE),implausible_or_missing=sum(missing|err>5),
        fit_seconds=sum(q$Runtime,na.rm=TRUE))
    }
  }
  saveRDS(all,file.path(w,'native-performance.rds'))
  category <- RIbench:::getBenchmarkResults(errorList=all,nameVec=methods,tableTCs=sets,cutoffZ=5,perfCombination='mean')
  saveRDS(category,file.path(w,'category-score.rds'))
  capture.output(print(category),file=file.path(w,'category-score.txt'))
  for(arm in c('armA','armB')) for(other in c('anchoredRI','anchorsOnly','refineR')) {
    x <- merge(scored[[arm]][,c('Index','zzDevAbs_Ov','reported','Distribution','fractionPathol')],
               scored[[other]][,c('Index','zzDevAbs_Ov','reported')],by='Index',suffixes=c('_arm','_other'))
    x$stratum <- paste(x$Distribution,ifelse(x$fractionPathol<=.2,'low <=20%','high >20%'))
    for(group in c('all',unique(x$stratum))) for(filter in c('before','after')) {
      q <- if(group=='all')x else x[x$stratum==group,]; attempted <- nrow(q)
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
write.csv(do.call(rbind,verdicts),file.path(root,'success-criteria.csv'),row.names=FALSE)
# Simulation comparisons retain frozen stress_metrics and paired 20-replicate MCSE.
sim <- list()
for(case in 1:45) for(replicate in 1:20) {
  key <- sprintf('s-%02d-%02d',case,replicate); d <- stress_data(case,replicate)
  r <- readRDS(file.path(root,'fits',paste0(key,'.rds')))
  old <- readRDS(file.path(root,'validation',paste0(key,'.rds')))
  r$anchoredRI <- list(fit=old$baseline,safeguard=old$safeguard,seconds=NA_real_)
  auto <- ao_auto(d$anchor)
  r$anchorsOnly <- list(fit=auto,safeguard=list(limits=auto$limits,reported=auto$limits,declined=FALSE,technical=FALSE),seconds=NA_real_)
  for(method in names(r)[names(r)%in%c('armA','armB','anchoredRI','anchorsOnly')]) {
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
legacy <- read.csv(cfg$simulation_archive)
legacy <- legacy[legacy$method=='refineR'&legacy$replicate<=20,]
stopifnot(nrow(legacy)==900,!anyDuplicated(paste(legacy$case,legacy$replicate)))
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
write.csv(do.call(rbind,regression),file.path(root,'simulation-paired.csv'),row.names=FALSE)
writeLines('Native RIbench scoring and simulation summaries complete. See success-criteria.csv and paired tables.',
           file.path(root,'SCORED.txt'))

text <- c('# Task E results', '',
  'Baseline commit 55eee86b3f4111197b0043cc062ef0b781dc2057. RIbench uses the frozen C2 one-group adapter. Simulation uses the original two-setting estimator.', '',
  'All 5,220 same-machine input and equivalence checks passed, and the amended cross-machine review gate passed, before either experimental arm ran. Native RIbench scoring uses the same 2,880 small-N and 1,440 N=50,000 datasets and independent anchor seeds as C2.', '',
  '| Stage | Subset | Method | Safeguard | Scored | Mean absolute z | Capped mean | Reported / attempted |',
  '|---|---|---|---|---:|---:|---:|---:|')
for(i in which(summary$group %in% c('all','lambda=0.432'))) {
  r <- summary[i,]
  text <- c(text,sprintf('| %s | %s | %s | %s | %d | %.4f | %.4f | %d / %d |',
    r$stage,r$group,r$method,r$filter,r$scored,r$mean_abs_z,r$capped_mean,r$reported,r$attempts))
}
text <- c(text,'','Success is judged before the safeguard. Missing fits receive a capped error of 5 in the all-attempt sensitivity analysis and count against the 1% implausible-result threshold. See success-criteria.csv for the separate arm verdicts.', '',
  'After-safeguard errors describe the intervals still reported. Paired comparisons after the safeguard use only pairs retained by both methods. Declined estimates are never replaced by an anchors-only interval.', '',
  'The fitted Box-Cox likelihood retains the original estimator’s untruncated transformed-normal convention. Nonzero shapes have a restricted inverse domain. Invalid endpoints are withheld without clipping.', '',
  'Heavily mixed populations with many mildly abnormal patients, including ICU-only populations, can still produce wrong intervals that pass both checks.', '',
  'Simulation replicates 1–20 are development and regression checks, not independent confirmation. A paired increase exceeding twice its Monte Carlo standard error flags a regression. Unflagged differences do not establish equivalence.', '',
  'The frozen anchoredRI baseline was refitted in the same runtime as both arms. Anchors-only and refineR were retained from the archives and can have small cross-platform numerical differences. Those differences were not measured by refitting. Validation and fit wall times are saved separately. Input generation reuses the hash-verified archived data, so no new generation time is attributed to this study.')
writeLines(text,file.path(root,'REPORT.md'))
