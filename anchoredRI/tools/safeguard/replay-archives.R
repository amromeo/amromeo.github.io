# Run from anchoredRI/ with the proposed package installed. This replays saved fits.
# Rscript tools/safeguard/replay-archives.R config.rds
library(RIbench)
.ri_checks <- getFromNamespace(".ri_checks", "anchoredRI")
.ri_flags <- getFromNamespace(".ri_flags", "anchoredRI")
.ri_reporting_decision <- getFromNamespace(".ri_reporting_decision", "anchoredRI")
source('tools/stress-design.R')
cfg <- readRDS(commandArgs(TRUE)[1])
dir.create(cfg$output,recursive=TRUE,showWarnings=FALSE)
all_rows <- list(); totals <- list()
b2 <- read.csv(gzfile(file.path(cfg$b2,'results','confirmation-datasets.csv.gz')))
b2 <- b2[b2$method=='lrt',]
stopifnot(nrow(b2)==900,all(b2$replicate %in% 21:40),
          !anyDuplicated(b2[,c('case','replicate')]))
for(i in seq_len(nrow(b2))) {
  r <- b2[i,];d <- stress_data(r$case,r$replicate)
  f <- list(limits=c(r$diagnostic_lower,r$diagnostic_upper),transform=r$transform,shift=r$offset)
  check <- .ri_checks(f,d$value,d$anchor)
  stopifnot(identical(check$lrt_flag,r$lrt_flag),!check$technical)
  decision <- .ri_reporting_decision(check,TRUE,length(d$value))
  all_rows[[length(all_rows)+1]] <- data.frame(stage='B2',key=paste(r$case,r$replicate,sep='-'),
    lrt_recomputed=TRUE,lrt_flag=check$lrt_flag,lower_flag=check$lower_flag,
    declined=decision$declined,lrt_p=check$lrt_p)
}
z <- do.call(rbind,all_rows)
stopifnot(sum(z$lrt_flag)==91,sum(z$lrt_flag[b2$error>10])==86,
          sum(z$lrt_flag[b2$error<5])==5,sum(z$lower_flag)==0,sum(z$declined)==91)
totals[[1]] <- data.frame(stage='B2',n=900,lrt=sum(z$lrt_flag),lower=sum(z$lower_flag),
                         either=sum(z$declined),both=sum(z$lrt_flag&z$lower_flag))
t <- loadTestsetDefinition()
for(stage in c('small','large')) {
  w <- if(stage=='small')cfg$small else cfg$large
  archive <- read.csv(if(stage=='small')cfg$small_diagnostics else cfg$large_diagnostics)
  if(stage=='large')archive <- archive[archive$method=='anchoredRI',]
  archive$lower_flag <- as.logical(toupper(as.character(archive$lower_flag)))
  custom <- do.call(rbind,lapply(list.files(w,pattern='^custom-[0-9]+[.]csv$',full.names=TRUE),read.csv))
  custom <- custom[custom$algo=='anchoredRI',]
  stopifnot(!anyDuplicated(archive$Index),!anyDuplicated(custom$Index),
            setequal(archive$Index,custom$Index))
  rows <- list()
  for(i in seq_len(nrow(archive))) {
    r <- archive[i,];s <- t[t$Index==r$Index,];c <- custom[custom$Index==r$Index,]
    set.seed(17000000L+s$Index)
    a <- RIbench:::invBoxCox(rnorm(40,s$nonp_mu,s$nonp_sigma),s$nonp_lambda)+s$nonp_shift
    a <- round(a[is.finite(a)&a>=0],s$decimals)
    if(s$Analyte!='CRP') {
      check <- .ri_checks(list(limits=c(r$LRL,r$URL),transform=c$model,shift=c$shift),
                          r$routine_min,a)
      stopifnot(identical(check$lrt_flag,unname(c$lrt_p<.01)),!check$technical)
    } else {
      # CRP artifacts saved only the upper endpoint, not both component
      # parameters. Replay its saved LRT diagnostic, never infer parameters
      # from a single limit. Public-entry tests separately recompute its LRT.
      check <- .ri_flags(c(upper=r$URL),c$lrt_statistic,r$routine_min)
      stopifnot(check$lower_check_skipped,!check$technical,
                isTRUE(all.equal(check$lrt_p,c$lrt_p,tolerance=1e-12)),
                identical(check$lrt_flag,unname(c$lrt_p<.01)))
    }
    stopifnot(identical(check$lower_flag,unname(r$lower_flag)))
    decision <- .ri_reporting_decision(check,TRUE,s$N)
    expected <- c$lrt_p<.01 || r$lower_flag
    stopifnot(identical(decision$declined,expected))
    rows[[i]] <- data.frame(stage=stage,key=r$Index,lrt_recomputed=s$Analyte!='CRP',
      lrt_flag=check$lrt_flag,lower_flag=check$lower_flag,declined=decision$declined,lrt_p=check$lrt_p)
  }
  z <- do.call(rbind,rows)
  if(stage=='small') {
    stopifnot(nrow(z)==2880,sum(z$lrt_flag)==533,sum(z$lower_flag)==9,sum(z$declined)==533,
      sum(z$declined[archive$zzDevAbs_Ov>5])==23,
      sum(z$declined[archive$zzDevAbs_Ov<.25])==224,
      sum(z$declined[archive$endpoint_extreme])==32)
  } else stopifnot(nrow(z)==1440,sum(z$lower_flag)==21,sum(z$lrt_flag)==402,
                   sum(z$declined)==402,sum(z$lrt_flag&z$lower_flag)==21)
  totals[[length(totals)+1]] <- data.frame(stage=stage,n=nrow(z),lrt=sum(z$lrt_flag),
    lower=sum(z$lower_flag),either=sum(z$declined),both=sum(z$lrt_flag&z$lower_flag))
  all_rows <- c(all_rows,rows)
}
write.csv(do.call(rbind,all_rows),file.path(cfg$output,'archived-decisions.csv'),row.names=FALSE)
write.csv(do.call(rbind,totals),file.path(cfg$output,'archived-counts.csv'),row.names=FALSE)
cat('All archived decisions reproduced exactly\n');print(do.call(rbind,totals))
