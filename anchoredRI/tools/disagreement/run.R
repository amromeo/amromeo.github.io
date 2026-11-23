source('tools/stress-design.R')
source('tools/anchors-only/estimators.R')
# Same grid likelihood, vectorized to avoid 60 million R function calls.
ao_fast <- ao_auto
body(ao_fast) <- parse(text=sub('vals <- vapply(grid, objective, numeric(1))',
 'yg <- log(outer(a, grid, "-")); mg <- colMeans(yg); vg <- colMeans((yg-rep(mg,each=length(a)))^2); vals <- -length(a)/2*(log(2*pi)+1+log(vg))-colSums(yg)',
 paste(deparse(body(ao_auto)),collapse='\n'),fixed=TRUE))[[1]]
check_interval <- function(limits,ci) {
 available <- stress_valid(limits)&&all(is.finite(ci))
 flag <- if(available)any(limits<ci[,1]|limits>ci[,2])else NA
 list(flag=flag,declined=!available||flag,reported=if(available&&!flag)limits else c(NA_real_,NA_real_))
}
lrt <- function(a,limits,model,shift) {
 trans <- function(x)switch(model,identity=x,log=log(x),shifted_log=log(x-shift),stop('Unknown transformation'))
 y<-suppressWarnings(trans(a));z<-suppressWarnings(trans(limits))
 if(any(!is.finite(c(y,z)))||!stress_valid(z))return(c(statistic=NA,p=NA))
 mu<-mean(z);sigma<-diff(z)/(2*qnorm(.975))
 own<-ao_normal(y);stat<-max(0,2*(own$loglik-sum(dnorm(y,mu,sigma,log=TRUE))))
 c(statistic=stat,p=pchisq(stat,2,lower.tail=FALSE))
}
args<-commandArgs(TRUE)
if(length(args)&&args[1]=='test') {
 stopifnot(!check_interval(c(1,2),matrix(c(1,2,1,2),2))$declined)
 stopifnot(all(is.na(check_interval(c(0,2),matrix(c(1,2,1,2),2))$reported)))
 stopifnot(check_interval(c(1,2),matrix(NA_real_,2,2))$declined)
 for(f in 1:3)for(r in 1:20){a<-stress_data(f,r)$anchor;set.seed(43+r);a<-sample(a,40,TRUE);x<-ao_auto(a);y<-ao_fast(a);stopifnot(x$transform==y$transform,max(abs(x$limits-y$limits))<1e-5,abs(x$bic-y$bic)<1e-6)}
 a<-1:40;fit<-ao_normal(a);stopifnot(lrt(a,fit$mu+qnorm(c(.025,.975))*fit$mle_sd,'identity',0)['p']> .99999)
 cat('Boundary, decline, LRT and 60 full-selection equivalence checks passed\n')
} else if(length(args)) {
 id<-as.integer(args[1]);f<-(id-1)%/%20+1;r<-(id-1)%%20+1
 a<-stress_data(f,r)$anchor;set.seed(19400000+10000*f+r)
 ix<-matrix(sample.int(40,40000,replace=TRUE),40,1000)
 rows<-lapply(1:1000,function(b){fit<-tryCatch(ao_fast(a[ix[,b]]),error=function(e)list(limits=c(NA_real_,NA_real_),transform='',shift=NA_real_));data.frame(bootstrap=b,lower=fit$limits[1],upper=fit$limits[2],model=fit$transform,shift=fit$shift)})
 boot<-do.call(rbind,rows);ok<-apply(boot[,c('lower','upper')],1,stress_valid)
 ci<-if(all(ok))t(sapply(boot[,c('lower','upper')],quantile,c(.025,.975),type=7))else matrix(NA_real_,2,2)
 saveRDS(list(anchor=a,bootstrap=boot,ci=ci,failed=sum(!ok),auto=ao_auto(a),seed=19400000+10000*f+r),sprintf('assets/disagreement/bootstrap-%02d.rds',id))
 cat(id,'completed\n')
} else {
 prior<-read.csv(unz('assets/comparison-study/summary/comparison-results.zip','all-results.csv'));prior<-prior[prior$method=='anchored40',]
 stopifnot(nrow(prior)==900,!anyDuplicated(prior[c('case','replicate')]))
 rows<-list();bootrows<-list()
 for(i in seq_len(nrow(prior))){p<-prior[i,];d<-stress_data(p$case,p$replicate);f<-match(d$family,c('normal','lognormal','shifted'));id<-(f-1)*20+p$replicate
 b<-readRDS(sprintf('assets/disagreement/bootstrap-%02d.rds',id));stopifnot(identical(d$anchor,b$anchor),d$seed==p$seed,max(abs(d$truth-c(p$true_lower,p$true_upper)))<1e-10)
 limits<-c(p$lower,p$upper);v<-check_interval(limits,b$ci);l<-lrt(d$anchor,limits,p$model,p$shift)
 metrics<-stress_metrics(v$reported,d$truth,d$family)
 rows[[i]]<-data.frame(case=p$case,family=p$family,overlap=p$overlap,affected=p$affected,replicate=p$replicate,original_error=p$width_error,
 diagnostic_anchored_lower=p$lower,diagnostic_anchored_upper=p$upper,diagnostic_auto_lower=b$auto$limits[1],diagnostic_auto_upper=b$auto$limits[2],
 lower_ci_low=b$ci[1,1],lower_ci_high=b$ci[1,2],upper_ci_low=b$ci[2,1],upper_ci_high=b$ci[2,2],bootstrap_failed=b$failed,
 primary_flag=v$flag,primary_declined=v$declined,reported_lower=v$reported[1],reported_upper=v$reported[2],reported_error=metrics$width_error,
 recommendation=if(v$declined)'At least 120 screened patients or verify a suitable published interval with 20'else'',
 lrt_statistic=l['statistic'],lrt_p=l['p'],secondary_flag=if(is.finite(l['p']))l['p']<.01 else NA,secondary_declined=!is.finite(l['p'])||l['p']<.01)
 }
 write.csv(do.call(rbind,rows),'assets/disagreement/datasets.csv',row.names=FALSE,na='')
 for(id in 1:60){b<-readRDS(sprintf('assets/disagreement/bootstrap-%02d.rds',id));bootrows[[id]]<-cbind(anchor_sample=id,seed=b$seed,b$bootstrap)}
 write.csv(do.call(rbind,bootrows),'assets/disagreement/bootstrap-estimates.csv',row.names=FALSE,na='')
 writeLines(c(capture.output(sessionInfo()),paste('Input archive MD5:',tools::md5sum('assets/comparison-study/summary/comparison-results.zip'))),'assets/disagreement/versions.txt')
 cat('Scored 900 datasets\n')
}
