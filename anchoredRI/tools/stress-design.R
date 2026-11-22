# Frozen design for the outside-favorable-conditions study, 22 November 2026.
# No method is tuned on these results. Source this file from anchoredRI/.
stress_design <- expand.grid(family=c("normal","lognormal","shifted"),
  overlap=c("separated","moderate","heavy"),affected=c(.1,.3,.5,.7,.9),
  stringsAsFactors=FALSE)
stress_design$case <- seq_len(nrow(stress_design))
stress_design$delta <- c(separated=4,moderate=2,heavy=.75)[stress_design$overlap]
stress_design$disease_sd_z <- 1.25
stress_design$p_a <- stress_design$affected-.05
stress_design$p_b <- stress_design$affected+.05
stress_design$n_routine <- 5000L
stress_design$n_screened <- 40L
stress_design$refine_model <- ifelse(stress_design$family=="shifted","modBoxCox","BoxCox")
stress_design$disease_inside_true_ri <- with(stress_design,
  pnorm(qnorm(.975),delta,disease_sd_z)-pnorm(qnorm(.025),delta,disease_sd_z))
stress_forward <- function(z,family) switch(family,
  normal=30+3*z,lognormal=exp(log(20)+.35*z),shifted=10+exp(log(15)+.5*z))
stress_cdf <- function(x,family) switch(family,
  normal=pnorm(x,30,3),lognormal=plnorm(pmax(x,0),log(20),.35),
  shifted=plnorm(pmax(x-10,0),log(15),.5))
stress_data <- function(case,replicate,n_routine=5000L) {
  stopifnot(case %in% stress_design$case,replicate>=1L,n_routine>=500L,n_routine%%2==0L)
  s <- stress_design[case,]
  # Fresh routine patients per cell and replicate. Anchors are reused across
  # cells of one family/replicate to hold their quality fixed across stresses.
  seed <- 8300000L+10000L*case+replicate
  set.seed(seed)
  setting <- rep(c("service_a","service_b"),each=n_routine/2L)
  affected <- runif(n_routine)<rep(c(s$p_a,s$p_b),each=n_routine/2L)
  z <- rnorm(n_routine)
  dz <- rnorm(n_routine,s$delta,s$disease_sd_z)
  z[affected] <- dz[affected]
  value <- round(stress_forward(z,s$family),2)
  set.seed(9300000L+10000L*match(s$family,c("normal","lognormal","shifted"))+replicate)
  screened <- round(stress_forward(rnorm(400L),s$family),2)
  list(value=value,setting=setting,anchor=head(screened,40L),screened=screened,
    truth=stress_forward(qnorm(c(.025,.975)),s$family),seed=seed,
    actual_affected=mean(affected),family=s$family,case=case,replicate=replicate,
    actual_a=mean(affected[setting=="service_a"]),actual_b=mean(affected[setting=="service_b"]))
}
stress_valid <- function(x) is.numeric(x)&&length(x)==2L&&all(is.finite(x))&&x[1]<x[2]
stress_metrics <- function(limits,truth,family) {
  ok <- stress_valid(limits)
  if(!ok)limits <- c(NA_real_,NA_real_)
  width_error <- 100*max(abs(limits-truth))/diff(truth)
  tails <- if(ok)stress_cdf(limits,family)else c(NA_real_,NA_real_)
  data.frame(lower=limits[1],upper=limits[2],true_lower=truth[1],true_upper=truth[2],
    lower_error=limits[1]-truth[1],upper_error=limits[2]-truth[2],
    width_error=width_error,mpe=100*mean(abs(limits-truth)/abs(truth)),
    below_reference=tails[1],above_reference=1-tails[2],reference_coverage=diff(tails),
    large_error=if(ok)width_error>10 else NA,status=if(ok)"ok"else"failed",row.names=NULL)
}
stress_row <- function(d,method,limits,model="",shift=NA_real_,components=NA_integer_,
                       at_bound=NA,error="",elapsed=NA_real_,ci=matrix(NA_real_,2,2),
                       boot_requested=0L,boot_successful=0L,boot_failed=0L) {
  s <- stress_design[d$case,]
  m <- stress_metrics(limits,d$truth,d$family)
  cbind(data.frame(case=d$case,family=s$family,overlap=s$overlap,affected=s$affected,
    replicate=d$replicate,method=method,seed=d$seed,
    n_routine=if(grepl("^direct",method))0L else length(d$value),
    n_screened=if(method=="refineR")0L else if(method=="anchored40")40L else as.integer(sub("direct","",method)),
    actual_affected=d$actual_affected,actual_a=d$actual_a,actual_b=d$actual_b,
    model=model,shift=shift,components=components,at_bound=at_bound,elapsed=elapsed,
    lower_ci_low=ci[1,1],lower_ci_high=ci[1,2],upper_ci_low=ci[2,1],upper_ci_high=ci[2,2],
    lower_ci_covers=if(all(is.finite(ci[1,])))ci[1,1]<=d$truth[1]&&ci[1,2]>=d$truth[1] else NA,
    upper_ci_covers=if(all(is.finite(ci[2,])))ci[2,1]<=d$truth[2]&&ci[2,2]>=d$truth[2] else NA,
    bootstrap_requested=boot_requested,bootstrap_successful=boot_successful,
    bootstrap_failed=boot_failed,error=error,stringsAsFactors=FALSE),m)
}
stress_hashes <- function() tools::md5sum(c(sort(list.files("R",full.names=TRUE)),
  "tools/stress-design.R","tools/run-stress-grid.R"))
