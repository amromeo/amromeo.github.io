# Run from anchoredRI: Rscript tools/validate-offset.R [workers=2] [reps=20]
# Fixed before generating the fresh datasets. This compares old/new anchoredRI,
# not repeated performance against refineR; no bootstrap-CI coverage claim.
for(path in list.files("R",pattern="[.]R$",full.names=TRUE))source(path)
args <- commandArgs(trailingOnly=TRUE)
workers <- if(length(args)>=1L)as.integer(args[1])else 2L
reps <- if(length(args)>=2L)as.integer(args[2])else 20L
stopifnot(workers>=1L,reps>=2L)
dir.create("assets/offset",showWarnings=FALSE)
dir.create("tools/offset-validation-cache",showWarnings=FALSE)
source_files <- sort(list.files("R",full.names=TRUE));hash <- tools::md5sum(source_files)
scenarios <- data.frame(
  id=c("normal","overlap","lognormal","shift10","shift3","shift20","skew03","skew075","anchor20","anchor80","shift_overlap"),
  title=c("Normal","Normal, overlapping disease","Lognormal, no shift","Shift 10, SDlog 0.50","Shift 3","Shift 20","Shift 10, SDlog 0.30","Shift 10, SDlog 0.75","Shift 10, 20 screened","Shift 10, 80 screened","Shift 10, overlapping disease"),
  family=c("normal","normal",rep("lognormal",9)),
  shift=c(0,0,0,10,3,20,10,10,10,10,10),
  median=c(30,30,20,rep(15,8)),
  sd=c(3,3,.35,.5,.5,.5,.3,.75,.5,.5,.5),
  disease_median=c(42,36,60,rep(65,7),30),
  disease_sd=c(4,4.5,.3,rep(.25,7),.5),
  p_a=c(.05,.2,rep(.05,8),.2),p_b=c(.15,.3,rep(.15,8),.3),
  n_anchor=c(rep(40L,8),20L,80L,40L))
write.csv(scenarios,"assets/offset/validation-design.csv",row.names=FALSE)
make_fresh <- function(s,r) {
  set.seed(7610000L+10000L*match(s$id,scenarios$id)+r)
  qref <- function(p)if(s$family=="normal")qnorm(p,s$median,s$sd)else s$shift+qlnorm(p,log(s$median),s$sd)
  cdf <- function(x)if(s$family=="normal")pnorm(x,s$median,s$sd)else plnorm(pmax(0,x-s$shift),log(s$median),s$sd)
  setting <- rep(c("a","b"),each=2500L)
  affected <- runif(5000L)<rep(c(s$p_a,s$p_b),each=2500L)
  value <- qref(runif(5000L))
  disease <- if(s$family=="normal")rnorm(5000L,s$disease_median,s$disease_sd)else
    s$shift+rlnorm(5000L,log(s$disease_median),s$disease_sd)
  value[affected] <- disease[affected]
  list(value=round(value,2),setting=setting,anchor=round(qref(runif(s$n_anchor)),2),
       truth=qref(c(.025,.975)),cdf=cdf)
}
grid <- expand.grid(scenario=seq_len(nrow(scenarios)),replicate=seq_len(reps))
run_one <- function(i) {
  s <- scenarios[grid$scenario[i],];r <- grid$replicate[i]
  cache <- sprintf("tools/offset-validation-cache/%s-%03d.rds",s$id,r)
  if(file.exists(cache)) {
    old <- readRDS(cache)
    if(!identical(old$hash,hash))stop("Cache uses a different estimator.")
    return(old$result)
  }
  d <- make_fresh(s,r)
  result <- do.call(rbind,lapply(c("original","auto"),function(method) {
    transform <- if(method=="auto")"auto"else if(s$family=="normal")"identity"else"log"
    start <- proc.time()[["elapsed"]]
    f <- tryCatch(fit_anchored_ri(d$value,d$setting,d$anchor,transform),error=identity)
    elapsed <- proc.time()[["elapsed"]]-start
    ok <- !inherits(f,"error") && all(is.finite(f$limits))
    lim <- if(ok)f$limits else c(NA_real_,NA_real_)
    data.frame(scenario=s$id,replicate=r,method=method,status=if(ok)"ok"else"failed",
      lower=lim[1],upper=lim[2],true_lower=d$truth[1],true_upper=d$truth[2],
      lower_error=lim[1]-d$truth[1],upper_error=lim[2]-d$truth[2],
      mpe=100*mean(abs(lim-d$truth)/abs(d$truth)),worst_error=max(abs(lim-d$truth)),
      reference_coverage=if(ok)diff(d$cdf(lim))else NA_real_,
      transform=if(ok)f$transform else NA_character_,shift=if(ok)f$shift else NA_real_,
      components=if(ok)f$components else NA_integer_,at_bound=if(ok)f$at_bound else NA,
      elapsed=elapsed,error=if(ok)""else conditionMessage(f))
  }))
  saveRDS(list(hash=hash,result=result),cache)
  message("Completed ",s$id," ",r,"/",reps)
  result
}
out <- parallel::mclapply(seq_len(nrow(grid)),run_one,mc.cores=workers,mc.set.seed=FALSE)
stopifnot(all(vapply(out,is.data.frame,logical(1))),identical(hash,tools::md5sum(source_files)))
out <- do.call(rbind,out)
write.csv(out,"assets/offset/validation-fits.csv",row.names=FALSE,na="")
summary <- do.call(rbind,lapply(scenarios$id,function(id) {
  x <- out[out$scenario==id & out$method=="original",];y <- out[out$scenario==id & out$method=="auto",]
  stopifnot(identical(x$replicate,y$replicate))
  good <- x$status=="ok" & y$status=="ok"
  delta <- y$mpe[good]-x$mpe[good]
  data.frame(scenario=id,title=scenarios$title[match(id,scenarios$id)],attempted=nrow(x),paired_successful=sum(good),
    original_failed=sum(x$status!="ok"),auto_failed=sum(y$status!="ok"),
    original_mpe=mean(x$mpe[good]),auto_mpe=mean(y$mpe[good]),
    paired_change_mpe=mean(delta),paired_change_mcse=sd(delta)/sqrt(length(delta)),
    auto_better=sum(delta < -1e-6),ties=sum(abs(delta)<=1e-6),auto_worse=sum(delta>1e-6),
    original_worst_error=mean(x$worst_error[good]),auto_worst_error=mean(y$worst_error[good]),
    original_reference_coverage=mean(x$reference_coverage[good]),auto_reference_coverage=mean(y$reference_coverage[good]),
    auto_identity=sum(y$transform=="identity",na.rm=TRUE),auto_log=sum(y$transform=="log",na.rm=TRUE),
    auto_shifted_log=sum(y$transform=="shifted_log",na.rm=TRUE),auto_boundary=sum(y$at_bound,na.rm=TRUE))
}))
write.csv(summary,"assets/offset/validation-summary.csv",row.names=FALSE,na="")
print(summary,row.names=FALSE)
