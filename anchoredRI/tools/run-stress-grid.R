# Rscript tools/run-stress-grid.R MODE CASE [workers=4] [reps=20]
# MODE: repeat (20 independent point fits), bootstrap (replicate 1, B=200),
# or both. CASE is 1..45, or 0 to run all cases sequentially.
source("tools/stress-design.R")
for(p in sort(list.files("R",pattern="[.]R$",full.names=TRUE)))source(p)
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args))args[1]else"both"
case <- if(length(args)>1L)as.integer(args[2])else 0L
workers <- if(length(args)>2L)as.integer(args[3])else 4L
reps <- if(length(args)>3L)as.integer(args[4])else 20L
stopifnot(mode%in%c("repeat","bootstrap","both"),case%in%0:45,workers>=1L,reps>=2L)
hash <- stress_hashes()
dir.create("assets/stress-grid/cases",recursive=TRUE,showWarnings=FALSE)
dir.create("tools/stress-grid-cache",recursive=TRUE,showWarnings=FALSE)
write.csv(stress_design,"assets/stress-grid/design.csv",row.names=FALSE)
write.csv(data.frame(path=names(hash),md5=unname(hash)),"assets/stress-grid/source-hashes.csv",row.names=FALSE)
save_case <- function(x,path) {
  stopifnot(identical(hash,stress_hashes()))
  write.csv(x,path,row.names=FALSE,na="")
}
point_one <- function(case,r) {
  path <- sprintf("tools/stress-grid-cache/point-%02d-%02d.rds",case,r)
  if(file.exists(path)) {
    old<-readRDS(path);stopifnot(identical(hash,old$hash));return(old$result)
  }
  d <- stress_data(case,r)
  begin <- proc.time()[["elapsed"]]
  f <- tryCatch(fit_anchored_ri(d$value,d$setting,d$anchor,transform="auto"),error=identity)
  elapsed <- proc.time()[["elapsed"]]-begin
  ok <- !inherits(f,"error")&&stress_valid(f$limits)
  out <- list(stress_row(d,"anchored40",if(ok)f$limits else c(NA_real_,NA_real_),
    model=if(ok)f$transform else "",shift=if(ok)f$shift else NA_real_,
    components=if(ok)f$components else NA_integer_,at_bound=if(ok)f$at_bound else NA,
    error=if(inherits(f,"error"))conditionMessage(f)else if(!ok)"Invalid limits"else"",elapsed=elapsed))
  for(n in c(40L,120L,400L))out[[length(out)+1L]] <- stress_row(d,paste0("direct",n),
    quantile(head(d$screened,n),c(.025,.975),type=7,names=FALSE),model="Empirical quantiles (type 7)")
  out<-do.call(rbind,out);saveRDS(list(hash=hash,result=out),path)
  message(sprintf("Point case %02d replicate %02d complete",case,r))
  out
}
bootstrap_one <- function(case) {
  stopifnot(requireNamespace("refineR",quietly=TRUE),as.character(utils::packageVersion("refineR"))=="2.0.0")
  d<-stress_data(case,1L);s<-stress_design[case,];B<-200L;out<-list()
  for(method in c("refineR","anchored40")) {
    path<-sprintf("tools/stress-grid-cache/boot-%02d-%s.rds",case,method)
    if(file.exists(path)) {
      old<-readRDS(path);stopifnot(identical(hash,old$hash));out[[method]]<-old$result;next
    }
    message("Bootstrap case ",case," ",method," started")
    begin<-proc.time()[["elapsed"]]
    f<-tryCatch(if(method=="refineR")refineR::findRI(d$value,model=s$refine_model,
      NBootstrap=B,NCores=workers,seed=1200000L+case)else
      bootstrap_anchored_ri(d$value,d$setting,d$anchor,transform="auto",n_boot=B,
        seed=2200000L+case,workers=workers),error=identity)
    elapsed<-proc.time()[["elapsed"]]-begin
    if(inherits(f,"error")) {
      row<-stress_row(d,method,c(NA_real_,NA_real_),error=conditionMessage(f),elapsed=elapsed,
        boot_requested=B,boot_failed=B)
    } else if(method=="refineR") {
      ri<-tryCatch(refineR::getRI(f,pointEst="medianBS",CIprop=.95,UCMargins=FALSE),error=identity)
      successful<-sum(is.finite(f$MuBS)&is.finite(f$SigmaBS)&is.finite(f$LambdaBS)&is.finite(f$ShiftBS))
      good<-!inherits(ri,"error")&&stress_valid(ri$PointEst)&&successful>0L
      row<-stress_row(d,method,if(good)ri$PointEst else c(NA_real_,NA_real_),
        model=s$refine_model,ci=if(good)cbind(ri$CILow,ri$CIHigh)else matrix(NA_real_,2,2),
        boot_requested=B,boot_successful=successful,boot_failed=B-successful,elapsed=elapsed,
        error=if(inherits(ri,"error"))conditionMessage(ri)else if(!good)"Invalid bootstrap estimate"else"")
      if(good) {
        grDevices::svg(sprintf("assets/stress-grid/cases/fit-%02d.svg",case),width=10,height=5.5)
        plot(f,pointEst="medianBS",showPathol=FALSE,showMargin=FALSE,
          title=sprintf("%s | %s overlap | %.0f%% affected",s$family,s$overlap,100*s$affected),xlab="Synthetic units")
        grDevices::dev.off()
      }
      save_case(data.frame(lambda=f$LambdaBS,mu=f$MuBS,sigma=f$SigmaBS,shift=f$ShiftBS,
        proportion=f$PBS,cost=f$CostBS),sprintf("assets/stress-grid/cases/refine-bootstrap-%02d.csv",case))
    } else {
      row<-stress_row(d,method,f$limits,model=f$fit$transform,shift=f$fit$shift,
        components=f$fit$components,at_bound=f$fit$at_bound,ci=f$ci,
        boot_requested=B,boot_successful=f$successful,boot_failed=f$failed,elapsed=elapsed)
      save_case(f$bootstrap,sprintf("assets/stress-grid/cases/anchor-bootstrap-%02d.csv",case))
    }
    saveRDS(list(hash=hash,result=row,fit=f),path);out[[method]]<-row
    message("Bootstrap case ",case," ",method," complete in ",round(elapsed)," seconds")
  }
  for(n in c(40L,120L,400L)) {
    a<-head(d$screened,n);set.seed(3200000L+case*1000L+n)
    bs<-replicate(B,quantile(sample(a,n,replace=TRUE),c(.025,.975),type=7,names=FALSE))
    ci<-t(apply(bs,1,quantile,probs=c(.025,.975),names=FALSE))
    out[[paste0("direct",n)]]<-stress_row(d,paste0("direct",n),
      quantile(a,c(.025,.975),type=7,names=FALSE),model="Empirical quantiles (type 7)",
      ci=ci,boot_requested=B,boot_successful=B)
  }
  do.call(rbind,out)
}
for(k in if(case==0L)seq_len(nrow(stress_design))else case) {
  if(mode%in%c("repeat","both")) {
    rows<-parallel::mclapply(seq_len(reps),function(r)point_one(k,r),mc.cores=workers,mc.set.seed=FALSE)
    stopifnot(all(vapply(rows,is.data.frame,logical(1))))
    save_case(do.call(rbind,rows),sprintf("assets/stress-grid/cases/repeat-%02d.csv",k))
  }
  if(mode%in%c("bootstrap","both"))save_case(bootstrap_one(k),sprintf("assets/stress-grid/cases/example-%02d.csv",k))
}
writeLines(c(capture.output(sessionInfo()),paste("refineR",if(requireNamespace("refineR",quietly=TRUE))
  as.character(utils::packageVersion("refineR"))else"not used")),sprintf("assets/stress-grid/cases/versions-%02d-%s.txt",case,mode))
