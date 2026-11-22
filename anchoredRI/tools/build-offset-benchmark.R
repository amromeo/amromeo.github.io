# Run from anchoredRI: Rscript tools/build-offset-benchmark.R [workers=6] [bootstrap=200]
# Reuse the exact four published datasets and comparator estimates.
s <- readLines("tools/build-guided-benchmark.R")
eval(parse(text=s[seq_len(which(grepl("^grid <-",s))-1L)]))
B <- if(length(args)>=2L)as.integer(args[2])else 200L
stopifnot(B>=200L)
dir.create("assets/offset",showWarnings=FALSE)
dir.create("tools/offset-cache",showWarnings=FALSE)
source_files <- sort(list.files("R",full.names=TRUE))
hash <- tools::md5sum(source_files)
write.csv(data.frame(file=names(hash),md5=unname(hash)),"assets/offset/source-hashes.csv",row.names=FALSE)
baseline <- read.csv("assets/guided/fits.csv",stringsAsFactors=FALSE)
out <- details <- list()
for(i in seq_len(nrow(profiles))) {
  id <- profiles$id[i]
  cache <- sprintf("tools/offset-cache/%s-B%d.rds",id,B)
  old <- if(file.exists(cache))readRDS(cache)else NULL
  if(!is.null(old) && !identical(old$hash,hash))stop("Cached results use a different estimator; move the old cache first.")
  if(is.null(old)) {
    message("Starting ",id,": ",B," bootstrap fits")
    d <- make_data(id,1L)
    f <- bootstrap_anchored_ri(d$value,d$setting,d$anchor[1:40],"auto",
      n_boot=B,seed=930000L+1000L*i,workers=workers)
    saveRDS(list(hash=hash,result=f),cache)
  } else f <- old$result
  d <- make_data(id,1L)
  row <- baseline[baseline$profile==id & baseline$method=="anchored40",,drop=FALSE]
  stopifnot(nrow(row)==1L)
  row$method <- "anchored_auto40"
  row$model <- paste("Automatic selection:",f$fit$transform)
  row$lower <- f$limits[1]; row$upper <- f$limits[2]
  row$lower_ci_low <- f$ci[1,1]; row$lower_ci_high <- f$ci[1,2]
  row$upper_ci_low <- f$ci[2,1]; row$upper_ci_high <- f$ci[2,2]
  row$lower_error <- f$limits[1]-d$truth[1]; row$upper_error <- f$limits[2]-d$truth[2]
  row$mpe <- 100*mean(abs(f$limits-d$truth)/d$truth)
  row$worst_error <- max(abs(f$limits-d$truth))
  row$reference_coverage <- diff(reference_cdf(f$limits,id))
  row$lower_ci_covers <- f$ci[1,1]<=d$truth[1] && f$ci[1,2]>=d$truth[1]
  row$upper_ci_covers <- f$ci[2,1]<=d$truth[2] && f$ci[2,2]>=d$truth[2]
  row$bootstrap_requested <- B;row$bootstrap_successful <- f$successful;row$bootstrap_failed <- f$failed
  row$at_bound <- f$fit$at_bound;row$status <- "ok";row$error <- ""
  out[[i]] <- row
  details[[i]] <- data.frame(profile=id,selected_transform=f$fit$transform,
    components=f$fit$components,estimated_shift=f$fit$shift,at_bound=f$fit$at_bound,
    shift_at_bound=isTRUE(f$fit$shift_at_bound),
    bootstrap_identity=sum(f$bootstrap$transform=="identity",na.rm=TRUE),
    bootstrap_log=sum(f$bootstrap$transform=="log",na.rm=TRUE),
    bootstrap_shifted_log=sum(f$bootstrap$transform=="shifted_log",na.rm=TRUE),
    bootstrap_boundary=f$boundary_fits,bootstrap_failed=f$failed)
  write.csv(f$bootstrap,paste0("assets/offset/",id,"-bootstrap.csv"),row.names=FALSE,na="")
  write.csv(f$fit$selection,paste0("assets/offset/",id,"-selection.csv"),row.names=FALSE,na="")
  message("Completed ",id,": ",paste(round(f$limits,3),collapse=" to "),"; MPE ",round(row$mpe,3),"%; ",f$successful,"/",B)
}
stopifnot(identical(hash,tools::md5sum(source_files)))
write.csv(rbind(baseline,do.call(rbind,out)),"assets/offset/fits.csv",row.names=FALSE,na="")
write.csv(do.call(rbind,details),"assets/offset/models.csv",row.names=FALSE,na="")
writeLines(c(R.version.string,"anchoredRI 0.2.0",paste("Bootstrap resamples:",B),
  "Original worked datasets and refineR/direct/anchoredRI 0.1 comparator rows preserved.",
  "Automatic family and component-count selection repeated in every bootstrap.",
  "Offset estimated from observations; no true parameters enter the fitter."),"assets/offset/versions.txt")
print(do.call(rbind,details),row.names=FALSE)
