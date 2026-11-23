# Inspect stored starts and boxes for observed inversions. No fit is repeated.
args <- commandArgs(TRUE);stopifnot(length(args)==4L)
pairs <- read.csv(gzfile(args[3])); pairs <- pairs[pairs$violation==1L,]
rows <- lapply(seq_len(nrow(pairs)),function(i) {
 v <- pairs[i,];dir <- if(v$runtime=='local')args[1] else args[2]
 fit <- readRDS(file.path(dir,sprintf('fit-%03d.rds',v$Index)))
 small <- fit$models[[paste(v$nested,v$nested_k)]]
 large <- fit$models[[paste(v$transform,v$k)]]
 nested_label <- if(v$comparison=='adjacent K')'exact widest duplicate' else if(v$transform=='shifted_log')'unshifted log' else
   switch(as.character(v$nested),log='fixed log',identity='fixed identity',shifted_log='fixed shifted log',boxcox_unshifted='unshifted Box-Cox')
 starts <- Filter(function(s)s$label==nested_label,large$starts)
 ans <- data.frame(Index=v$Index,runtime=v$runtime,method=v$method,transform=v$transform,k=v$k,
   nested_transform=v$nested,nested_k=v$nested_k,loglik_loss=small$ll-large$ll,
   nested_start=nested_label,start_saved=length(starts)>0,min_nested_sd=min(small$p$sigma),
   minimum_allowed_sd=NA_real_,max_nested_sd=max(small$p$sigma),maximum_allowed_sd=NA_real_,
   projected=NA,loss_from_projection=NA_real_,sd_outside_box=NA)
 if(length(starts)) {
  s <- starts[[1]]
  ans$minimum_allowed_sd <- exp(s$lo[v$k+1])*s$scale
  ans$maximum_allowed_sd <- exp(s$hi[v$k+1])*s$scale
  ans$projected <- s$projected
  ans$loss_from_projection <- small$ll-s$initial
  ans$sd_outside_box <- ans$min_nested_sd<ans$minimum_allowed_sd||ans$max_nested_sd>ans$maximum_allowed_sd
 }
 ans
})
if(length(rows))write.csv(do.call(rbind,rows),args[4],row.names=FALSE)
