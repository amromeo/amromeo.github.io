# Audit both binary replay artifacts and export comparable candidate, limit and likelihood tables.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
work <- args[1]
round3 <- args[2]
out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
hosts <- c("local","Actions")
dirs <- setNames(file.path(work,c("local","actions")),hosts)
old_dirs <- setNames(file.path(round3,c("local","actions/phase1")),hosts)
stopifnot(all(file.exists(file.path(dirs,"DONE.txt"))),
          grepl("C2 image match TRUE",readLines(file.path(dirs["Actions"],"runtime.txt"),n=1L),fixed=TRUE))
sources <- lapply(dirs,function(p) read.csv(file.path(p,"sources.csv")))
stopifnot(identical(sources$local,sources$Actions))
write.csv(sources$local,file.path(out,"source-hashes.csv"),row.names=FALSE)
for (host in hosts) {
  writeLines(sub("[[:space:]]+$","",readLines(file.path(dirs[[host]],"runtime.txt"))),
             file.path(out,paste0(host,"-runtime.txt")))
}
save_table <- function(tab,name) {
  for (key in names(tab)[vapply(tab,is.double,logical(1))]) {
    values <- tab[[key]]
    tab[[key]] <- ifelse(is.na(values), NA_character_, sprintf("%.17g",values))
  }
  write.csv(tab,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
}
fit_path <- function(p,i) file.path(p,sprintf("fit-%03d.rds",i))
collect <- function(x) do.call(rbind,x)
pairs <- nesting <- checked <- candidates <- zero <- again <- list()
for (i in seq_len(156L)) {
  new <- lapply(dirs,function(path) readRDS(fit_path(path,i)))
  prior <- lapply(old_dirs,function(path) readRDS(fit_path(path,i)))
  sha <- function(record) c(record$routine,record$anchor)
  stopifnot(all(vapply(c(new,prior),function(f) identical(sha(f),sha(prior$local)),logical(1))))
  model <- list()
  for (h in hosts) {
    fit <- new[[h]]
    stopifnot(fit$runtime==h,fit$Index==i,identical(unname(fit$manifest),sources[[h]]$sha256))
    t <- fit$table
    key <- paste(i,h)
    eligible <- t$eligible & is.finite(t$bic)
    pick <- t[t$selected,]
    stopifnot(nrow(t)==15L, nrow(pick)==1L, pick$bic==min(t$bic[eligible]), pick$k<=5L,
              all(t$minimum_sd_to_floor>=1-1e-11))
    model[[h]] <- pick
    candidates[[key]] <- cbind(Index=i,runtime=h,t)
    ch <- fit$nesting
    expect <- with(ch,status=="checked" & nested_loglik-loglik>roundoff)
    stopifnot(!anyNA(expect),identical(expect,ch$violation))
    nesting[[key]] <- cbind(Index=i,runtime=h,ch)
    checked[[key]] <- data.frame(Index=i,runtime=h,checked=sum(ch$status=="checked"),
            missing=sum(ch$status=="unavailable"),excluded=sum(ch$status=="structurally excluded"),
            violations=sum(ch$violation),min_sd_to_floor=min(t$minimum_sd_to_floor))
    zero[[key]] <- cbind(Index=i,runtime=h,fit$zero)
    again[[key]] <- fit$repeatability
  }
  l <- model$local
  a <- model$Actions
  width <- a$upper-a$lower
  finite <- all(is.finite(c(l$lower,l$upper,a$lower,a$upper))) && width>0
  largest <- if(finite) max(abs(l$lower-a$lower),abs(l$upper-a$upper)) else NA_real_
  pairs[[i]] <- data.frame(Index=i,routine_hash=new$local$routine,anchor_hash=new$local$anchor,
        local_transform=l$transform,local_components=l$k,local_loglik=l$loglik,local_bic=l$bic,
        local_lower=l$lower,local_upper=l$upper,actions_transform=a$transform,actions_components=a$k,
        actions_loglik=a$loglik,actions_bic=a$bic,actions_lower=a$lower,actions_upper=a$upper,
        width_actions=width,absolute_difference_lower=abs(l$lower-a$lower),absolute_difference_upper=abs(l$upper-a$upper),
        max_absolute_difference=largest,fraction_of_width=largest/width,finite_pair=finite,
        model_match=l$transform==a$transform && l$k==a$k,
        limits_within_0_5pct=finite && largest<=.005*width)
}
by_dataset <- collect(pairs)
stopifnot(identical(by_dataset$Index,seq_len(156L)))
save_table(by_dataset,"arm-b-pairs")
exceptions <- subset(by_dataset,!model_match | !limits_within_0_5pct)
save_table(exceptions,"arm-b-exceptions")
audit <- collect(checked)
save_table(audit,"arm-b-checks-by-dataset")
tab <- collect(candidates)
nc <- collect(nesting)
ze <- collect(zero)
for (h in hosts) {
  save_table(tab[tab$runtime==h,],paste0("arm-b-candidates-",h))
  save_table(nc[nc$runtime==h,],paste0("arm-b-nesting-",h))
}
save_table(ze,"zero-shift-vectors")
save_table(collect(again),"baseline-repeatability")
summary <- collect(lapply(hosts,function(h) {
  k <- audit[audit$runtime==h,]
  z <- ze[ze$runtime==h,]
  stopifnot(nrow(k)==156L, all(z$absolute_difference<=z$roundoff))
  data.frame(runtime=h,paired_datasets=nrow(by_dataset),model_matches=sum(by_dataset$model_match),
       limit_matches=sum(by_dataset$limits_within_0_5pct),missing_limit_pairs=sum(!by_dataset$finite_pair),
       eligible_nesting_checks=sum(k$checked),excluded=sum(k$excluded),unavailable=sum(k$missing),
       violations=sum(k$violations),flagged_datasets=sum(k$violations>0L),
       zero_vectors=nrow(z),max_zero_loglik_difference=max(z$absolute_difference),
       max_zero_floor_difference=max(z$physical_floor_difference))
}))
passed <- all(summary$model_matches>=152L & summary$limit_matches>=149L &
              summary$violations==0L & summary$unavailable==0L)
summary$decision <- if(passed) "passed" else "stopped"
save_table(summary,"arm-b-gate")
# Hb 142 baseline analysis. Use the saved round 3 baselines to preserve the original selections.
old142 <- lapply(old_dirs,function(path) readRDS(fit_path(path,142L))$base_tab)
start_tab <- collect(lapply(hosts,function(h) cbind(Index=142L,runtime=h,old142[[h]])))
save_table(start_tab,"hb142-original-candidates")
raw_cross <- lapply(dirs,function(path) readRDS(fit_path(path,142L))$cross)
solutions <- collect(lapply(raw_cross,`[[`,"fits"))
blocks <- collect(lapply(raw_cross,`[[`,"blocks"))
keys <- with(solutions,paste(receiver,donor,transform))
stopifnot(nrow(solutions)==8L,!anyDuplicated(keys),all(solutions$k==3L),
          all(c("identity","log") %in% solutions$transform))
for (idx in seq_len(nrow(solutions))) {
  s <- solutions[idx,]
  original <- old142[[s$donor]]
  start <- original[original$transform==s$transform & original$k==s$k,]
  stopifnot(nrow(start)==1,abs(start$loglik-s$initial_loglik)<1e-10,
            abs(start$lower-s$initial_lower)<1e-10,abs(start$upper-s$initial_upper)<1e-10)
}
save_table(solutions,"hb142-cross-starts")
save_table(blocks,"hb142-optimizer-blocks")
final <- collect(lapply(hosts,function(h) {
  base <- old142[[h]][,c("transform","k","loglik","bic","lower","upper")]
  base$source <- "original round 3"
  x <- subset(solutions,receiver==h & converged)
  crossed <- data.frame(transform=x$transform,k=x$k,loglik=x$best_loglik,bic=x$bic,
                        lower=x$lower,upper=x$upper,source=paste(x$donor,"start"))
  all <- rbind(base,crossed)
  pick <- all[which.min(all$bic),]
  cbind(Index=142L,runtime=h,pick)
}))
save_table(final,"hb142-best-witnessed")
print(summary,row.names=FALSE)
print(final,row.names=FALSE,digits=10)
