# Fit Arm B from scratch and cross-start the two original Hb 142 baseline winners.
args <- commandArgs(TRUE)
stopifnot(length(args)==7L, args[1] %in% c("local","Actions"))
runtime <- args[1]
data_dir <- args[2]
archive <- args[3]
old_dir <- args[4]
out <- args[5]
workers <- as.integer(args[6])
so <- args[7]
study <- "tools/boxcox/post-round3"
source("tools/boxcox/reliability-v3/common.R")
source(file.path(study,"hb142.R"))
dyn.load(so)
sources <- r3_pedigree(runtime,archive,out)
hash_files <- function(x) vapply(x,function(p) digest::digest(file=p,algo="sha256"),character(1))
frozen <- read.csv("tools/boxcox/reliability-v3/results/source-hashes.csv")
stopifnot(all(hash_files(frozen$file)==frozen$sha256))
native <- read.csv(paste0(so,".sources.csv"))
stopifnot(all(native$sha256==frozen$sha256[match(native$file,frozen$file)]))
extra <- file.path(study,c("hb142.R","hb142-starts.csv","replay.R"))
manifest <- c(setNames(frozen$sha256,frozen$file),hash_files(extra))
write.csv(data.frame(file=names(manifest),sha256=unname(manifest)),file.path(out,"sources.csv"),row.names=FALSE)
defs <- loadTestsetDefinition()
hashes <- read.csv(file.path(archive,"anchor-sample-hashes.csv"))
donors <- read.csv(file.path(study,"hb142-starts.csv"),colClasses="character")
donors$k <- as.integer(donors$k)
run_one <- function(i) {
  path <- file.path(out,sprintf("fit-%03d.rds",i))
  s <- defs[defs$Index==i,]
  stopifnot(nrow(s)==1L, s$Analyte=="Hb", s$N<=5000)
  input <- r3_inputs(s,data_dir,hashes)
  prior <- readRDS(file.path(old_dir,sprintf("fit-%03d.rds",i)))
  stopifnot(prior$runtime==runtime, prior$Index==i, identical(prior$sources,sources),
            prior$routine==input$routine, prior$anchor==input$anchor)
  if (file.exists(path)) {
    v <- readRDS(path)
    stopifnot(identical(v$manifest,manifest), v$runtime==runtime,
              v$routine==input$routine, v$anchor==input$anchor, v$Index==i)
    return(TRUE)
  }
  d <- input$d
  models <- r3_baseline(d,i,s$decimals,kmax=5L)
  table <- r3_bt(models,d,r3_rule(d,s$decimals),input$probs)
  table$method <- "armB"
  stopifnot(nrow(table)==15L, sum(table$selected)==1L)
  subset <- table[table$k<=3L,]
  compare <- merge(prior$base_tab,subset,by=c("transform","k"),suffixes=c("_round3","_replay"))
  stopifnot(nrow(compare)==9L)
  repeated <- data.frame(Index=i,runtime=runtime,
    same_baseline_winner=paste(subset$transform[which.min(subset$bic)],subset$k[which.min(subset$bic)])==
                         paste(prior$base_tab$transform[prior$base_tab$selected],prior$base_tab$k[prior$base_tab$selected]),
    max_change_loglik=max(abs(compare$loglik_round3-compare$loglik_replay)),
    max_change_endpoint=max(abs(compare$lower_round3-compare$lower_replay),abs(compare$upper_round3-compare$upper_replay)))
  cross <- if (i==142L) hb142_cross(d,input$routine,input$anchor,s$decimals,runtime,donors)
  z <- list(Index=i,runtime=runtime,manifest=manifest,routine=input$routine,anchor=input$anchor,
            models=models,table=table,nesting=r3_nesting(models,d),zero=r3_zero(d,models,r3_rule(d,s$decimals)),
            repeatability=repeated,cross=cross)
  r3_save(z,path)
  cat("saved",runtime,i,"\n"); flush.console()
  TRUE
}
r3_pool(seq_len(156L),workers,run_one)
fits <- lapply(seq_len(156L),function(i) readRDS(file.path(out,sprintf("fit-%03d.rds",i))))
for (name in c("table","nesting","zero")) {
  rows <- do.call(rbind,lapply(fits,function(x) cbind(Index=x$Index,runtime=x$runtime,x[[name]])))
  write.csv(rows,file.path(out,paste0(name,".csv")),row.names=FALSE)
}
write.csv(do.call(rbind,lapply(fits,`[[`,"repeatability")),file.path(out,"baseline-repeatability.csv"),row.names=FALSE)
meta <- do.call(rbind,lapply(fits,function(x) data.frame(Index=x$Index,runtime=x$runtime,
                                               routine_hash=x$routine,anchor_hash=x$anchor)))
write.csv(meta,file.path(out,"input-hashes.csv"),row.names=FALSE)
for (field in c("fits","blocks"))
  write.csv(fits[[142]]$cross[[field]],file.path(out,paste0("hb142-cross-",field,".csv")),row.names=FALSE)
writeLines("156 Arm B fits and Hb 142 reciprocal starts complete",file.path(out,"DONE.txt"))
