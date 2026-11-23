# Load the pinned fitting program through its run_one definition.
program <- parse('tools/boxcox/accuracy/fit.R')
n <- length(program)
stopifnot(identical(program[[n-1]],quote(r3_pool(keys,workers,run_one))))
for (expr in program[seq_len(n-2)]) eval(expr,.GlobalEnv)
ci <- Sys.getenv('CI_FILES')
plan <- read.csv(file.path(ci,'manifest.csv'))
stopifnot(setequal(keys,plan$key),!anyDuplicated(plan$key))
shard <- as.integer(Sys.getenv('BATCH'))
retry_local <- Sys.getenv('PHASE')=='fallback'
use <- if(retry_local) plan$local & !plan$audit else !plan$local | plan$audit
selected <- plan$key[use & plan$shard==shard]
run_measured <- function(key) {
  ok <- run_one(key)
  rss <- grep('^VmHWM:',readLines('/proc/self/status'),value=TRUE)
  writeLines(c(key,rss),file.path(out,paste0(key,'.memory')))
  ok
}
r3_pool(selected,workers,run_measured)
write.csv(data.frame(key=selected,runtime='Actions',run=Sys.getenv('GITHUB_RUN_ID'),attempt=Sys.getenv('GITHUB_RUN_ATTEMPT'),shard),file.path(out,'provenance.csv'),row.names=FALSE)
if(!retry_local) {
 a <- read.csv(file.path(ci,'audit.csv'))
 a <- a[a$key %in% selected,]
 a$pass <- logical(nrow(a))
 for(i in seq_len(nrow(a))){
   x<-a[i,];z<-readRDS(file.path(out,'fits',paste0(x$key,'.rds')))
   stopifnot(z$routine==x$routine,z$anchor==x$anchor)
   tab<-v2_table(Filter(function(m)m$k<=x$cap,z$models),if(x$one_sided).95 else c(.025,.975))
   y<-subset(tab,selected==TRUE)
   span<-if(x$one_sided)x$upper else x$upper-x$lower
   a$pass[i]<-isTRUE(y$transform==x$transform && y$k==x$components && abs(y$upper-x$upper)<=.005*span && (x$one_sided || abs(y$lower-x$lower)<=.005*span))
   a$actions_lower[i]<-y$lower;a$actions_upper[i]<-y$upper;a$actions_transform[i]<-y$transform;a$actions_k[i]<-y$k
 }
 write.csv(a,file.path(out,'audit-result.csv'),row.names=FALSE)
}
