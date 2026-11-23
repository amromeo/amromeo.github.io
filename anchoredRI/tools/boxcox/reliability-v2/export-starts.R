# Per-start convergence records. The fit RDS files also retain parameter vectors and boxes.
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
files <- sort(list.files(args[1],'^fit-[0-9]+[.]rds$',full.names=TRUE))
stopifnot(length(files)==156)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
for(block in 1:26) {
 rows <- list()
 for(file in files[(block-1)*6+1:6]) {
  r <- readRDS(file)
  for(m in r$models) {
   rows[[length(rows)+1L]] <- data.frame(Index=r$Index,runtime=r$runtime,transform=m$transform,k=m$k,
    seed=m$seed,start=seq_along(m$starts),chosen=seq_along(m$starts)==m$origin,
    source=vapply(m$starts,`[[`,'','label'),initial_loglik=vapply(m$starts,`[[`,0.0,'initial'),
    optimizer_loglik=vapply(m$starts,`[[`,0.0,'final'),retained_loglik=vapply(m$starts,`[[`,0.0,'retained'),
    convergence=vapply(m$starts,`[[`,0L,'code'),message=vapply(m$starts,`[[`,'','message'),
    projected=vapply(m$starts,`[[`,FALSE,'projected'),initial_retained=vapply(m$starts,`[[`,FALSE,'initial_retained'),
    at_bound=vapply(m$starts,`[[`,FALSE,'at_bound'))
  }
 }
 name <- file.path(args[2],sprintf('part-%02d.csv.gz',block)); handle <- gzfile(name,'wt')
 write.csv(do.call(rbind,rows),handle,row.names=FALSE);close(handle)
}
