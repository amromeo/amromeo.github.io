library(anchoredRI)
library(digest)
args <- commandArgs(TRUE)
mode <- args[1]
group <- as.integer(args[2])
root <- "source/results/stage3-confirmation-45463873162"
idx <- read.csv(file.path(root,"task-g-index.csv"))
design <- read.csv("source/stage3/design-confirmation.csv")
audit <- read.csv("audit.csv")$j
hash <- function(x) digest(x,algo="sha256",serializeVersion=2L)
archives <- sort(unique(idx$archive))
chosen <- archives[(seq_along(archives)-1L) %% 32L == group]
out <- list()
for (archive in chosen) {
  i <- idx[idx$archive==archive,]
  if(mode=="audit") i <- i[i$j %in% audit,]
  if(!nrow(i)) next
  path <- file.path(root,archive)
  member <- unique(i$member)
  f <- unzip(path,files=member,exdir=tempdir())
  manifest <- jsonlite::fromJSON(file.path(root,"task-g-manifest.json"))
  expected <- manifest$primary_rds
  stopifnot(length(member)==1L,identical(digest(file=f,algo="sha256"),expected$rds_sha256[expected$archive==archive]))
  shard <- readRDS(f)
  for(k in seq_len(nrow(i))) {
    row <- i[k,]
    case <- shard$cases[[row$rds_position]]
    d <- design[match(row$j,design$j),]
    x <- case$routine
    a <- case$panels$guidance
    stopifnot(case$j==row$j,hash(x)==row$routine_hash,hash(a)==row$anchor_hash,length(a)==40L)
    z <- tryCatch(fit_anchored_ri(x,rep("routine",length(x)),a,rounding_unit=d$h,fit_seed=1L),error=identity)
    rec <- data.frame(j=row$j,routine_hash=row$routine_hash,anchor_hash=row$anchor_hash,h=d$h,seed=1L,
       status="fit_failed",reason="",family=NA_character_,components=NA_integer_,optimizer=NA_character_,
       lrt=NA,lower_flag=NA,raw_lower=NA_real_,raw_upper=NA_real_,lower=NA_real_,upper=NA_real_)
    if(inherits(z,"error")) rec$reason <- conditionMessage(z) else {
      rec$status <- z$status
      rec$reason <- paste(z$reason,collapse=" | ")
      rec$family <- z$transform
      rec$components <- z$components
      rec$optimizer <- z$optimizer_status
      rec$lrt <- z$safeguard$lrt_flag
      rec$lower_flag <- z$safeguard$lower_flag
      rec[1,c("raw_lower","raw_upper")] <- as.list(unname(z$diagnostics$anchored$limits))
      rec[1,c("lower","upper")] <- as.list(unname(z$limits))
    }
    out[[length(out)+1L]] <- rec
  }
  if(mode=="audit") {
    chk <- subset(shard$checks,j %in% i$j & arm %in% c("w1","t0.25","baseline") & panel=="guidance" & seed_kind=="primary")
    write.csv(chk,file.path("out",paste0("checks-",basename(archive),".csv")),row.names=FALSE)
  }
}
write.csv(do.call(rbind,out),sprintf("out/%s-%02d.csv",mode,group),row.names=FALSE)
capture.output(sessionInfo(),file="out/environment.txt")
