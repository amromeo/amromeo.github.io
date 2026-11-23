# Materialize baseline rows from existing local fits. No refitting or truth scoring.
args <- commandArgs(TRUE)
stopifnot(length(args)==4L)
config <- readRDS(args[1])
arm_b <- args[2]
independent_baseline <- args[3]
out <- args[4]
source("tools/boxcox/reliability-v3/common.R")
source("tools/boxcox/post-round3/comparator.R")
defs <- loadTestsetDefinition()
sources <- r3_sources()
hashes <- do.call(rbind,lapply(c(config$small_archive,config$large_archive),
    function(p) read.csv(file.path(p,"anchor-sample-hashes.csv"))))
saved <- read.csv(file.path(arm_b,"datasets.csv"))
stopifnot(file.exists(file.path(arm_b,"DONE.txt")), !anyDuplicated(saved$Index),nrow(saved)==540L)
rows <- lapply(sort(saved$Index),function(i) {
  fit <- readRDS(file.path(arm_b,sprintf("fit-%04d.rds",i)))
  s <- defs[defs$Index==i,]
  stopifnot(fit$Index==i, identical(fit$sources[names(sources)],sources),nrow(s)==1L)
  directory <- if (s$N<=5000) config$small_data else config$large_data
  d <- bc_ribench_input(s,directory,hashes$sha256[match(i,hashes$Index)])
  routine <- digest::digest(d$value,algo="sha256")
  anchor <- digest::digest(d$anchor,algo="sha256")
  stopifnot(routine==fit$meta$routine_hash,anchor==fit$meta$anchor_hash)
  meta <- data.frame(schema_version=1L,dataset_id=d$key,benchmark="RIbench",
        stage=if (s$N<=5000) "small" else "50000",Index=i,case=NA_integer_,replicate=NA_integer_,
        Analyte=s$Analyte,N=s$N,n_routine=length(d$value),n_anchor=length(d$anchor),
        distribution=s$Distribution,pathological_fraction=s$fractionPathol,scenario=s$Scenario,
        rounding_step=10^(-s$decimals),routine_seed=s$startSeed,anchor_seed=17000000L+i,
        routine_hash=routine,anchor_hash=anchor,ribench_hashcode=s$HashCode,stringsAsFactors=FALSE)
  record <- comparator_row(fit$models,meta,d$value,d$anchor,d$probs)
  if (i %in% seq_len(156L)) {
    original <- readRDS(file.path(independent_baseline,sprintf("fit-%03d.rds",i)))
    original_winner <- original$base_tab[original$base_tab$selected,]
    stopifnot(original$routine==routine,original$anchor==anchor,
              record$transform==original_winner$transform,record$components==original_winner$k,
              record$lower==original_winner$lower,record$upper==original_winner$upper,
              record$loglik==original_winner$loglik,record$bic==original_winner$bic)
  }
  record
})
write_comparator(rows,out)
expected <- do.call(rbind,rows)
check <- read.csv(out,na.strings="",stringsAsFactors=FALSE)
for (key in c("lower","upper","reported_lower","reported_upper","loglik","bic","lrt_p","lrt_statistic"))
  stopifnot(identical(as.double(check[[key]]),as.double(expected[[key]])))
stopifnot(sum(check$stage=="small")==396L,sum(check$stage=="50000")==144L,
          all(is.na(check$reported_upper)==(check$declined|check$technical)))
cat(nrow(check),"unique comparator rows. Exact double-precision CSV round-trip verified.\n")
