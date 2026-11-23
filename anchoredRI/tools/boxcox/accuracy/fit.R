args <- commandArgs(TRUE)
stopifnot(length(args)==5L)
cfg <- readRDS(args[1])
out <- args[2]
cache <- args[3]
workers <- as.integer(args[4])
source("tools/boxcox/reliability-v3/common.R")
source("tools/stress-design.R")
dyn.load(args[5])
dir.create(file.path(out,"fits"),recursive=TRUE,showWarnings=FALSE)
defs <- loadTestsetDefinition()
defs <- defs[defs$N<=5000,]
hashes <- read.csv(file.path(cfg$small_archive,"anchor-sample-hashes.csv"))
core <- r3_sources()
old <- read.csv("tools/boxcox/reliability-v3/results/source-hashes.csv")
stopifnot(all(vapply(old$file,function(f) digest::digest(file=f,algo="sha256"),"")==old$sha256))
manifest <- c(core,design=digest::digest(file="tools/stress-design.R",algo="sha256"),
              run=digest::digest(file="tools/boxcox/accuracy/fit.R",algo="sha256"))
grid <- expand.grid(replicate=1:20,case=1:45)
keys <- c(paste0("r-",defs$Index),sprintf("s-%02d-%02d",grid$case,grid$replicate))
stopifnot(length(unique(keys))==3780L)
run_one <- function(key) {
  f <- file.path(out,"fits",paste0(key,".rds"))
  if(file.exists(f)) {
    z <- readRDS(f)
    stopifnot(z$key==key,identical(z$manifest,manifest))
    return(TRUE)
  }
  begin <- as.numeric(Sys.time())
  fields <- strsplit(key,"-",fixed=TRUE)[[1]]
  cached <- FALSE
  if(fields[1]=="r") {
    i <- as.integer(fields[2]); s <- defs[defs$Index==i,]
    stopifnot(nrow(s)==1L,s$N<=5000)
    d <- bc_ribench_input(s,cfg$small_data,hashes$sha256[match(i,hashes$Index)])
    decimals <- s$decimals; index <- i
    previous <- file.path(cache,sprintf("fit-%04d.rds",i))
    cached <- file.exists(previous)
  } else {
    d <- stress_data(as.integer(fields[2]),as.integer(fields[3]))
    decimals <- 2L; index <- d$seed
  }
  routine <- digest::digest(d$value,algo="sha256"); anchor <- digest::digest(d$anchor,algo="sha256")
  t0 <- proc.time()
  if(cached) {
    r <- readRDS(previous)
    stopifnot(r$Index==i,identical(r$sources[names(core)],core),
              r$meta$routine_hash==routine,r$meta$anchor_hash==anchor)
    models <- r$models
  } else {
    model_data <- v2_data(d$value,d$setting,d$anchor)
    models <- r3_baseline(model_data,index,decimals,5L)
  }
  dt <- proc.time()-t0
  stopifnot(max(vapply(models,function(m) as.integer(m$k),1L))==5L)
  z <- list(key=key,manifest=manifest,routine=routine,anchor=anchor,models=models,
            cached=cached,began=begin,ended=as.numeric(Sys.time()),elapsed=unname(dt["elapsed"]),
            cpu=unname(dt["user.self"]+dt["sys.self"]))
  r3_save(z,f)
  TRUE
}
set.seed(3245);keys <- sample(keys)
r3_pool(keys,workers,run_one)
writeLines("3780 checkpoints verified",file.path(out,"FITTED.txt"))
