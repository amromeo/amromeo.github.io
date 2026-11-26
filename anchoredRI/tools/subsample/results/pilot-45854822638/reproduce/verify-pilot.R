# Independent readback, R 4.6.x. Uses no fitting/sampling code or external packages.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
artifacts <- normalizePath(args[1])
input <- normalizePath(args[2])
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
repo <- normalizePath(".")
study <- file.path(repo,"anchoredRI/tools/subsample")
roster <- read.csv(file.path(study,"pilot-roster.csv"),stringsAsFactors=FALSE)
index <- read.csv(file.path(input,"task-g-index.csv"),stringsAsFactors=FALSE)
design <- read.csv(file.path(input,"design-large.csv"),stringsAsFactors=FALSE)
readtab <- function(d, name) read.csv(file.path(d,paste0(name,".csv")),stringsAsFactors=FALSE,
                                      na.strings=c("","NA"),check.names=FALSE)
digest_v2 <- function(x) unname(tools::sha256sum(bytes=serialize(x,NULL,version=2)[-(1:14)]))
# CSVs were exported with %.17g. Numeric readback across R versions may differ by an ulp.
same <- function(x,y,eps=5e-13) {
  if (!length(x)) x <- NA
  if (!length(y)) y <- NA
  if (length(x)!=length(y) || !identical(as.vector(is.na(x)),as.vector(is.na(y)))) return(FALSE)
  i <- !is.na(x)
  if (!any(i)) return(TRUE)
  if (is.numeric(x) || is.numeric(y))
    return(all(abs(as.numeric(x[i])-as.numeric(y[i])) <= eps * pmax(1,abs(as.numeric(y[i])))))
  all(as.character(x[i]) == as.character(y[i]))
}
fields <- function(a,b,keys) all(vapply(keys,function(k) same(a[[k]],b[[k]]),logical(1)))
valid <- function(l,u,crp) is.finite(u) & (crp | (is.finite(l) & l < u))
writetab <- function(rows,path) write.csv(rows,file.path(out,path),row.names=FALSE,na="")
markers <- list.files(artifacts,pattern="^COMPLETE$",recursive=TRUE,full.names=TRUE)
dirs <- dirname(markers)
ids <- suppressWarnings(as.integer(basename(dirs)))
stopifnot(length(ids)==20, !anyDuplicated(ids), setequal(ids,roster$Index),
          !anyDuplicated(roster[c("Analyte","pathology_stratum")]),all(roster$eligible_datasets == 72))
names(dirs) <- ids
freeze <- lapply(c("43c8834","f1f230a"),function(rev) {
  lines <- system2("git", c("-C",repo,"show",paste0(rev,":anchoredRI/tools/subsample/SOURCE.sha256")),stdout=TRUE)
  parts <- strsplit(lines,"  ",fixed=TRUE)
  setNames(vapply(parts,`[`,character(1),1L),vapply(parts,`[`,character(1),2L))
})
records <- anomalies <- failures <- list()
context <- ""
check <- function(pass,what) {
  if (!isTRUE(pass)) stop(context," | ",what,call.=FALSE)
}
odd <- function(id,kind,count,detail="") {
  if (count) anomalies[[length(anomalies)+1L]] <<- data.frame(Index=id,kind=kind,count=count,detail=detail)
}
for (id in roster$Index) {
  context <- paste("Index",id)
  result <- tryCatch({
    d <- dirs[as.character(id)]
    meta <- design[match(id,design$Index),]
    i <- index[match(id,index$Index),]
    r <- roster[match(id,roster$Index),]
    src <- readRDS(file.path(input,sprintf("taskg-shard-%d.rds",i$shard)))$cases[[i$rds_position]]
    routine <- src$routine
    anchor <- src$panels$guidance
    n <- length(routine)
    rh <- digest_v2(routine)
    ah <- digest_v2(anchor)
    crp <- meta$Analyte == "CRP"
    check(src$Index == id && src$j == id && n == i$routine_n && length(anchor)==40L &&
          rh == i$routine_hash && rh == meta$routine_hash && rh == r$routine_hash &&
          ah == i$anchor_hash && ah == meta$anchor_hash && ah == r$anchor_hash,"archive identities")
    version <- readRDS(file.path(d,"input-id.rds"))
    check(version$Index == id && version$routine_hash == rh && version$anchor_hash == ah &&
          startsWith(version$runtime,"R version 4.5.0") && version$anchoredRI=="0.4.0","runtime receipt")
    source_files <- file.path("anchoredRI/tools/subsample",c("pilot.R","aggregate.R","inputs.R","mixture.R"))
    check(identical(names(version$code),source_files) &&
          any(vapply(freeze,function(manifest) all(manifest[source_files]==version$code),logical(1))),
          "four pilot sources jointly match a frozen implementation")
    f <- readtab(d,"fits"); dr <- readtab(d,"draws")
    m <- readtab(d,"medians"); v <- readtab(d,"decisions"); tm <- readtab(d,"timing")
    models <- readRDS(file.path(d,"models.rds"))
    pos <- readRDS(file.path(d,"positions.rds"))
    check(all(vapply(list(f,dr,m,v,tm),function(x) all(x$Index==id),logical(1))),"dataset IDs across CSVs")
    check(length(models)==204L && nrow(f)==204L && length(pos)==100L && nrow(dr)==100L &&
            nrow(m)==8L && nrow(v)==48L,"export sizes")
    expected_draws <- as.vector(outer(c("A","B"),sprintf("%02d",1:50),paste0))
    labels <- function(x) sprintf("%s%02d",x$seed_set,x$draw)
    check(!anyDuplicated(names(pos)) && !anyDuplicated(labels(dr)) && !anyDuplicated(dr$position_hash) &&
          setequal(names(pos),expected_draws) && setequal(labels(dr),expected_draws),"100 distinct draws")
    dr <- dr[match(names(pos),labels(dr)),]
    for (j in seq_along(pos)) {
      p <- pos[[j]]
      expected_seed <- (if(dr$seed_set[j]=="A") 100000000L else 200000000L) + 100L*id + dr$draw[j]
      check(is.integer(p) && length(p)==5000L && !is.unsorted(p, strictly=TRUE) &&
            all(p>=1L & p<=n) && digest_v2(p)==dr$position_hash[j] &&
            dr$sampling_seed[j]==expected_seed && dr$positions_n[j]==5000,"position order/bounds/hash/draw seed")
    }
    fit_key <- function(t) with(t,paste(cap,fit_type,seed_set,draw,sep="/"))
    k <- fit_key(f)
    check(!anyDuplicated(k) && setequal(f$cap,c(3L,5L)) && all(f$routine_hash==rh & f$anchor_hash==ah) &&
            all(f$fit_seed == 1L) && all(f$rounding_unit==meta$h),"fit keys/identity")
    used <- character()
    p_error <- sum_error <- 0
    for(j in seq_along(models)) {
      entry <- models[[j]]
      key <- fit_key(entry$key)
      row <- f[match(key,k),]
      fit <- entry$fit
      context <- paste("Index",id,key)
      check(!is.na(row$Index) && fields(entry$key,row,names(entry$key)),"model key and CSV")
      used <- c(used,key)
      sub <- row$fit_type=="subsample"
      obs <- if(sub) 5000L else n
      q <- if(row$fit_type=="W") 5000/n else 1
      est_n <- if(row$fit_type=="full") n+40L else 5040L
      x <- routine
      if(sub) {
        draw_key <- labels(row)
        matching <- dr[match(draw_key,names(pos)),]
        check(same(matching$sampling_seed,row$sampling_seed) &&
              same(matching$position_hash,row$position_hash),"subsample identity against original positions")
        x <- x[pos[[draw_key]]]
      } else {
        check(is.na(row$sampling_seed) && is.na(row$position_hash) && is.na(row$seed_set) &&
              is.na(row$draw) && row$fit_type %in% c("W","full"),"full/weighted row identity")
      }
      check(same(obs,row$routine_n) && same(q,row$routine_weight) && same(est_n,row$bic_n) &&
            row$cap %in% c(3L,5L),"sample mass, cap and BIC sample size")
      ft <- fit$diagnostics$anchored$limits
      expected <- list(first_status=fit$anchoredRI_status,status=fit$status,method=fit$method,
                       components=fit$components,transform=fit$transform,shift=fit$shift,
                       mu=fit$reference$mu,sigma=fit$reference$sigma,
                       p=fit$safeguard$lrt_p,lrt_flag=fit$safeguard$lrt_flag,lower_flag=fit$safeguard$lower_flag,
                       lower_check_skipped=fit$safeguard$lower_check_skipped,routine_min=fit$safeguard$routine_min,
                       raw_lower=unname(ft["lower"]),raw_upper=unname(ft["upper"]),
                       lower=unname(fit$limits["lower"]),upper=unname(fit$limits["upper"]))
      check(fields(expected,row,names(expected)),"saved model vs published status, checks, parameters and endpoints")
      if (length(fit$weights)) {
        w <- as.vector(fit$weights)
        check(length(w)==row$components && row$components>=1 && row$components<=row$cap &&
              all(is.finite(w) & w>=0) && abs(sum(w)-1)<=1e-10 &&
              same(w[1],row$share) && same(sum(w),row$weight_sum),"component shares")
        sum_error <- max(sum_error,abs(sum(w)-1))
      } else check(is.na(row$share) && !is.na(row$share_reason),"missing share reason")
      if(row$first_status!="withheld") {
        y <- switch(row$transform,identity=anchor,log=log(anchor),shifted_log=log(anchor-row$shift))
        sig <- sqrt(mean((y-mean(y))^2))
        st <- max(0,2*sum(dnorm(y,mean(y),sig,log=TRUE)-dnorm(y,row$mu,row$sigma,log=TRUE)))
        p <- pchisq(st,df=2,lower.tail=FALSE)
        p_error <- max(p_error,abs(p-row$p))
        lower_flag <- !crp && row$raw_lower < min(x)
        flag <- p < .01
        check(abs(p-row$p)<1e-9 && same(row$lrt_statistic,st,1e-10) && flag==row$lrt_flag &&
              lower_flag==row$lower_flag && row$lower_check_skipped==crp && same(min(x),row$routine_min) &&
              row$first_status==if(flag||lower_flag) "declined" else "reported","recomputed anchor LRT, lower safeguard or status")
      }
      if(row$status=="reported") check(valid(row$lower,row$upper,crp),"released limits")
      else check(is.na(row$lower) && is.na(row$upper),"withheld limits")
    }
    check(!anyDuplicated(used) && setequal(used,k),"model-to-row bijection")
    standard <- f[f$cap==5 & f$fit_type=="full",]
    check(nrow(standard)==1L && nrow(m[!duplicated(m[c("cap","seed_set","definition")]),])==8L &&
          !anyDuplicated(v[c("cap","seed_set","definition","estimator","rule")]),"decision keys")
    for (cap in c(3L,5L)) for (bank in c("A","B")) for (definition in c("all","passing")) {
      context <- paste("Index",id,cap,bank,definition)
      subset_fit <- f[f$cap==cap,]
      zz <- subset_fit[subset_fit$fit_type=="subsample" & subset_fit$seed_set %in% bank,]
      full <- subset_fit[subset_fit$fit_type=="full",]; w <- subset_fit[subset_fit$fit_type=="W",]
      check(nrow(full)==1L && nrow(w)==1L && identical(sort(zz$draw),1:50),"fit and decision coverage")
      med <- m[m$cap==cap & m$seed_set==bank & m$definition==definition,]
      decisions <- v[v$cap==cap & v$seed_set==bank & v$definition==definition,]
      ndecl <- sum(zz$first_status=="declined"); ntech <- sum(zz$first_status=="withheld")
      gate <- if (ndecl>=26) "declined" else if (ndecl+ntech<=25) "reported" else "withheld"
      check(nrow(med)==1 && nrow(decisions)==6 && med$votes_decline==ndecl && med$votes_technical==ntech &&
            med$votes_pass==50-ndecl-ntech && all(unlist(med[paste0("k",1:5)])==tabulate(zz$components,5)),"votes and component histogram")
      keep <- if(definition=="all") rep(TRUE,50L) else zz$first_status=="reported"
      use <- keep & valid(zz$raw_lower,zz$raw_upper,crp)
      shares <- keep & is.finite(zz$share) & zz$share>=0 & zz$share<=1
      matched <- use & shares
      mid <- function(x) if(length(x)) median(x) else NA_real_
      vals <- list(upper=mid(zz$raw_upper[use]),lower=if(crp) NA_real_ else mid(zz$raw_lower[use]),
                   contributors=sum(use),share=mid(zz$share[shares]),share_n=sum(shares),
                   matched_share=mid(zz$share[matched]),matched_share_n=sum(matched),
                   matched_lower=if(crp) NA_real_ else mid(zz$raw_lower[matched]),matched_upper=mid(zz$raw_upper[matched]))
      check(fields(vals,med,names(vals)),"median share, limits, usable/matched/pass contributions")
      for (j in seq_len(nrow(decisions))) {
        drow <- decisions[j,]
        l <- if(drow$estimator=="median") vals$lower else w$raw_lower
        u <- if(drow$estimator=="median") vals$upper else w$raw_upper
        status <- switch(drow$rule,own=w$first_status,full=full$first_status,
                         majority=gate,majority_or_package=gate,stop("unrecognized rule"))
        if(!valid(l,u,crp)) status <- "withheld"
        route <- drow$estimator
        if(drow$rule=="majority_or_package" && status!="reported") {
          status <- standard$status; l<-standard$lower; u<-standard$upper
          route <- if(status!="reported") "full_nonreturn" else if(standard$method=="anchoredRI") "full_first_stage" else "full_fallback"
        }
        check(same(drow$raw_lower,l) && same(drow$raw_upper,u) && drow$route==route &&
              drow$status==status && same(drow$lower,if(status=="reported") l else NA) &&
              same(drow$upper,if(status=="reported") u else NA),"release rule, endpoints and routing")
      }
    }
    context <- paste("Index",id)
    time <- tm[tm$phase=="dataset",]
    check(nrow(time)==1L && time$wall_seconds>0 && time$cpu_seconds>0,"dataset timing")
    unusual_optimizer <- !is.na(f$optimizer_code) & f$optimizer_code!=0
    odd(id,"technical first-stage fit",sum(f$first_status=="withheld"))
    odd(id,"fit with missing share",sum(!is.finite(f$share)))
    odd(id,"optimizer nonzero exit",sum(unusual_optimizer))
    odd(id,"zero passing fits",sum(m$definition=="passing" & m$contributors==0))
    odd(id,"incomplete all-fit median",sum(m$definition=="all" & m$contributors<50))
    odd(id,"public fallback withheld",sum(!is.na(standard$fallback_status) & standard$status!="reported"))
    odd(id,"majority-changing bank",sum(with(m[m$definition=="all",],tapply(votes_decline>=26,cap,function(x)x[1]!=x[2]))))
    bank_meet <- vapply(1:50,function(b) length(intersect(pos[[sprintf("A%02d",b)]],pos[[sprintf("B%02d",b)]])),integer(1))
    rec <- data.frame(Index=id,Analyte=r$Analyte,stratum=r$pathology_stratum,routine_n=n,share_rows=sum(is.finite(f$share)),
           models=nrow(f),draws=nrow(dr),position_indices=sum(lengths(pos)),sample_seed_checks=nrow(dr),medians=nrow(m),
           decisions=nrow(v),max_recomputed_p_diff=p_error,max_abs_component_sum_error=sum_error,
           max_set_overlap=max(bank_meet),min_set_overlap=min(bank_meet),iteration_exits=sum(unusual_optimizer),
           full_first=standard$first_status,full_final=standard$status,full_method=standard$method,
           fallback_status=standard$fallback_status,technical=sum(f$first_status=="withheld"),total_cpu_s=time$cpu_seconds,
           total_wall_s=time$wall_seconds)
    records[[length(records)+1L]] <- rec
    message("Verified ",id)
    TRUE
  },error=function(e){
    message(conditionMessage(e))
    failures[[length(failures)+1L]] <<- data.frame(Index=id,detail=conditionMessage(e))
    FALSE
  })
}
empty <- data.frame(Index=integer(),kind=character(),count=integer(),detail=character())
writetab(if(length(anomalies)) do.call(rbind,anomalies) else empty, "observations.csv")
writetab(if(length(failures)) do.call(rbind,failures) else empty[c("Index","detail")],"failures.csv")
if(length(records)) writetab(do.call(rbind,records),"verified-datasets.csv")
if(length(failures)) stop(length(failures)," pilot datasets failed independent readback")
cat("Verified all",length(records),"datasets\n")
