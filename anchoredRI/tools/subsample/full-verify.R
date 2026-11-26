# SPDX-License-Identifier: GPL-3.0-only
# Independent readback of saved Task H records, adapted from the frozen pilot
# verifier. The script uses no fit code. Replaying the sampler checks the
# actual selected positions, including their order and ties.
args <- commandArgs(TRUE)
if(length(args)!=5L) stop("Usage: full-verify.R RECORDS_ROOT INPUTS_DIR OUT_DIR pilot|full ID[,ID...]")
records_root <- normalizePath(args[[1L]],mustWork=TRUE)
input <- normalizePath(args[[2L]],mustWork=TRUE)
out <- args[[3L]]
source_type <- match.arg(args[[4L]],c("full","pilot"))
ids <- as.integer(strsplit(args[[5L]],",",fixed=TRUE)[[1L]])
stopifnot(length(ids)>0,!anyNA(ids),!anyDuplicated(ids))
dir.create(out, recursive=TRUE, showWarnings=FALSE)
study <- "anchoredRI/tools/subsample"
roster <- read.csv(file.path(study,"pilot-roster.csv"),stringsAsFactors=FALSE)
index <- read.csv(file.path(input,"task-g-index.csv"),stringsAsFactors=FALSE)
design <- read.csv(file.path(input,"design-large.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(index)==1440L,nrow(design)==1440L,!anyDuplicated(index$Index),
          !anyDuplicated(design$Index),setequal(index$Index,design$Index),nrow(roster)==20L,
          !anyDuplicated(roster$Index),all(ids %in% design$Index),
          all((ids %in% roster$Index) == (source_type=="pilot")))
source_files <- file.path(study,c("pilot.R","aggregate.R","inputs.R","mixture.R"))
image_digest <- "sha256:5e8963d73b87865e0ad779df4202555fb1b436ee95dc7d7442450308b8ea77f3"
readtab <- function(d, name) read.csv(file.path(d,paste0(name,".csv")),stringsAsFactors=FALSE,
                                      na.strings=c("","NA"),check.names=FALSE)
# The Task G digest uses R serialization v2 minus the standard 14-byte header.
# R 4.5 uses its locked digest library. Base R gained the same hash in 4.6.
file_hash <- if (requireNamespace("digest", quietly=TRUE))
  function(path) digest::digest(file=path, algo="sha256") else
  function(path) unname(tools::sha256sum(path))
digest_v2 <- if (requireNamespace("digest", quietly=TRUE))
  function(x) digest::digest(x, algo="sha256",serializeVersion=2L) else
  function(x) unname(tools::sha256sum(bytes=serialize(x,NULL,version=2)[-(1:14)]))
# CSVs were exported with %.17g. Numeric readback across R versions may differ by an ulp.
same <- function(x,y,eps=5e-13) {
  if (!length(x)) x <- NA
  if (!length(y)) y <- NA
  if (length(x)!=length(y) || !identical(as.vector(is.na(x)),as.vector(is.na(y)))) return(FALSE)
  i <- !is.na(x)
  if (!any(i)) return(TRUE)
  if (is.numeric(x) || is.numeric(y)) {
    left <- as.numeric(x[i]); right <- as.numeric(y[i])
    finite <- is.finite(left) & is.finite(right)
    return(all(left[!finite]==right[!finite]) &&
           all(abs(left[finite]-right[finite]) <= eps*pmax(1,abs(right[finite]))))
  }
  all(as.character(x[i]) == as.character(y[i]))
}
fields <- function(a,b,keys) all(vapply(keys,function(k) same(a[[k]],b[[k]]),logical(1)))
valid <- function(l,u,crp) is.finite(u) & (crp | (is.finite(l) & l < u))
writetab <- function(rows,path) write.csv(rows,file.path(out,path),row.names=FALSE,na="")

# The two approved pilot versions differed only in fallback timing
# instrumentation. Their source manifests are permanent on the pilot branch.
pilot_common <- c("32e0ead6fec27a443eb292a0c787d0d222f7af528bd48ad49f80190e0f113683",
                  "46035378ab8d0717fc3cc97fa6facf0596966ba94431f07f69e241d7b7504da4",
                  "c389f872173bf54693712389befcf5afe23ccfaefdaea7ed2ed962d7b735eb16")
pilot_main <- c("7ec50d0551dbf170f9830b33d4dd56ac9ce4f957bb2ae209dd2694a553f5af01",
                "72fb53503a5908d878054c2e0cd6acdeb09f40069e0ec617410545dd15beee27")
current_sha <- vapply(source_files,file_hash,"")
required <- c("COMPLETE","decisions.csv","draws.csv","fits.csv","input-id.rds",
              "medians.csv","models.rds","positions.rds","timing.csv")
hashes <- list()
records <- anomalies <- failures <- list()
context <- ""
check <- function(pass,what) {
  if (!isTRUE(pass)) stop(context," | ",what,call.=FALSE)
}
odd <- function(id,kind,count,detail="") {
  if (count) anomalies[[length(anomalies)+1L]] <<- data.frame(Index=id,kind=kind,count=count,detail=detail)
}
for (id in ids) {
  context <- paste("Index",id)
  result <- tryCatch({
    d <- file.path(records_root,as.character(id))
    check(all(file.exists(file.path(d,required))),"required record files")
    check(identical(readLines(file.path(d,"COMPLETE")),
                    sprintf("Index %d; 204 fits, 100 draws",id)),"completed record marker")
    meta <- design[match(id,design$Index),]
    i <- index[match(id,index$Index),]
    src <- readRDS(file.path(input,sprintf("taskg-shard-%d.rds",i$shard)))$cases[[i$rds_position]]
    routine <- src$routine
    anchor <- src$panels$guidance
    n <- length(routine)
    rh <- digest_v2(routine)
    ah <- digest_v2(anchor)
    crp <- meta$Analyte == "CRP"
    check(src$Index == id && src$j == id && n == i$routine_n && length(anchor)==40L &&
          rh == i$routine_hash && rh == meta$routine_hash &&
          ah == i$anchor_hash && ah == meta$anchor_hash,"archive identities")
    check(meta$N==50000 && i$N==50000 && meta$anchor_n==40L &&
            meta$routine_n==n && meta$anchor_seed==17000000L+id &&
            meta$pit_primary_seed==38000000L+id && isTRUE(meta$requested_lower == !crp) &&
            meta$upper_probability==(if(crp) .95 else .975),"design parameters and seeds")
    if(source_type=="pilot") {
      pilot_meta <- roster[match(id,roster$Index),]
      check(rh==pilot_meta$routine_hash && ah==pilot_meta$anchor_hash,"published pilot hashes")
    }
    version <- readRDS(file.path(d,"input-id.rds"))
    check(version$Index == id && version$routine_hash == rh && version$anchor_hash == ah &&
          startsWith(version$runtime,"R version 4.5.0") && version$anchoredRI=="0.4.0","runtime receipt")
    code <- version$code
    check(identical(names(code),source_files),"source files in runtime receipt")
    if(source_type=="pilot") {
      check(unname(code[1L]) %in% pilot_main && identical(unname(code[-1L]),pilot_common),
            "pinned pilot source hashes")
      # Original pilot receipts left the container field empty. Their separately
      # published image.json binds those records to the immutable image digest.
      check(version$container %in% c("", image_digest),"pilot runtime receipt")
    } else {
      check(identical(code,current_sha) && identical(version$container,image_digest),
            "full run source files and immutable container")
    }
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
      set.seed(expected_seed, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
      replay <- sort(sample.int(n,5000L,replace=FALSE))
      check(identical(p,replay),"replayed 5,000 sampled positions")
      check(is.integer(p) && length(p)==5000L && !is.unsorted(p, strictly=TRUE) &&
            all(p>=1L & p<=n) && digest_v2(p)==dr$position_hash[j] &&
            dr$sampling_seed[j]==expected_seed && dr$positions_n[j]==5000,"position order/bounds/hash/draw seed")
    }
    fit_key <- function(t) with(t,paste(cap,fit_type,seed_set,draw,sep="/"))
    k <- fit_key(f)
    check(!anyDuplicated(k) && setequal(f$cap,c(3L,5L)) && all(f$routine_hash==rh & f$anchor_hash==ah) &&
            all(f$fit_seed == 1L) && all(f$rounding_unit==meta$h),"fit keys/identity")
    check(all(f$first_status %in% c("reported","declined","withheld")) &&
            all(f$status %in% c("reported","declined","withheld")),"known fit statuses")
    check(all(is.na(f$fallback_status) | (f$cap==5L & f$fit_type=="full")),
          "only standard public full fit invokes fallback")
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
      fallback <- fit$diagnostics$anchorRefine
      reason <- function(z) {
        text <- paste(z,collapse=" | ")
        if(nzchar(text)) text else NA_character_
      }
      check(same(row$optimizer_code,fit$optimizer_code) &&
            same(row$optimizer_status,fit$optimizer_status) &&
            same(row$fit_error,reason(fit$diagnostics$anchored$error)) &&
            same(row$check_error,reason(fit$safeguard$error)) &&
            same(row$first_reason,reason(fit$anchoredRI_reason)) &&
            same(row$reason,reason(fit$reason)) && same(row$fallback_reason,reason(fallback$reason)),
            "optimizer exits, technical errors and recorded withholding reasons")
      if (row$fit_type=="full" && row$cap==5L) {
        called <- row$first_status!="reported"
        check((!is.null(fallback))==called,"public fallback invocation")
      } else check(is.null(fallback),"absence of subsample, W or cap-3 fallback")
      released_fallback <- identical(fallback$status,"reported")
      final_status <- if(released_fallback) "reported" else row$first_status
      method <- if(released_fallback) "anchorRefine" else if(final_status=="reported") "anchoredRI" else NA
      check(same(row$fallback_status,fallback$status) && row$status==final_status && same(row$method,method),
            "first-stage/fallback/final routing and method")
      if(released_fallback) check(same(fallback$limits,fit$limits),"returned fallback endpoints")
      if (length(fit$weights)) {
        choices <- fit$selection
        ord <- order(choices$bic_original,match(choices$transform,c("identity","log","shifted_log")),
                     choices$components,na.last=TRUE)
        winner <- choices[ord[1L],]
        pars <- 3L*row$components-1L+(row$transform=="shifted_log")
        check(nrow(choices)>0L && all(choices$components>=1L & choices$components<=row$cap) &&
                winner$components==row$components && winner$transform==row$transform &&
                same(winner$bic_original,row$bic) && same(fit$selected_bic,row$bic) &&
                same((pars*log(est_n)-row$bic)/2,row$loglik),"original-scale BIC winner and likelihood")
        probs <- if(crp) .95 else c(.025,.975)
        yhat <- row$mu+qnorm(probs)*row$sigma
        limits <- if(row$transform=="identity") yhat else row$shift+exp(yhat)
        check(same(limits,ft,1e-10),"limits calculated independently from selected healthy component")
        if(row$first_status=="reported") check(same(ft,fit$limits),"released first-stage endpoints")
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
      expected_policies <- c("median full","median majority","W own","W full",
                             "W majority","median majority_or_package")
      check(setequal(paste(decisions$estimator,decisions$rule),expected_policies),
            "six unique release policies")
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
    file_paths <- file.path(d,required)
    hashes[[length(hashes)+1L]] <- data.frame(Index=id,source=source_type,file=required,
      bytes=file.size(file_paths),sha256=unname(vapply(file_paths,file_hash,"")))
    rec <- data.frame(Index=id,source=source_type,Analyte=meta$Analyte,
           stratum=if(meta$fractionPathol<=.2) "<=20%" else ">20%",routine_n=n,
           routine_hash=rh,anchor_hash=ah,source_pilot_R=unname(code[1L]),container=version$container,
           share_rows=sum(is.finite(f$share)),
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
if(length(hashes)) writetab(do.call(rbind,hashes),"file-hashes.csv")
if(length(failures)) stop(length(failures)," datasets failed independent readback")
cat("Verified all",length(records),"datasets\n")
