args <- commandArgs(TRUE); cfg <- readRDS(args[1]); root <- args[2]
out <- file.path(root,"tables"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
source("tools/boxcox/reliability-v3/common.R")
source("tools/boxcox/post-round3/comparator.R")
source("tools/stress-design.R")
defs <- loadTestsetDefinition(); defs <- defs[defs$N<=5000,]
hashes <- read.csv(file.path(cfg$small_archive,"anchor-sample-hashes.csv"))
custom <- do.call(rbind,lapply(list.files(cfg$small_archive,"^custom-[0-9]+[.]csv$",full.names=TRUE),read.csv))
custom <- custom[custom$algo=="anchoredRI",]
arc <- read.csv(file.path(cfg$small_archive,"anchoredRI-all-performance.csv"))
archive_grid <- read.csv(unz("assets/disagreement/disagreement-results.zip","datasets.csv"))
read_arch <- function(zip,member) { con<-unz(zip,member,open="r");on.exit(close(con));read.csv(con) }
legacy <- list(anchorsOnly=read_arch("assets/anchors-only/anchors-only-results.zip","estimates.csv"),
               refineR=read_arch("assets/comparison-study/summary/comparison-results.zip","all-results.csv"))
legacy$anchorsOnly <- subset(legacy$anchorsOnly,method=="anchors_auto")
legacy$refineR <- subset(legacy$refineR,method=="refineR")
old_csv <- "tools/boxcox/anchorRefine-comparator/anchoredRI-comparator.csv"
files <- list.files(file.path(root,"fits"),"[.]rds$",full.names=TRUE)
partial <- length(args)>=3L && args[3]=="preview"
if(!partial)stopifnot(file.exists(file.path(root,"FITTED.txt")),length(files)==3780L)
parameters <- diagnostics <- exit_details <- timings <- candidates <- simulations <- baselines <- list()
for (j in seq_along(files)) {
  z <- readRDS(files[j]);key <- z$key;fields <- strsplit(key,"-",fixed=TRUE)[[1]]
  small <- fields[1]=="r"
  if(small) {
    i <- as.integer(fields[2]);s <- defs[defs$Index==i,]
    d <- bc_ribench_input(s,cfg$small_data,hashes$sha256[match(i,hashes$Index)])
    meta <- data.frame(schema_version=1L,dataset_id=key,benchmark="RIbench",stage="small",Index=i,
          case=NA_integer_,replicate=NA_integer_,Analyte=s$Analyte,N=s$N,n_routine=length(d$value),
          n_anchor=length(d$anchor),distribution=s$Distribution,pathological_fraction=s$fractionPathol,
          scenario=s$Scenario,rounding_step=10^(-s$decimals),routine_seed=s$startSeed,anchor_seed=17000000L+i,
          routine_hash=z$routine,anchor_hash=z$anchor,ribench_hashcode=s$HashCode)
  } else {
    d <- stress_data(as.integer(fields[2]),as.integer(fields[3]))
    d$probs <- c(.025,.975);s <- stress_design[d$case,]
    meta <- data.frame(schema_version=1L,dataset_id=key,benchmark="simulation",stage="grid",Index=NA_integer_,
         case=d$case,replicate=d$replicate,Analyte=NA_character_,N=length(d$value),n_routine=length(d$value),
         n_anchor=length(d$anchor),distribution=d$family,pathological_fraction=s$affected,scenario=s$overlap,
         rounding_step=.01,routine_seed=d$seed,anchor_seed=9300000L+10000L*match(d$family,c("normal","lognormal","shifted"))+d$replicate,
         routine_hash=z$routine,anchor_hash=z$anchor,ribench_hashcode=NA_character_)
  }
  stopifnot(z$routine==digest::digest(d$value,algo="sha256"),z$anchor==digest::digest(d$anchor,algo="sha256"))
  full <- v2_table(z$models,d$probs)
  pair <- lapply(c(3L,5L),function(k) comparator_row(z$models,meta,d$value,d$anchor,d$probs,k,
                                                if(k==3) "anchoredRI" else "armB"))
  baselines[[j]]<-pair[[1]]
  for(a in pair) {
    m<-z$models[[paste(a$transform,a$components)]]
    par<-m$p;settings<-unique(as.character(d$setting))
    stopifnot(nrow(par$w)==m$k,ncol(par$w)==length(settings))
    ix<-expand.grid(component=seq_len(m$k),setting_index=seq_along(settings))
    parameters[[length(parameters)+1L]]<-data.frame(dataset_id=key,stage=meta$stage,method=a$method,
       transform=m$transform,lambda=par$lambda,shift=par$shift,component=ix$component,
       screened_component=ix$component==1L,setting=settings[ix$setting_index],
       transformed_mean=par$mu[ix$component]+if(m$transform=="identity")1 else 0,
       transformed_sd=par$sigma[ix$component],routine_weight=par$w[as.matrix(ix)])
  }
  timings[[j]]<-data.frame(dataset_id=key,stage=meta$stage,analyte=meta$Analyte,N=meta$N,
       unique_values=length(unique(c(d$value,d$anchor))),cached=z$cached,began=z$began,ended=z$ended,cpu=z$cpu,fit_seconds=z$elapsed)
  full$selected_baseline<-full$transform==pair[[1]]$transform & full$k==pair[[1]]$components
  candidates[[j]]<-cbind(dataset_id=key,full)
  old <- readRDS(file.path(cfg$output,"validation",paste0(key,".rds")))
  record <- function(method,lower,upper,p=NA_real_,lower_flag=FALSE,declined=FALSE,technical=FALSE,
                     transform=NA_character_,components=NA_integer_,origin="",code=NA_integer_,shift=NA_real_,bound=NA) {
    data.frame(dataset_id=key,stage=meta$stage,Index=meta$Index,case=meta$case,replicate=meta$replicate,
      distribution=meta$distribution,fraction=meta$pathological_fraction,scenario=meta$scenario,Analyte=meta$Analyte,
      method=method,lower=lower,upper=upper,lrt_p=p,lower_flag=lower_flag,declined=declined,technical=technical,
      reported=!(declined||technical),transform=transform,components=components,origin=origin,code=code,shift=shift,at_bound=bound)
  }
  rr <- lapply(pair,function(a){
    m <- z$models[[paste(a$transform,a$components)]]
    st <- m$starts[[m$origin]]
    record(a$method,a$lower,a$upper,a$lrt_p,a$lower_flag,a$declined,a$technical,a$transform,a$components,
           st$label,m$code,a$shift,m$at_bound)
  })
  g<-old$safeguard;f<-old$baseline
  rr[[3]]<-record("frozenLocal",if(length(g$limits)==2)g$limits[1] else NA_real_,tail(g$limits,1),
                  g$lrt_p,g$lower_flag,g$declined,g$technical,f$transform,f$components,shift=f$shift,bound=f$at_bound)
  if(small) {
    a<-arc[match(i,arc$Index),];c2<-custom[match(i,custom$Index),]
    lo<-if(s$Analyte=="CRP")NA_real_ else a$LRL;hi<-a$URL
    p<-c2$lrt_p;model<-c2$model;k<-c2$components
  } else {
    a<-archive_grid[archive_grid$case==d$case & archive_grid$replicate==d$replicate,]
    stopifnot(nrow(a)==1)
    lo<-a$diagnostic_anchored_lower;hi<-a$diagnostic_anchored_upper;p<-a$lrt_p;model<-NA_character_;k<-NA_integer_
  }
  lf<-length(d$probs)==2 && is.finite(lo) && lo<min(d$value)
  tech<-!is.finite(p)||!is.finite(hi)||(length(d$probs)==2 && !is.finite(lo))
  rr[[4]]<-record(if(small)"archivedC2" else "archivedGrid",lo,hi,p,lf,!tech&&(p<.01||lf),tech,model,k)
  if(!small)for(method in names(legacy)) {
    a<-legacy[[method]];a<-a[a$case==d$case & a$replicate==d$replicate,]
    stopifnot(nrow(a)==1L,a$seed==d$seed)
    rr[[length(rr)+1L]]<-record(method,a$lower,a$upper,technical=!stress_valid(c(a$lower,a$upper)))
  }
  rows<-do.call(rbind,rr);diagnostics[[j]]<-rows
  for(a in pair) {
    m<-z$models[[paste(a$transform,a$components)]];winning<-m$starts[[m$origin]]
    line_search<-grepl("line[ -]?search|lnsrch",winning$message,ignore.case=TRUE)
    iteration_limit<-a$convergence_code==1L || grepl("iterat.*limit|limit.*iterat",winning$message,ignore.case=TRUE)
    if(!iteration_limit && !line_search && a$convergence_code!=52L)next
    exit_details[[length(exit_details)+1L]]<-data.frame(dataset_id=key,stage=meta$stage,method=a$method,
       transform=m$transform,components=m$k,code=m$code,origin=winning$label,message=winning$message,iteration_limit=iteration_limit,line_search=line_search,
       initial_retained=m$initial_retained,initial_ll=winning$initial,final_ll=winning$final,
       retained_ll=m$ll,declined=a$declined,reported=a$report_status=="reported")
  }
  if(!small) {
    gridrows<-list()
    for(k in seq_len(nrow(rows))) for(post in c(FALSE,TRUE)) {
      row<-rows[k,];lim<-c(row$lower,row$upper)
      if(post && !row$reported)lim[]<-NA_real_
      met<-suppressWarnings(stress_row(d,row$method,lim))
      met$n_screened<-if(row$method=="refineR")0L else length(d$anchor)
      gridrows[[length(gridrows)+1L]]<-cbind(dataset_id=key,filter=if(post)"after"else"before",
                 declined=row$declined,technical=row$technical,reported=row$reported,met)
    }
    simulations[[j]]<-do.call(rbind,gridrows)
  }
}
bind <- function(x)do.call(rbind,x)
save_outputs<-function(x,n)write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE,na="")
diag<-bind(diagnostics);save_outputs(diag,"diagnostics")
exits<-if(length(exit_details))bind(exit_details)else data.frame(dataset_id=character(),stage=character(),method=character(),code=integer())
save_outputs(exits,"winning-optimizer-exits")
save_outputs(bind(parameters),"selected-component-parameters");save_outputs(bind(timings),"timings");save_outputs(bind(candidates),"candidates")
save_outputs(bind(simulations),"grid-outcomes")
write_comparator(baselines,file.path(out,"stage-comparator.csv"))
if(!partial) {
    existing<-read.csv(old_csv,stringsAsFactors=FALSE,na.strings="")
    new<-read.csv(file.path(out,"stage-comparator.csv"),stringsAsFactors=FALSE,na.strings="")
    ix<-match(existing$dataset_id,new$dataset_id); keep<-existing$stage=="small"
    stopifnot(sum(keep) %in% c(396L,2880L),all(!is.na(ix[keep])),identical(names(new),names(existing)))
    for(col in c("routine_hash","anchor_hash","lower","upper","lrt_p","components","transform"))
      stopifnot(isTRUE(all.equal(existing[[col]][keep],new[[col]][ix[keep]],tolerance=0)))
    # Keep the existing 144 large-N rows byte-for-byte. No large-N fit is performed here.
    old_lines<-readLines(old_csv);new_lines<-readLines(file.path(out,"stage-comparator.csv"))
    writeLines(c(new_lines,old_lines[-1][existing$stage=="50000"]),file.path(out,"anchoredRI-comparator.csv"))
}
