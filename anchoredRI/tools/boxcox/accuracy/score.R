args<-commandArgs(TRUE);cfg<-readRDS(args[1]);root<-args[2]; out<-file.path(root,"tables")
library(RIbench)
source("tools/boxcox/native-scores.R")
save_tab<-function(d,name)write.csv(d,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
diag<-read.csv(file.path(out,"diagnostics.csv"),stringsAsFactors=FALSE)
large <- length(args) >= 3L && args[3] == "50000"
stage <- if (large) "50000" else "small"
archive <- if (large) cfg$large_archive else cfg$small_archive
defs<-loadTestsetDefinition()
defs<-defs[defs$Index %in% diag$Index[diag$stage==stage],]
preview<-length(args)>=3L && args[3]=="preview"
if(large) {
  stopifnot(nrow(defs)==1440L,all(defs$N==50000L),!anyDuplicated(diag[c("Index","method")]),
            nrow(diag)==4L*nrow(defs))
} else if(!preview) stopifnot(nrow(defs)==2880L,nrow(diag[diag$stage=="grid",])==5400L)
methods<-c("anchoredRI","armB",if(large)"frozenActions" else "frozenLocal","archivedC2","anchorsOnly","refineR")
read_native<-function(m)readRDS(file.path(archive,paste0(m,"-native-performance.rds")))$resRIs
template<-read_native("anchoredRI")
markers<-unique(defs$Analyte)
native<-scores<-categories<-list()
for(m in methods) {
  archived<-m %in% c("archivedC2","anchorsOnly","refineR")
  source_results<-if(archived && m!="archivedC2")read_native(m)else template
  d<-diag[diag$stage==stage & diag$method==m,]
  ris<-lapply(markers,function(analyte){
    a<-source_results[[analyte]]
    s<-defs[defs$Analyte==analyte,]
    a<-a[match(s$Index,a$Index),seq_len(match("Error_LRL",names(a))-1L),drop=FALSE]
    stopifnot(identical(as.integer(s$Index),as.integer(a$Index)))
    a$Algorithm<-m
    if(!archived) {
      q<-d[match(a$Index,d$Index),]
      stopifnot(all(q$Index==a$Index))
      a$LRL<-if(analyte=="CRP")0 else q$lower;a$URL<-q$upper;a$Runtime<-NA_real_
      for(k in grep("^Est_",names(a)))a[[k]]<-NA_real_
    }
    a
  })
  native[[m]]<-RIbench:::computePerfMeasAll(as.character(markers),m,ris,defs,cutoffZ=5)
  p<-do.call(rbind,native[[m]]$resRIs)
  if(archived) {
    old<-do.call(rbind,source_results)
    for(k in c("zzDevAbs_Ov","zzDevAbsCutoff_Ov","zzDev_LRL","zzDev_URL"))
      stopifnot(isTRUE(all.equal(p[[k]],old[[k]][match(p$Index,old$Index)],tolerance=1e-10)))
  }
  s<-defs[match(p$Index,defs$Index),]
  if(nrow(d)) {
     q<-d[match(p$Index,d$Index),];stopifnot(all(q$Index==p$Index))
     p$declined<-q$declined;p$technical<-q$technical
     p$lrt_p<-q$lrt_p;p$lower_flag<-q$lower_flag
     p$transform<-q$transform;p$components<-q$components;p$at_bound<-q$at_bound
  } else {
     p$declined<-FALSE;p$technical<-!is.finite(p$URL)|(p$Analyte!="CRP"&!is.finite(p$LRL))
     p$lrt_p<-NA_real_;p$lower_flag<-FALSE;p$transform<-NA_character_;p$components<-NA_integer_;p$at_bound<-NA
  }
  p$reported<-!p$declined & !p$technical;p$distribution<-s$Distribution;p$fraction<-s$fractionPathol
  p$stratum<-paste(s$Distribution,ifelse(s$fractionPathol<=.2,"low <=20%","high >20%"));p$method<-m
  for(view in c("before","after")) {
    keep<-if(view=="before")rep(TRUE,nrow(p))else p$reported
    pkg<-list(resRIs=split(p[keep,],p$Analyte[keep]))
    cats<-RIbench:::getBenchmarkResults(errorList=list(pkg),nameVec=m,tableTCs=defs[defs$Index %in% p$Index[keep],],cutoffZ=5,perfCombination="mean")
    names(cats)[3]<-"value";categories[[paste(m,view)]]<-cbind(method=m,filter=view,cats)
  }
  scores[[m]]<-p
}
save_tab(do.call(rbind,categories),"ribench-category-scores")
cols<-unique(unlist(lapply(scores,names)))
all<-do.call(rbind,lapply(scores,function(x){for(col in setdiff(cols,names(x)))x[[col]]<-NA;x[,cols]}))
save_tab(all,"ribench-outcomes")
bad<-(all$Analyte!="CRP" & abs(all$zzDev_LRL)>5) | abs(all$zzDev_URL)>5
save_tab(all[which(bad),],"ribench-extreme-endpoints")
save_tab(all[which(all$zzDevAbs_Ov>5),],"ribench-implausible")
# Withheld intervals have no post-safeguard endpoints. Raw diagnostics above retain both endpoints.
reported<-all; hide<-!reported$reported
reported[hide,c("LRL","URL",grep("^(Error_|PercError_|AbsError|AbsPercError|zz)",names(reported),value=TRUE))]<-NA
save_tab(reported,"ribench-after-safeguard")
groups<-function(p)c(list(all=rep(TRUE,nrow(p)),`lambda=0.432`=abs(p$GT_Lambda-.432)<1e-12),
 setNames(lapply(unique(p$distribution),function(s)p$distribution==s),unique(p$distribution)),
 setNames(lapply(unique(p$stratum),function(s)p$stratum==s),unique(p$stratum)))
rows<-paired<-flagrows<-list()
for(m in methods) {
 p<-scores[[m]]
 for(g in names(groups(p)))for(view in c("before","after")) {
   q<-p[groups(p)[[g]],];r<-if(view=="before")q else q[q$reported,]
   e<-r$zzDevAbs_Ov;ok<-is.finite(e);cs<-bc_cutoff_scores(e,r$zzDevAbsCutoff_Ov)
   rows[[length(rows)+1L]]<-data.frame(method=m,group=g,filter=view,attempts=nrow(q),
      reported=sum(q$reported),declined=sum(q$declined),technical=sum(q$technical),decline_rate=mean(q$declined),
      missing=sum(!ok),scored=sum(ok),mean_abs_z=mean(e[ok]),capped_mean=cs$native,capped_n=cs$native_n,
      all_attempt_capped=cs$all_attempt,implausible=sum(e>5,na.rm=TRUE),implausible_or_missing=sum(!ok|e>5))
 }
 for(group in c("all","below_0.25","above_5","unscored")){
   e<-p$zzDevAbs_Ov;keep<-switch(group,all=rep(TRUE,nrow(p)),below_0.25=is.finite(e)&e<.25,
                                 above_5=is.finite(e)&e>5,unscored=!is.finite(e))
   q<-p[keep,]
   flagrows[[length(flagrows)+1L]]<-data.frame(method=m,group=group,attempts=nrow(q),
       lrt_flags=sum(q$lrt_p<.01,na.rm=TRUE),lower_flags=sum(q$lower_flag,na.rm=TRUE),
       declined=sum(q$declined),technical=sum(q$technical),reported=sum(q$reported))
 }
}
for(m in c("anchoredRI","armB"))for(other in setdiff(methods,m)) {
 p<-scores[[m]];r<-scores[[other]][match(p$Index,scores[[other]]$Index),]
 stopifnot(all(p$Index==r$Index))
 for(g in names(groups(p)))for(view in c("before","after")){
   mask<-groups(p)[[g]]; keep<-if(view=="before") mask else mask&p$reported&r$reported
   delta<-p$zzDevAbs_Ov[keep]-r$zzDevAbs_Ov[keep];good<-is.finite(delta);delta<-delta[good]
   bothcut<-keep&is.finite(p$zzDevAbsCutoff_Ov)&is.finite(r$zzDevAbsCutoff_Ov)
   paired[[length(paired)+1L]]<-data.frame(arm=m,comparator=other,stratum=g,filter=view,attempts=sum(mask),
       retained=sum(keep),complete=sum(good),mean_difference=mean(delta),wins=sum(delta<0),ties=sum(delta==0),
       losses=sum(delta>0),both_cutoff=sum(bothcut),paired_cutoff_difference=mean(p$zzDevAbsCutoff_Ov[bothcut]-r$zzDevAbsCutoff_Ov[bothcut]))
 }
}
summ<-do.call(rbind,rows);save_tab(summ,"ribench-summary")
save_tab(do.call(rbind,paired),"ribench-pairs");save_tab(do.call(rbind,flagrows),"ribench-safeguard-by-error")
if (large) {
  source("tools/boxcox/accuracy/large/criteria.R", local = TRUE)
  writeLines("Large-N scoring complete", file.path(root, "SCORED.txt"))
  quit(save = "no", status = 0L)
}
before<-subset(summ,group=="all"&filter=="before");thr<-before$capped_mean[before$method=="frozenLocal"]
criteria<-do.call(rbind,lapply(c("anchoredRI","armB"),function(m){
 a<-before[before$method==m,];data.frame(method=m,small_native_capped=a$capped_mean,threshold=thr,
    all_attempt_capped=a$all_attempt_capped,missing=a$missing,passed=bc_error_criterion(a,thr))
}));save_tab(criteria,"small-success-criteria")
grid<-read.csv(file.path(out,"grid-outcomes.csv"));gr<-gp<-list()
for(m in unique(grid$method))for(view in c("before","after")) {
 p<-subset(grid,method==m & filter==view)
 for(g in c("all",unique(p$family),as.character(sort(unique(p$case))))) {
   q<-if(g=="all")p else p[p$family==g|as.character(p$case)==g,]
   e<-q$width_error;gr[[length(gr)+1L]]<-data.frame(method=m,group=g,filter=view,attempts=nrow(q),
      reported=sum(q$reported),declined=sum(q$declined),technical=sum(q$technical),decline_rate=mean(q$declined),
      scored=sum(is.finite(e)),mean_width_error=mean(e,na.rm=TRUE),median_width_error=median(e,na.rm=TRUE),
      mean_mpe=mean(q$mpe,na.rm=TRUE),mean_coverage=mean(q$reference_coverage,na.rm=TRUE),
      errors_above_10=sum(e>10,na.rm=TRUE))
 }
}
for(m in c("anchoredRI","armB"))for(other in setdiff(unique(grid$method),m))for(view in c("before","after")) {
 p<-subset(grid,method==m&filter==view);r<-subset(grid,method==other&filter==view)
 r<-r[match(p$dataset_id,r$dataset_id),];stopifnot(all(p$dataset_id==r$dataset_id))
 for(g in c("all",as.character(sort(unique(p$case))))) {
    keep<-if(g=="all")rep(TRUE,nrow(p))else as.character(p$case)==g
    delta<-p$width_error[keep]-r$width_error[keep];delta<-delta[is.finite(delta)]
    case<-if(g=="all")NA_integer_ else as.integer(g);ix<-which(keep)[1]
    mc<-if(!is.na(case)&&length(delta)==20L)sd(delta)/sqrt(20) else NA_real_
    gp[[length(gp)+1L]]<-data.frame(arm=m,comparator=other,filter=view,group=g,case=case,
     family=if(g!="all")p$family[ix]else NA,overlap=if(g!="all")p$overlap[ix]else NA,
     affected=if(g!="all")p$affected[ix]else NA,attempts=sum(keep),complete=length(delta),
     mean_difference=mean(delta),MCSE=mc,wins=sum(delta<0),ties=sum(delta==0),losses=sum(delta>0),
     above_two_mcse=mean(delta)>2*mc)
 }
}
save_tab(do.call(rbind,gr),"grid-summary")
grid_pairs <- do.call(rbind, gp)
save_tab(grid_pairs, "grid-pairs")
reviewed <- subset(grid_pairs, arm == "anchoredRI" & comparator == "frozenLocal" &
                     filter == "before" & above_two_mcse == TRUE & case %in% c(2L, 19L, 20L) & abs(mean_difference) < .001)
reviewed <- reviewed[order(reviewed$case), c("case", "family", "overlap", "affected", "complete", "mean_difference", "MCSE")]
reviewed$arms <- rep("improved baseline and Arm B", nrow(reviewed))
reviewed$review <- rep("numerically negligible", nrow(reviewed))
save_tab(reviewed, "grid-review")
if(!preview)writeLines("Small-N and grid complete",file.path(root,"SCORED.txt"))
