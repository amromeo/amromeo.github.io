# Phase 1. Baseline and continued Arm A winners, with no reference-truth scoring.
args<-commandArgs(TRUE);stopifnot(length(args)==7L,args[1] %in% c('local','Actions'))
runtime<-args[1];data_dir<-args[2];archive<-args[3];prior<-args[4];out<-args[5];workers<-as.integer(args[6])
source('tools/boxcox/reliability-v3/common.R');dyn.load(args[7])
sources<-r3_pedigree(runtime,archive,out); hashes<-read.csv(file.path(archive,'anchor-sample-hashes.csv'))
defs<-loadTestsetDefinition();dir.create(file.path(out,'exchange'),showWarnings=FALSE)
fit_one<-function(i) {
 file<-file.path(out,sprintf('fit-%03d.rds',i))
 s<-defs[defs$Index==i,];stopifnot(nrow(s)==1L,s$Analyte=='Hb')
 input<-r3_inputs(s,data_dir,hashes);d<-input$d
 if(file.exists(file)) {
  z<-readRDS(file);stopifnot(identical(z$sources,sources),z$routine==input$routine,z$anchor==input$anchor,z$runtime==runtime);return(TRUE)
 }
 old<-readRDS(file.path(prior,sprintf('fit-%03d.rds',i)))
 stopifnot(old$Index==i,old$routine_hash==input$routine,old$anchor_hash==input$anchor,old$runtime==runtime)
 prev<-vapply(names(old$sources),function(f)digest::digest(file=f,algo='sha256'),character(1))
 stopifnot(identical(prev,old$sources))
 base<-r3_baseline(d,i,s$decimals);rule<-r3_rule(d,s$decimals)
 arm<-Filter(function(m)v2_shape(m$transform),old$models);seeds<-lapply(arm,r3_seed)
 before<-r3_winner(seeds,d);sel<-old$table[old$table$method=='armA' & old$table$selected,]
 stopifnot(before==paste(sel$transform,sel$k))
 traces<-list(); initial<-r3_row(seeds[[before]],d)
 repeat {
  winner<-r3_winner(seeds,d)
  if(seeds[[winner]]$code==0L)break
  stopifnot(seeds[[winner]]$code %in% c(1L,52L))
  m<-arm[[winner]]; last<-m$starts[[m$origin]]; stopifnot(!m$initial_retained)
  seed<-seeds[[winner]];seed$par<-last$final_par;seed$ll<-last$final
  continued<-r3_continue(d,seed);seeds[[winner]]<-continued$seed
  traces[[winner]]<-r3_traces(continued$blocks,winner)
 }
 arm_tab<-do.call(rbind,lapply(seeds,r3_row,d=d));winner<-r3_winner(seeds,d)
 r<-list(Index=i,runtime=runtime,routine=input$routine,anchor=input$anchor,sources=sources,models=base,
         base_tab=r3_bt(base,d,rule),zero=r3_zero(d,base,rule),nesting=r3_nesting(base,d),
         seeds=seeds,arm_tab=arm_tab,winner=winner,old_winner=before,initial=initial,traces=traces)
 r3_save(r,file);cat('saved',runtime,i,'\n');flush.console();TRUE
}
r3_pool(1:156,workers,fit_one)
rs<-lapply(1:156,function(i)readRDS(file.path(out,sprintf('fit-%03d.rds',i))))
bind<-function(field)do.call(rbind,lapply(rs,function(r)cbind(Index=r$Index,runtime=r$runtime,r[[field]])))
for(f in c('base_tab','zero','nesting','arm_tab','initial'))write.csv(bind(f),file.path(out,paste0(f,'.csv')),row.names=FALSE)
blocks<-do.call(rbind,lapply(rs,function(r)if(length(r$traces))cbind(Index=r$Index,runtime,do.call(rbind,r$traces))))
write.csv(blocks,file.path(out,'continuation-blocks.csv'),row.names=FALSE)
x<-do.call(rbind,lapply(rs,r3_exchange))
for(f in c('boxcox_shifted','boxcox_unshifted'))for(k in 1:3) {
 tab<-x[x$transform==f & x$k==k,];stopifnot(nrow(tab)==156)
 write.csv(tab,file.path(out,'exchange',sprintf('%s-K%d.csv',f,k)),row.names=FALSE)
}
writeLines('156 baseline replays and 156 Arm A winners complete',file.path(out,'DONE.txt'))
