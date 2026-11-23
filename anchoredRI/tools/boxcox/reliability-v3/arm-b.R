# Refit exactly the saved Arm B dataset keys. No accuracy measure is read or scored.
args<-commandArgs(TRUE);stopifnot(length(args)==5L)
cfg<-readRDS(args[1]);base_dir<-args[2];out<-args[3];workers<-as.integer(args[4])
source('tools/boxcox/reliability-v3/common.R');dyn.load(args[5]);dir.create(out,recursive=TRUE,showWarnings=FALSE)
core<-r3_sources();sources<-c(core,arm_b=digest::digest(file='tools/boxcox/reliability-v3/arm-b.R',algo='sha256'))
defs<-loadTestsetDefinition()
archives<-lapply(c(cfg$small_archive,cfg$large_archive),function(p)read.csv(file.path(p,'anchor-sample-hashes.csv')))
hashes<-do.call(rbind,archives);old_dir<-file.path(cfg$output,'fits')
files<-list.files(old_dir,'^r-[0-9]+[.]rds$'); ids<-as.integer(sub('^r-([0-9]+)[.]rds$','\\1',files))
saved<-do.call(rbind,lapply(list.files('tools/boxcox/reliability-v2/saved-audit/fit-decisions',
          '^local-armB-',full.names=TRUE),read.csv))
stopifnot(length(ids)==540L,identical(sort(paste0('r-',ids)),sort(saved$key)))
ids<-c(sort(setdiff(ids,1:156)),sort(intersect(ids,1:156)))
fit_one<-function(i) {
 out_file<-file.path(out,sprintf('fit-%04d.rds',i))
 if(file.exists(out_file)) {z<-readRDS(out_file);stopifnot(z$Index==i,identical(sources,z$sources));return(TRUE)}
 s<-defs[defs$Index==i,];stopifnot(nrow(s)==1,s$N<=50000)
 input<-r3_inputs(s,if(s$N<=5000)cfg$small_data else cfg$large_data,hashes);d<-input$d
 rb<-list()
 if(i<=156) {
  cache<-readRDS(file.path(base_dir,sprintf('fit-%03d.rds',i)))
  stopifnot(identical(cache$sources,core),cache$routine==input$routine,cache$anchor==input$anchor)
  rb<-cache$models
 }
 m<-r3_baseline(d,i,s$decimals,5L,rb);rules<-r3_rule(d,s$decimals);tab<-r3_bt(m,d,rules,input$probs)
 tab$method<-'armB';checks<-r3_nesting(m,d);before<-readRDS(file.path(old_dir,paste0('r-',i,'.rds')))$armB$candidate_table
 old<-before[before$selected,];now<-tab[tab$selected,];stopifnot(nrow(old)==1L,nrow(now)==1L)
 meta<-data.frame(Index=i,analyte=s$Analyte,stage=if(s$N<=5000)'small' else '50000',nominal_n=s$N,
   actual_n=d$n-length(d$anchor),anchor_n=length(d$anchor),rounding_unit=rules[1],
   old_transform=old$transform,old_k=old$components,new_transform=now$transform,new_k=now$k,
   changed=old$transform!=now$transform||old$components!=now$k,old_violations=saved$violations[match(paste0('r-',i),saved$key)],
   violations=sum(checks$violation),checked=sum(checks$status=='checked'),excluded=sum(checks$status=='structurally excluded'),
   unavailable=sum(checks$status=='unavailable'),min_sd_to_floor=min(tab$minimum_sd_to_floor),routine_hash=input$routine,anchor_hash=input$anchor)
 z<-list(Index=i,sources=sources,meta=meta,models=m,candidates=tab,nesting=checks,zero=r3_zero(d,m,rules))
 r3_save(z,out_file);cat('saved Arm B',i,'\n');flush.console();TRUE
}
r3_pool(ids,workers,fit_one)
rs<-lapply(sort(ids),function(i)readRDS(file.path(out,sprintf('fit-%04d.rds',i))))
write.csv(do.call(rbind,lapply(rs,`[[`,'meta')),file.path(out,'datasets.csv'),row.names=FALSE)
for(f in c('candidates','nesting','zero')) {
 z<-do.call(rbind,lapply(rs,function(r)if(!is.null(r[[f]]))cbind(Index=r$Index,r[[f]])))
 write.csv(z,file.path(out,paste0(f,'.csv')),row.names=FALSE)
}
writeLines('All 540 saved Arm B datasets completed',file.path(out,'DONE.txt'))
