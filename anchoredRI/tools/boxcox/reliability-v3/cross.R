# Cross-test both lineages and each distinct competing winning family/K on each receiving host.
args<-commandArgs(TRUE);stopifnot(length(args)==8L,args[1] %in% c('local','Actions'))
runtime<-args[1];data_dir<-args[2];archive<-args[3];local_dir<-args[4];actions_dir<-args[5]
out<-args[6];workers<-as.integer(args[7]);source('tools/boxcox/reliability-v3/common.R');dyn.load(args[8])
sources<-r3_pedigree(runtime,archive,out)
extra<-c('tools/boxcox/reliability-v3/cross.R','tools/boxcox/reliability-v3/native-recovery.cpp',
         'tools/boxcox/reliability-v3/build-production.R')
sources<-c(sources,vapply(extra,function(f)digest::digest(file=f,algo='sha256'),character(1)))
write.csv(data.frame(file=names(sources),sha256=unname(sources)),file.path(out,'sources.csv'),row.names=FALSE)
hashes<-read.csv(file.path(archive,'anchor-sample-hashes.csv'));defs<-loadTestsetDefinition()
l<-r3_read_exchange(local_dir);a<-r3_read_exchange(actions_dir)
stopifnot(all(l$runtime=='local'),all(a$runtime=='Actions'))
pair<-merge(l[l$selected,],a[a$selected,],by='Index',suffixes=c('_local','_Actions'),sort=TRUE)
stopifnot(identical(pair$Index,1:156),all(pair$routine_hash_local==pair$routine_hash_Actions),
   all(pair$anchor_hash_local==pair$anchor_hash_Actions),all(pair$code_local=='0'),all(pair$code_Actions=='0'))
width<-pair$upper_Actions-pair$lower_Actions
finite<-is.finite(width)&width>0&is.finite(pair$upper_local)&is.finite(pair$lower_local)
match_model<-with(pair,transform_local==transform_Actions & k_local==k_Actions)
same_limits<-finite & pmax(abs(pair$upper_local-pair$upper_Actions),abs(pair$lower_local-pair$lower_Actions))<=.005*width
ids<-pair$Index[!match_model|!same_limits]
write.csv(data.frame(Index=pair$Index,model_match=match_model,limit_match=same_limits,cross_tested=pair$Index %in%ids),
          file.path(out,'membership.csv'),row.names=FALSE)
cross_one<-function(i) {
 file<-file.path(out,sprintf('cross-%03d.rds',i));s<-defs[defs$Index==i,]; input<-r3_inputs(s,data_dir,hashes);d<-input$d
 ls<-l[l$Index==i,];as<-a[a$Index==i,];stopifnot(all(c(ls$routine_hash,as$routine_hash)==input$routine),
             all(c(ls$anchor_hash,as$anchor_hash)==input$anchor))
 if(file.exists(file)){r<-readRDS(file);stopifnot(identical(r$sources,sources),r$runtime==runtime,r$Index==i);return(TRUE)}
 specs<-unique(rbind(ls[ls$selected,c('transform','k')],as[as$selected,c('transform','k')]))
 results<-list();rows<-list();traces<-list()
 for(j in seq_len(nrow(specs)))for(from in c('local','Actions')) {
  ins<-if(from=='local')ls else as
  r<-ins[ins$transform==specs$transform[j] & ins$k==specs$k[j],];seed<-r3_from_row(r);key<-paste(seed$f,seed$k,from)
  v<-r3_continue(d,seed,force=TRUE);p<-v$seed;stopifnot(p$code==0L);row<-r3_row(p,d)
  rows[[key]]<-cbind(Index=i,receiver=runtime,donor=from,input_code=seed$code,input_loglik=seed$ll,
    start_error=v$transport_error,blocks=length(v$blocks),initial_lower=r$lower,initial_upper=r$upper,row)
  traces[[key]]<-r3_traces(v$blocks,key);results[[key]]<-v
 }
 z<-list(Index=i,runtime=runtime,sources=sources,results=results,tab=do.call(rbind,rows),traces=do.call(rbind,traces))
 r3_save(z,file);cat('crossed',runtime,i,'\n');flush.console();TRUE
}
if(length(ids)) {
 r3_pool(ids,workers,cross_one);rs<-lapply(ids,function(i)readRDS(file.path(out,sprintf('cross-%03d.rds',i))))
 write.csv(do.call(rbind,lapply(rs,`[[`,'tab')),file.path(out,'outcomes.csv'),row.names=FALSE)
 write.csv(do.call(rbind,lapply(rs,function(r)cbind(Index=r$Index,receiver=runtime,r$traces))),
           file.path(out,'blocks.csv'),row.names=FALSE)
}
writeLines(sprintf('%d discrepant datasets cross-tested',length(ids)),file.path(out,'DONE.txt'))
