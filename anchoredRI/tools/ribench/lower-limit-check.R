library(RIbench)
args <- commandArgs(TRUE); w <- normalizePath(args[1]); data_dir <- normalizePath(args[2]); out <- args[3]
t <- loadTestsetDefinition(); t <- t[t$N<=5000,]
p <- read.csv(file.path(w,'anchoredRI-all-performance.csv'))
d <- do.call(rbind,lapply(list.files(w,pattern='^custom-[0-9]+[.]csv$',full.names=TRUE),read.csv));d<-d[d$algo=='anchoredRI',]
p<-merge(p,d[,c('Index','lrt_p')],by='Index');p<-merge(p,t[,c('Index','startSeed','fractionPathol')],by='Index')
p$routine_min<-vapply(seq_len(nrow(p)),function(i){s<-p[i,];x<-data.table::fread(file.path(data_dir,sprintf('Data/%s/%s_%s_seed_%s.csv',s$Analyte,s$Index,s$Analyte,s$startSeed)))$V1;min(x)},numeric(1))
p$lower_flag<-p$LRL<p$routine_min;p$lower_flag[p$Analyte=='CRP']<-FALSE
p$lrt_flag<-p$lrt_p<.01;p$either_flag<-p$lower_flag|p$lrt_flag;p$both_flag<-p$lower_flag&p$lrt_flag
p$endpoint_extreme<-abs(p$zzDev_LRL)>5|abs(p$zzDev_URL)>5;p$endpoint_extreme[p$Analyte=='CRP']<-abs(p$zzDev_URL[p$Analyte=='CRP'])>5
dir.create(out,showWarnings=FALSE)
write.csv(p,file.path(out,'lower-limit-diagnostics.csv'),row.names=FALSE)
rows<-list();groups<-list(all=rep(TRUE,nrow(p)),mean_abs_z_gt5=p$zzDevAbs_Ov>5,mean_abs_z_lt025=p$zzDevAbs_Ov<.25,either_endpoint_gt5=p$endpoint_extreme)
for(g in names(groups))for(f in c('lower_flag','lrt_flag','either_flag','both_flag')){q<-p[groups[[g]],];rows[[length(rows)+1]]<-data.frame(group=g,rule=f,n=nrow(q),flagged=sum(q[[f]]),rate=mean(q[[f]]))}
s<-do.call(rbind,rows);write.csv(s,file.path(out,'lower-limit-summary.csv'),row.names=FALSE);print(s)
