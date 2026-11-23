args<-commandArgs(TRUE);stopifnot(length(args)==1L);dir.create(dirname(args[1]),recursive=TRUE,showWarnings=FALSE)
base<-'tools/boxcox/reliability-v3/'
modules<-c('tools/boxcox/reliability-v2/likelihood.cpp',paste0(base,c('native.cpp','native-recovery.cpp','build-production.R')))
stopifnot(system2(Sys.which('g++'),c('-std=c++17','-O2','-ffp-contract=off','-fPIC','-shared',
   paste0('-I',shQuote(R.home('include'))),shQuote(paste0(base,'native-recovery.cpp')),'-o',shQuote(args[1])))==0L)
manifest<-vapply(modules,function(f)digest::digest(file=f,algo='sha256'),character(1))
write.csv(data.frame(file=names(manifest),sha256=unname(manifest)),paste0(args[1],'.sources.csv'),row.names=FALSE)
