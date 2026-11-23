args<-commandArgs(TRUE);stopifnot(length(args)==1L)
dir.create(dirname(args[1]),recursive=TRUE,showWarnings=FALSE)
stopifnot(system2(Sys.which('g++'),c('-std=c++17','-O2','-ffp-contract=off','-fPIC','-shared',
  paste0('-I',shQuote(R.home('include'))),'tools/boxcox/reliability-v3/native.cpp','-o',shQuote(args[1])))==0L)
