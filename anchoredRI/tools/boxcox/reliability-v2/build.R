args <- commandArgs(TRUE)
stopifnot(length(args)==1L)
src <- 'tools/boxcox/reliability-v2/likelihood.cpp'
dir.create(dirname(args[1]),recursive=TRUE,showWarnings=FALSE)
compiler <- Sys.which('g++')
stopifnot(nzchar(compiler),system2(compiler,c('-std=c++17','-O2','-ffp-contract=off','-fPIC','-shared',
  paste0('-I',shQuote(R.home('include'))),shQuote(src),'-o',shQuote(args[1])))==0L)
