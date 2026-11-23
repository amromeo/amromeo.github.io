# Produce the Actions input only after the entire local validation has passed.
args <- commandArgs(TRUE);stopifnot(length(args)==2)
root <- args[1];proof <- read.csv(file.path(root,'cross-machine-summary.csv'))
keys <- list.files(file.path(root,'validation'),pattern='[.]rds$')
stopifnot(file.exists(file.path(root,'VALIDATED.txt')),!file.exists(file.path(root,'HALT.txt')),
          length(keys)==5220L,nrow(proof)==1L,proof$datasets==4320L,
          proof$model_changes<=.05*proof$datasets)
proof$exact_checked <- length(keys)
write.csv(proof,args[2],row.names=FALSE)
