# Same-runtime check against main 55eee86. No archived optimizer comparison.
for(f in list.files('R',full.names=TRUE))source(f)
source('tools/stress-design.R')
frozen <- new.env()
for(f in list.files('tools/safeguard/frozen',full.names=TRUE))sys.source(f,frozen)
rows <- list()
for(case in c(1,2,3)) {
  d <- stress_data(case,21)
  old <- frozen$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')
  new <- fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto',safeguard=FALSE)
  delta <- max(abs(new$limits-old$limits))
  stopifnot(delta<=1e-10,identical(new$transform,old$transform),
    identical(new$components,old$components),identical(new$shift,old$shift),
    identical(new$selection,old$selection),identical(new$weights,old$weights))
  rows[[length(rows)+1]] <- data.frame(case=case,replicate=21,max_limit_difference=delta,
    transform=new$transform,components=new$components)
}
dir.create('tools/safeguard/validation',showWarnings=FALSE)
write.csv(do.call(rbind,rows),'tools/safeguard/validation/same-runtime-parity.csv',row.names=FALSE)
cat('Frozen-source parity passed\n')
