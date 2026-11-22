source("tools/stress-design.R")
stopifnot(nrow(stress_design)==45L,!anyDuplicated(stress_design[c("family","overlap","affected")]))
for(k in 1:45) {
  d<-stress_data(k,1L);s<-stress_design[k,]
  stopifnot(identical(d,stress_data(k,1L)),length(d$value)==5000L,length(d$anchor)==40L,
    identical(d$anchor,head(d$screened,40L)),all(table(d$setting)==2500L),
    max(abs(stress_cdf(d$truth,s$family)-c(.025,.975)))<1e-10,
    abs(d$actual_affected-s$affected)<.04)
  m<-stress_metrics(d$truth,d$truth,s$family)
  stopifnot(m$width_error==0,m$mpe==0,!m$large_error,abs(m$reference_coverage-.95)<1e-12)
  bad<-stress_metrics(c(NA,NA),d$truth,s$family)
  stopifnot(bad$status=="failed",is.na(bad$large_error))
  reversed<-stress_metrics(rev(d$truth),d$truth,s$family)
  stopifnot(reversed$status=="failed")
  wrong<-stress_metrics(d$truth+c(0,.11*diff(d$truth)),d$truth,s$family)
  stopifnot(wrong$large_error,abs(wrong$width_error-11)<1e-10)
  same_family<-which(stress_design$family==s$family)[1]
  stopifnot(identical(d$screened,stress_data(same_family,1L)$screened),
    !identical(d$screened,stress_data(k,2L)$screened))
}
stopifnot(stress_design$disease_inside_true_ri[1]<.06,
  stress_design$disease_inside_true_ri[4]>.45,stress_design$disease_inside_true_ri[7]>.8)
message("Stress design, truth, independent replicate seeds, paired samples and failure metrics passed.")
