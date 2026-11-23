source('tools/boxcox/models.R')
source('tools/stress-design.R')
rows <- list()
for(G in 1:2) for(k in 1:3) for(lambda in c(-1,-.000001,0,.000001,.432,1,2)) {
  x <- c(2,3,5,8); a <- c(2.1,3.1,4.1)
  counts <- matrix(seq_len(4*G),4,G)
  par <- c(seq(.1,.5,length.out=k),rep(log(1.3),k),rep(.1,G*(k-1)),lambda,.2)
  fn <- bc_objective(x,a,counts,k,center=1,scale=3)
  analytic <- fn(par,TRUE)
  numeric <- vapply(seq_along(par),function(j) {
    plus <- minus <- par; h <- 1e-6
    plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
    (fn(plus)-fn(minus))/(2*h)
  },numeric(1))
  delta <- max(abs(analytic-numeric)/(1+abs(numeric)))
  stopifnot(delta<2e-6)
  rows[[length(rows)+1]] <- data.frame(G=G,k=k,lambda=lambda,max_relative_gradient_error=delta)
}
# Compression and Jacobians must agree with the expanded observations.
x <- c(1.2,2,4); a <- c(1.1,2.2); cnt <- matrix(c(3,2,4),3,1)
p <- c(.5,log(1.1),.432,.1)
f1 <- bc_objective(x,a,cnt,1)
f2 <- bc_objective(rep(x,cnt),a,matrix(1,sum(cnt),1),1)
stopifnot(abs(f1(p)-f2(p))<1e-12)
for(lambda in c(-2,-1,0,.432,1,2)) {
  x <- c(.1,1,2,5)
  stopifnot(max(abs(bc_inverse(bc_transform(x,lambda),lambda)-x))<1e-10)
}
stopifnot(is.na(bc_inverse(-2,1)))
for(case in c(1,2,3)) {
  d <- stress_data(case,1)
  frozen <- bc_frozen(FALSE); cap <- bc_capacity(3)
  f <- frozen$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')
  b <- cap$fit_anchored_ri(d$value,d$setting,d$anchor,transform='auto')
  stopifnot(identical(f$limits,b$limits),identical(f$selection,b$selection))
  for(lambda in c(0,1)) {
    a <- bc_fixed(d$value,d$setting,d$anchor,lambda)
    z <- frozen$fit_anchored_ri(d$value,d$setting,d$anchor,transform=if(lambda==1)'identity' else 'log')
    stopifnot(identical(a$candidate_limits,z$candidate_limits),identical(a$bic_original,z$bic_original))
  }
}
dir.create('tools/boxcox/validation',showWarnings=FALSE)
write.csv(do.call(rbind,rows),'tools/boxcox/validation/gradient-checks.csv',row.names=FALSE)
writeLines('Gradient, multiplicity, inverse-domain, frozen fit and fixed-shape objective checks passed.',
           'tools/boxcox/validation/unit-checks.txt')
cat('All model checks passed\n')
