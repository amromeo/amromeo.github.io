source('tools/anchors-only/estimators.R')
a<-seq(12,30,length.out=40)
f<-ao_oracle(a,'normal');stopifnot(max(abs(f$limits-(mean(a)+c(-1,1)*1.96*sd(a))))<1e-12)
f<-ao_oracle(a,'shifted');stopifnot(max(abs(f$limits-(10+exp(mean(log(a-10))+c(-1,1)*1.96*sd(log(a-10))))))<1e-12)
f<-ao_auto(a);d<-f$diagnostics;stopifnot(nrow(d)==3L,all(d$status=='ok'))
# Independent original-scale lognormal likelihood and BIC.
y<-log(a);sig<-sqrt(sum((y-mean(y))^2)/40)
ll<-sum(dlnorm(a,mean(y),sig,log=TRUE))
stopifnot(abs(d$bic[d$transform=='log']-(-2*ll+2*log(40)))<1e-10)
c<-d$shift[d$transform=='shifted_log'];y<-log(a-c);sig<-sqrt(mean((y-mean(y))^2))
ll<-sum(dlnorm(a-c,mean(y),sig,log=TRUE))
stopifnot(abs(d$bic[d$transform=='shifted_log']-(-2*ll+3*log(40)))<1e-10)
b<-min(a)+mad(a)*c(-5,-.001);stopifnot(c>=b[1],c<=b[2])
# Location equivariance of profile optimum and resulting interval.
g<-ao_auto(a+100)$diagnostics
stopifnot(abs(g$shift[3]-d$shift[3]-100)<1e-4,abs(g$bic[3]-d$bic[3])<1e-6)
f<-ao_auto(rep(2,40));stopifnot(all(is.na(f$limits)),all(f$diagnostics$status=='failed'))
f<-ao_auto(seq(-5,5,length.out=40));stopifnot(f$diagnostics$status[2]=='failed',all(is.finite(f$limits)))
cat('Estimator, Jacobian, parameter penalty, profile bounds and failure tests passed.\n')
