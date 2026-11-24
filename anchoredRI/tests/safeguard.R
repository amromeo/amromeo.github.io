library(anchoredRI)
check <- getFromNamespace('.ri_checks','anchoredRI')
apply <- getFromNamespace('.apply_ri_safeguard','anchoredRI')
decide <- getFromNamespace('.ri_reporting_decision','anchoredRI')
raw <- function(mu,sd=1)list(limits=setNames(mu+qnorm(c(.025,.975))*sd,c('lower','upper')),
                           transform='identity',shift=0,components=1)
a <- qnorm(seq(.02,.98,length.out=40))
x <- seq(-3,3,length.out=80)
# A mismatched component is declined and never replaced with the anchor interval.
f <- apply(raw(4),x,a)
stopifnot(f$status=='declined',all(is.na(f$limits)),
  'routine data disagree with the screened sample' %in% f$reason,
  all(is.finite(f$diagnostics$anchored$limits)),all(is.finite(f$diagnostics$anchors_only$limits)),
  grepl('120',f$recommendation),grepl('20 screened',f$recommendation))
g <- apply(raw(0),seq(-1,3,length.out=80),a)
stopifnot(g$declined,g$safeguard$lower_flag,!g$safeguard$lrt_flag,
  identical(g$reason,'estimated lower limit is below every observed result'))
both <- apply(raw(-4),seq(-1,3,length.out=80),a)
stopifnot(length(both$reason)==2L,all(is.na(both$limits)))
upper <- apply(raw(0),seq(-1,3,length.out=80),a,probs=.95)
stopifnot(upper$status=='reported',length(upper$limits)==1,
  upper$safeguard$lower_check_skipped,!upper$safeguard$lower_flag,
  identical(upper$safeguard$lrt_p,g$safeguard$lrt_p))
off <- apply(raw(4),x,a,safeguard=FALSE)
stopifnot(off$status=='reported',!off$safeguard$enabled,!off$declined,
  identical(off$limits,raw(4)$limits),!length(off$reason),
  is.na(off$safeguard$lrt_p),is.na(off$safeguard$lower_flag))
large <- apply(raw(4),rep(x,length.out=50000),a)
small <- apply(raw(4),rep(x,length.out=49999),a)
stopifnot(grepl('refineR',large$recommendation),grepl('anchorRefine',large$recommendation),
          !grepl('refineR',small$recommendation))
unavailable <- apply(raw(0),x,rep(0,40))
stopifnot(unavailable$status=='withheld',unavailable$technical_withholding,all(is.na(unavailable$limits)))
# Strict inequalities and one-sided skip, including lower one-sided requests.
boundary <- check(raw(0),c(raw(0)$limits[1],0,3),a)
stopifnot(!boundary$lower_flag,!check(raw(0),x,a,probs=.05)$lower_flag)
# Public entry points preserve raw fits when research reporting is selected.
set.seed(907)
v <- rnorm(80,0,1); a <- rnorm(40,0,1); setting <- rep(c('a','b'),each=40)
f <- fit_anchored_ri(v,setting,a,safeguard=FALSE)
g <- fit_anchored_ri(v,setting,a)
stopifnot(identical(f$diagnostics$anchored,g$diagnostics$anchored),
  identical(f$selected_bic,g$selected_bic),isTRUE(g$safeguard$enabled))
one <- fit_anchored_ri(v,rep('one',80),a,probs=.95)
stopifnot(length(one$limits)==1,ncol(one$weights)==1,one$safeguard$lower_check_skipped)
rng <- .Random.seed
b <- suppressWarnings(bootstrap_anchored_ri(v,setting,a,transform='identity',n_boot=4,seed=17))
stopifnot(identical(rng,.Random.seed),nrow(b$bootstrap)==4,
  b$successful+b$failed+b$declined_draws+b$withheld_draws==4,
  all(c('lrt_p','lower_flag','diagnostic_lower','anchors_lower') %in% names(b$bootstrap)),
  length(b$diagnostics$bootstrap)==4)
# A declined point fit must not expose a confidence interval even if draws pass.
if(b$status!='reported')stopifnot(all(is.na(b$limits)),all(is.na(b$ci)))
b1 <- suppressWarnings(bootstrap_anchored_ri(v,rep('one',80),a,transform='identity',
  n_boot=4,seed=18,probs=.95,safeguard=FALSE))
stopifnot(identical(dim(b1$ci),c(1L,2L)),!b1$ci_conditional_on_safeguard,
  all(b1$bootstrap$lower_check_skipped),all(is.na(b1$bootstrap$lrt_p)),b1$declined_draws==0)
# Sequential and parallel safeguard decisions must agree.
if(.Platform$OS.type!='windows') {
  b2 <- suppressWarnings(bootstrap_anchored_ri(v,setting,a,transform='identity',n_boot=4,seed=17,workers=2))
  stopifnot(identical(b$bootstrap,b2$bootstrap),identical(b$ci,b2$ci))
}
# Force a lower-limit decline and verify that bootstrap cannot restore it.
blocked <- suppressWarnings(bootstrap_anchored_ri(seq(10,11,length.out=40),rep('one',40),
  seq(-.1,.1,length.out=40),transform='identity',n_boot=2,seed=19))
stopifnot(blocked$declined,all(is.na(blocked$limits)),all(is.na(blocked$ci)),
  all(is.finite(blocked$diagnostics$point$anchored$limits)))
