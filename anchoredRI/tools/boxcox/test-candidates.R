source('tools/boxcox/models.R');source('tools/boxcox/candidates.R')
# Two converged Arm A candidates can have the same component count but differ in support.
make <- function(bic,shifted)list(k=1L,shifted=shifted,bic=bic,best=list(at_bound=FALSE,
  parameters=list(mu=0,sigma=.1,shift=if(shifted).5 else 0,lambda=0)))
candidates <- c(list(make(100,FALSE)),lapply(2:3,function(k)list(k=k,shifted=FALSE,bic=Inf)),
                list(make(101.5,TRUE)),lapply(2:3,function(k)list(k=k,shifted=TRUE,bic=Inf)))
f <- list(transform='boxcox',shift=0,lambda=0,components=1L,candidates=candidates,selected_bic=100)
t <- bc_collect_candidates(f,'armA',c(.025,.975),TRUE,'test');g <- bc_runner_up(t)
stopifnot(nrow(t)==6,sum(t$status=='failed')==4,g$within_2,g$gap==1.5,
          g$first_transform=='boxcox_unshifted',g$second_transform=='boxcox_shifted',
          abs(g$lower_difference-.5)<1e-12,abs(g$upper_difference-.5)<1e-12)
t1 <- bc_collect_candidates(f,'armA',.95,TRUE,'test');g1 <- bc_runner_up(t1)
stopifnot(is.na(g1$lower_first),g1$upper_first==g1$selected_width)
failed <- bc_collect_candidates(list(error='failure'),'armB',c(.025,.975),FALSE,'test')
missing <- bc_runner_up(failed)
stopifnot(nrow(failed)==15,sum(failed$status=='ineligible')==5,is.na(missing$gap),is.na(missing$within_2))
stat <- bc_near_tie_counts(rbind(g,missing));stopifnot(stat$attempts==2,stat$share_all==.5,stat$share_finite==1)
