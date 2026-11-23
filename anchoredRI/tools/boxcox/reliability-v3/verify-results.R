# Reconcile the published CSV decisions with the original binary64 checkpoint values.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
work <- args[1]
out <- args[2]
read_tab <- function(f) read.csv(file.path(out, f), stringsAsFactors=FALSE)
base <- read_tab('baseline-paired.csv')
continued <- read_tab('arm-a-continued-pairs.csv')
decisions <- read_tab('arm-a-cross-decisions.csv')
retained <- read_tab('arm-a-after-retaining-best.csv')
all_solutions <- read_tab('arm-a-cross-solutions.csv')
arms <- read_tab('arm-b-datasets.csv')
scales <- c('boxcox_unshifted', 'boxcox_shifted')
tol <- .005
phase <- c(local=file.path(work,'local'), Actions=file.path(work,'actions','phase1'))
cross <- c(local=file.path(work,'cross-local'), Actions=file.path(work,'actions','cross'))
fit <- function(path, i, prefix) readRDS(file.path(path, sprintf('%s-%03d.rds', prefix, i)))
family_key <- function(r) paste(r$transform, r$k)
choose <- function(t) t[order(t$bic, match(t$transform, scales), t$k)[1], , drop=FALSE]
check_pair <- function(l, a, row) {
  width <- a$upper-a$lower
  delta <- max(abs(l$lower-a$lower), abs(l$upper-a$upper))
  stopifnot(nrow(row) == 1L, is.finite(width), width>0,
            family_key(l) == paste(row$local_transform, row$local_k),
            family_key(a) == paste(row$actions_transform, row$actions_k),
            (family_key(l)==family_key(a)) == as.logical(row$model_match),
            (delta<=tol*width) == as.logical(row$limit_match),
            abs(delta/width-row$fraction_of_width)<1e-10)
}
for (i in seq_len(156)) {
  saved <- lapply(phase, fit, i=i, prefix='fit')
  stopifnot(identical(saved$local$routine,saved$Actions$routine),
            identical(saved$local$anchor,saved$Actions$anchor))
  bas <- lapply(saved, function(z) z$base_tab[z$base_tab$selected, , drop=FALSE])
  check_pair(bas$local, bas$Actions, base[base$Index==i,])
  winners <- lapply(saved, function(z) z$arm_tab[z$winner, , drop=FALSE])
  stopifnot(all(vapply(winners, function(r) r$code==0L, logical(1))))
  c_row <- continued[continued$Index==i,]
  check_pair(winners$local, winners$Actions, c_row)
  run <- !(as.logical(c_row$model_match) && as.logical(c_row$limit_match))
  stopifnot(run == (i %in% decisions$Index))
  if (run) {
    reciprocal <- lapply(cross, fit, i=i, prefix='cross')
    raw <- do.call(rbind,lapply(reciprocal, function(z) z$tab))
    published <- all_solutions[all_solutions$Index==i,]
    key <- function(z) paste(z$receiver,z$donor,z$transform,z$k)
    stopifnot(setequal(key(raw),key(published)), all(raw$code==0L),
              max(abs(raw$ll-published$ll[match(key(raw),key(published))]))<1e-8)
    leading <- choose(raw)
    nominated <- raw[family_key(raw)==family_key(leading) & raw$donor!=raw$receiver,]
    l <- nominated[nominated$receiver=='local',]
    a <- nominated[nominated$receiver=='Actions',]
    row <- decisions[decisions$Index==i,]
    check_pair(l,a,row)
    dll <- abs(l$ll-a$ll)
    category <- if (dll>1) 'missed better fit' else if(
      max(abs(l$lower-a$lower),abs(l$upper-a$upper))>tol*(a$upper-a$lower)
    ) 'weak identification' else 'concordant'
    stopifnot(category==row$classification, abs(dll-row$abs_loglik_difference)<1e-8)
    fields <- c('transform','k','bic','lower','upper')
    for (h in names(winners)) {
      winners[[h]] <- choose(rbind(saved[[h]]$arm_tab[fields],raw[raw$receiver==h,fields]))
    }
  }
  check_pair(winners$local,winners$Actions,retained[retained$Index==i,])
}
for (i in arms$Index) {
  r <- readRDS(file.path(work,'arm-b',sprintf('fit-%04d.rds',i)))
  tab <- r$nesting
  eligible <- tab[tab$status=='checked',]
  row <- arms[arms$Index==i,]
  stopifnot(nrow(row)==1L,row$checked==nrow(eligible),
            sum(eligible$nested_loglik-eligible$loglik>eligible$roundoff)==row$violations,
            all(r$candidates$minimum_sd_to_floor>=1-1e-11),
            all(abs(r$zero$fixed_loglik-r$zero$zero_shift_loglik)<=r$zero$roundoff))
}
cat('Validated 156 baseline and continued pairs,', nrow(decisions), 'cross-start classifications,',
    nrow(all_solutions), 'terminal solutions and', nrow(arms), 'Arm B checkpoints.\n')
