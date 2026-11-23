# Independent C2 anchors, with the published order of exclusion and rounding.
bc_anchors <- function(s) {
  set.seed(17000000L+s$Index)
  a <- RIbench:::invBoxCox(rnorm(40,s$nonp_mu,s$nonp_sigma),s$nonp_lambda)+s$nonp_shift
  round(a[is.finite(a)&a>=0],s$decimals)
}

bc_ribench_input <- function(s, directory, expected_anchor_hash) {
  file <- sprintf('%s_%s_seed_%s.csv',s$Index,s$Analyte,s$startSeed)
  x <- data.table::fread(file.path(directory,'Data',s$Analyte,file))$V1+0L
  a <- bc_anchors(s)
  stopifnot(identical(digest::digest(x,algo='md5',serializeVersion=2,ascii=TRUE),s$HashCode),
            identical(digest::digest(a,algo='sha256'),expected_anchor_hash))
  list(value=x,setting=rep('pooled',length(x)),anchor=a,key=paste0('r-',s$Index),s=s,
       stage=if(s$N<=5000)'small' else 'large',probs=if(s$Analyte=='CRP').95 else c(.025,.975))
}
