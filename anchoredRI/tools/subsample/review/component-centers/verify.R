# SPDX-License-Identifier: GPL-3.0-only
# Independent check using RIbench's inverse transform and the healthy quartiles.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, packageVersion("RIbench") == "1.0.2")
inputs <- args[1]
out <- args[2]
design <- read.csv(file.path(inputs, "testsets.csv"))
params <- subset(read.csv(file.path(inputs, "selected-component-parameters.csv")), method == "armB")
scores <- subset(read.csv(file.path(inputs, "ribench-outcomes.csv")), method == "armB" & components == 5)
archived <- subset(read.csv(file.path(inputs, "armB-fits.csv")), components == 5)
x <- merge(archived[c("Index", "dataset_id")], params, by = "dataset_id")
screened <- subset(x, screened_component)
x <- merge(subset(x, !screened_component), design, by = "Index")

# The Task E export is in measurement units for identity. Its log-family
# locations use log(x-shift). RIbench adds its own shift after inversion.
cx <- ifelse(x$transform == "identity", x$transformed_mean, x$shift + exp(x$transformed_mean))
quartiles <- function(i) {
    lam <- x$nonp_lambda[i]
    shift <- x$nonp_shift[i]
    min.y <- if (lam == 1) -1-shift else if (lam == 0) -Inf else -1/lam
    mu <- x$nonp_mu[i]
    sigma <- x$nonp_sigma[i]
    left <- pnorm(min.y, mu, sigma)
    y <- qnorm(left+(1-left)*c(.25,.75), mu, sigma)
    RIbench:::invBoxCox(y, lam) + shift
}
q <- vapply(seq_len(nrow(x)), quartiles, numeric(2))
b <- data.frame(Index=x$Index, component=x$component, center=cx, weight=x$routine_weight,
                reference=(x$Analyte == "CRP" | cx >= x$GT_LRL) & cx <= x$GT_URL,
                central50=cx >= q[1,] & cx <= q[2,])
parts <- merge(b, read.csv(file.path(out, "components.csv")), by = c("Index", "component"))
stopifnot(nrow(parts) == 3764, max(abs(parts$center - parts$center_original)) < 1e-10,
          all(parts$reference == as.logical(parts$in_reference)),
          all(parts$central50 == as.logical(parts$in_central50)))
per.fit <- read.csv(file.path(out, "datasets.csv"))
grouped <- aggregate(cbind(reference_n=as.integer(b$reference), central50_n=as.integer(b$central50),
                           reference_weight_pct=100*b$weight*b$reference,
                           central50_weight_pct=100*b$weight*b$central50) ~ Index,
                     data = b, FUN = sum)
grouped <- grouped[match(per.fit$Index, grouped$Index),]
own <- scores[match(per.fit$Index, scores$Index),]
ref <- screened[match(per.fit$Index, screened$Index),]
truth <- design[match(per.fit$Index, design$Index),]
shortfall <- 100*(1-truth$fractionPathol-truth$bg_fraction-ref$routine_weight)
stopifnot(nrow(per.fit) == 941, all(own$Index == per.fit$Index), all(ref$Index == per.fit$Index),
          max(abs(per.fit$d-own$zzDevAbs_Ov)) < 1e-12, max(abs(shortfall - per.fit$shortfall_pp)) < 1e-10)
for (fld in setdiff(names(grouped), "Index"))
    stopifnot(max(abs(per.fit[[fld]] - grouped[[fld]])) < 1e-10)
s <- subset(read.csv(file.path(out,"summary.csv")), axis == "all")
for (i in seq_len(nrow(s))) {
    weights <- grouped[[paste0(s$region[i], "_weight_pct")]]
    n <- grouped[[paste0(s$region[i], "_n")]]
    stopifnot(sum(n>0) == s$fits_with_center[i], sum(n) == s$centers_inside[i],
              abs(cor(weights, per.fit$d, method="spearman") - s$rho_weight_d[i]) < 1e-10,
              abs(cor(weights, per.fit$shortfall_pp, method="spearman") - s$rho_weight_shortfall[i]) < 1e-10)
}
cat("Matched both region indicators for all 3,764 unscreened centers and both weight sums on 941 fits.\n")
cat("All paired shortfalls, raw errors, cohort counts and four pooled Spearman correlations match.\n")
