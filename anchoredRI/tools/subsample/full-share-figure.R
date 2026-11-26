# Plot pre-safeguard errors on the same matched share and endpoint contributors.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- args[1L]
path <- file.path(out, "paired-fits.csv.gz")
names_only <- names(read.csv(path,nrows=0L))
keep <- c("fit_type","seed_set","cap","definition","group_components","matched_gap_pp",
          "matched_d","extreme","status")
classes <- setNames(rep("NULL",length(names_only)),names_only)
classes[c("fit_type","seed_set","definition","status")] <- "character"
classes[c("cap","group_components")] <- "integer"
classes[c("matched_gap_pp","matched_d")] <- "numeric"
# Python exports booleans as True and False; handle the missing outcome explicitly.
classes["extreme"] <- "character"
stopifnot(all(keep %in% names_only))
s <- read.csv(path,colClasses=classes,na.strings=c("","NA"),stringsAsFactors=FALSE)
s <- s[is.finite(s$matched_gap_pp) & is.finite(s$matched_d),]
stopifnot(nrow(s)>0L)
ey <- function(d) log1p(d / .01)
ticks <- c(0,.01,.1,1,10,100,1000,10000)
max_y <- max(s$matched_d)
ticks <- ticks[ticks <= max_y]
x_range <- range(c(0,s$matched_gap_pp)) + c(-1,1)
y_range <- c(0,ey(max_y)*1.04)
for (cap in c(3,5)) for (view in c("all","passing")) {
  tab <- s[s$cap==cap & (s$fit_type=="median" & s$definition==view |
             s$fit_type!="median" & (view=="all" | s$status=="reported")),]
  file <- file.path(out,sprintf("share-gap-cap%d-%s.png",cap,view))
  png(file,width=1900,height=310*cap+120,res=126,type="cairo")
  par(mfrow=c(cap,6),mar=c(1.65,2.35,1.5,.45),oma=c(2.4,2.6,2.8,.6),
      mgp=c(1.25,.3,0),tcl=-.17,cex.axis=.74)
  methods <- list(c("full",""),c("W",""),c("subsample","A"),c("median","A"),
                  c("subsample","B"),c("median","B"))
  labels <- c("Full routine","Arm W","A · draws","A · median","B · draws","B · median")
  for (k in seq_len(cap)) for (j in seq_along(methods)) {
    kind <- methods[[j]][1]; bank <- methods[[j]][2]
    group <- tab[!is.na(tab$group_components) & tab$group_components==k & tab$fit_type==kind &
                 (bank=="" | tab$seed_set %in% bank),]
    plot(NA,xlim=x_range,ylim=y_range,xaxs="i",yaxs="i",xlab="",ylab="",axes=FALSE)
    if (k==cap) axis(1)
    if (j==1) axis(2,at=ey(ticks),labels=format(ticks,trim=TRUE,scientific=FALSE),las=1)
    box(col="#aeb7bb")
    abline(v=0,col="#aeb7bb",lty=3)
    flagged <- group$extreme %in% c("TRUE","True")
    points(group$matched_gap_pp[!flagged],ey(group$matched_d[!flagged]),pch=16,cex=.32,col="#27597268")
    points(group$matched_gap_pp[flagged],ey(group$matched_d[flagged]),pch=17,cex=.48,col="#ae441abb")
    if (k==1L) title(main=labels[j],cex.main=1.05)
    mtext(sprintf("K = %d   n = %s",k,format(nrow(group),big.mark=",")),side=3,line=.02,cex=.58,adj=1)
  }
  mtext(sprintf("Component 1 share and reference-limit error · cap %d · %s",cap,
                if(view=="all") "all eligible fits" else "passing fits"),side=3,outer=TRUE,line=1,
        adj=0,font=2,cex=1.1)
  mtext("Orange triangles mark outcomes with at least one requested endpoint |z| > 5.",
        side=3,outer=TRUE,line=-.1,adj=0,cex=.73)
  mtext("Mean absolute endpoint z error",side=2,outer=TRUE,line=.9,cex=.92)
  mtext("Fitted share minus designed physiological share (percentage points)",side=1,outer=TRUE,line=.8,cex=.91)
  dev.off()
}
