args <- commandArgs(trailingOnly = TRUE)
stopifnot(getRversion() == "4.5.0", packageVersion("anchoredRI") == "0.4.0")
lock <- read.csv(args[1], colClasses = "character")
observed <- vapply(lock$package, function(pkg) unname(read.dcf(system.file("DESCRIPTION", package = pkg), fields = "Version")[1L]), "")
stopifnot(all(lock$version == observed))
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "Built", "LibPath")]),
          file.path(args[2], "installed-packages.csv"), row.names = FALSE)
capture.output(sessionInfo(), extSoftVersion(), RNGkind(), file = file.path(args[2], "R-session.txt"))
