# Export both original Hb 142 candidate optima on each host for deterministic transport.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
dirs <- setNames(args[1:2], c("local","Actions"))
original <- lapply(dirs, function(path) readRDS(file.path(path, "fit-142.rds")))
stopifnot(original$local$routine == original$Actions$routine, original$local$anchor == original$Actions$anchor)
winners <- unique(unlist(lapply(original, function(fit) {
  r <- fit$base_tab[fit$base_tab$selected,]
  paste(r$transform, r$k)
})))
stopifnot(setequal(winners, c("identity 3","log 3")))
precise <- function(x) sprintf("%.17g", x)
rows <- list()
for (host in names(original)) {
  fit <- original[[host]]
  stopifnot(fit$runtime == host, fit$Index == 142L)
  for (key in winners) {
    m <- fit$models[[key]]
    row <- data.frame(Index=142L, runtime=host, transform=m$transform, k=m$k,
                      n=length(m$par), center=precise(m$cs[1]), scale=precise(m$cs[2]),
                      code=m$code, loglik=precise(m$ll), routine_hash=fit$routine, anchor_hash=fit$anchor)
    for (param in c("par","lo","hi")) {
      vals <- if (param == "par") m$par else m$box[[param]]
      for (j in seq_along(vals)) row[[sprintf("%s%02d",param,j)]] <- precise(vals[j])
    }
    rows[[length(rows)+1L]] <- row
  }
}
dir.create(dirname(args[3]), recursive=TRUE, showWarnings=FALSE)
write.csv(do.call(rbind,rows), args[3], row.names=FALSE)
