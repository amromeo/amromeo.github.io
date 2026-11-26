# Exercise the production timer without fitting, including an early return and
# a fallback that overwrites its own cleanup handler.
source("anchoredRI/tools/subsample/pilot.R")
fixture <- new.env(parent = baseenv())
evalq({
  cleaned <- FALSE
  .fit_anchor_refine <- function() {
    on.exit(assign("cleaned", TRUE, envir = parent.env(environment())))
    return("reported")
  }
  public <- function(use) {
    fit <- list(diagnostics = list(anchorRefine = NULL))
    if (use) fit$diagnostics$anchorRefine <- .fit_anchor_refine()
    fit
  }
}, fixture)
expected <- fixture$public(TRUE)
fixture$cleaned <- FALSE
called <- h_observe_fallback(fixture$public(TRUE), fixture)
skipped <- h_observe_fallback(fixture$public(FALSE), fixture)
checks <- c(public_result = identical(called$fits[["5"]], expected), cleanup = fixture$cleaned,
            fallback_clock = length(called$fallback_time) == 2L && all(is.finite(called$fallback_time)) &&
              all(called$fallback_time >= 0),
            skipped_clock = length(skipped$fallback_time) == 0L,
            skipped_diagnostics = is.null(skipped$fits[["5"]]$diagnostics$anchorRefine))
stopifnot(all(checks))
write.csv(data.frame(check = names(checks), passed = unname(checks)), commandArgs(TRUE)[1], row.names = FALSE)
