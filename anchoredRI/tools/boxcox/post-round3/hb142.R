# Evaluate and restart the round 3 constrained likelihood. The production fitter is untouched.
hb142_cross <- function(data, sha_routine, sha_anchor, decimals, runtime, donors) {
  stopifnot(nrow(donors)==4L, all(donors$Index==142), all(donors$routine_hash==sha_routine),
            all(donors$anchor_hash==sha_anchor))
  rules <- r3_rule(data, decimals)
  blocks <- list()
  rows <- vector("list", nrow(donors))
  for (i in seq_len(nrow(donors))) {
    row <- donors[i,]
    m <- r3_from_row(row)
    model <- m$f
    evaluate <- function(par, mode) r3_eval(data,model,m$k,m$cs,par,m$box,rules,mode=mode,maxit=10000L)
    start <- evaluate(m$par, 0L)
    tol <- 128*.Machine$double.eps*max(1,data$n,abs(m$ll))
    stopifnot(abs(start$final-m$ll)<=tol, all(abs(start$par-m$par)<=1e-11*pmax(1,abs(m$par))))
    start_ll <- m$ll
    start_est <- r3_row(m, data)
    best <- m
    par <- m$par
    converged <- FALSE
    for (j in seq_len(100L)) {
      fit <- evaluate(par, 1L)
      blocks[[length(blocks)+1L]] <- data.frame(Index=142L, receiver=runtime, donor=row$runtime,
                        transform=model, k=m$k, block=j, initial_loglik=fit$initial,
                        final_loglik=fit$final, convergence=fit$code, message=fit$message)
      if (is.finite(fit$final) && fit$final > best$ll) {
        best$ll <- fit$final
        best$par <- fit$par
      }
      if (fit$code==0L && is.finite(fit$final) && fit$final>=best$ll-tol) {
        converged <- TRUE
        break
      }
      par <- best$par
    }
    best_est <- r3_row(best,data)
    rows[[i]] <- data.frame(Index=142L, receiver=runtime, donor=row$runtime,
        transform=model,k=m$k, initial_loglik=start_ll, recomputed_start_loglik=start$final,
        transfer_error=start$final-start_ll, initial_lower=start_est$lower, initial_upper=start_est$upper,
        terminal_loglik=fit$final, code=fit$code, converged=converged, blocks=j,
        best_loglik=best$ll, bic=best_est$bic, lower=best_est$lower, upper=best_est$upper,
        gain=best$ll-start_ll)
  }
  list(fits=do.call(rbind,rows),blocks=do.call(rbind,blocks))
}
