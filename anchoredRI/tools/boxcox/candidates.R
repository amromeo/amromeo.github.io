# Tabulate the fitted candidates before the safeguard. All BICs use the original measurement scale.
bc_collect_candidates <- function(fit, method, probs, positive, key) {
  is_boxcox <- method == 'armA'
  families <- if (is_boxcox) c('boxcox_unshifted','boxcox_shifted') else c('identity','log','shifted_log')
  max_k <- if (method == 'armB') 5L else 3L
  tab <- expand.grid(components=seq_len(max_k), transform=families, stringsAsFactors=FALSE)
  tab$eligible <- positive | !tab$transform %in% c('log','boxcox_unshifted')
  tab$lambda <- tab$shift <- tab$bic <- tab$lower <- tab$upper <- NA_real_
  tab$at_bound <- NA
  tab$status <- ifelse(tab$eligible,'unavailable','ineligible')
  one_sided <- length(probs)==1L
  for (i in which(tab$eligible)) {
    k <- tab$components[i]; family <- tab$transform[i]
    if (!is.null(fit$error)) next
    if (is_boxcox) {
      j <- which(vapply(fit$candidates,function(c) c$k==k && identical(c$shifted, family=='boxcox_shifted'),logical(1)))
      stopifnot(length(j)==1L)
      candidate <- fit$candidates[[j]]; tab$bic[i] <- candidate$bic
      if (is.null(candidate$best)) {tab$status[i] <- 'failed'; next}
      p <- candidate$best$parameters; tab$lambda[i] <- p$lambda; tab$shift[i] <- p$shift
      limits <- p$shift + bc_inverse(p$mu[1]+qnorm(probs)*p$sigma[1],p$lambda)
      tab$at_bound[i] <- candidate$best$at_bound
    } else {
      j <- which(fit$selection$transform==family & fit$selection$components==k)
      stopifnot(length(j)==1L)
      tab$bic[i] <- fit$selection$bic_original[j]; tab$at_bound[i] <- fit$selection$at_bound[j]
      tab$lambda[i] <- as.numeric(family=='identity'); tab$shift[i] <- fit$selection$shift[j]
      if(!is.finite(tab$bic[i])) {tab$status[i] <- 'failed'; next}
      details <- fit$family_candidates[[family]]
      stopifnot(!is.null(details),nrow(details$limits)>=k)
      limits <- details$limits[k,]
      if(one_sided) {
        p <- details$parameters[[k]]; stopifnot(!is.null(p))
        y <- p$mu[1]+qnorm(probs)*p$sigma[1]
        limits <- if(family=='identity') y else tab$shift[i]+exp(y)
      }
    }
    tab$status[i] <- if(all(is.finite(limits))) 'ok' else 'unavailable_endpoint'
    if(!one_sided) tab$lower[i] <- limits[1]
    tab$upper[i] <- tail(limits,1)
  }
  valid <- which(is.finite(tab$bic)); ranked <- valid[order(tab$bic[valid],valid)]
  tab$rank <- NA_integer_; tab$rank[ranked] <- seq_along(ranked)
  tab$selected <- !is.na(tab$rank) & tab$rank==1L
  if(is.null(fit$error)) {
    winner <- tab[tab$selected,]; stopifnot(nrow(winner)==1L, winner$components==fit$components,
                                            abs(winner$bic-fit$selected_bic)<=1e-10)
    if(!is_boxcox) stopifnot(winner$transform==fit$transform) else
      stopifnot(isTRUE(all.equal(as.numeric(winner$lambda),as.numeric(fit$lambda),tolerance=0)),
                isTRUE(all.equal(as.numeric(winner$shift),as.numeric(fit$shift),tolerance=0)))
  }
  tab$delta_bic <- if(length(ranked)) tab$bic-tab$bic[ranked[1L]] else NA_real_
  data.frame(key=key, method=method, one_sided=one_sided,tab,row.names=NULL)
}

bc_runner_up <- function(tab) {
  ranks <- match(1:2,tab$rank); top <- tab[ranks,,drop=FALSE]
  gap <- diff(top$bic); lower <- diff(top$lower); upper <- diff(top$upper)
  width <- if(tab$one_sided[1])top$upper[1] else top$upper[1]-top$lower[1]
  endpoints <- if(tab$one_sided[1]) upper else c(lower,upper)
  max_delta <- if(all(is.finite(endpoints)))max(abs(endpoints)) else NA_real_
  data.frame(key=tab$key[1],method=tab$method[1],one_sided=tab$one_sided[1],
    first_transform=top$transform[1],first_k=top$components[1],first_bic=top$bic[1],
    second_transform=top$transform[2],second_k=top$components[2],second_bic=top$bic[2],
    gap=gap,within_2=if(is.finite(gap))gap<=2 else NA,
    lower_first=top$lower[1],upper_first=top$upper[1],lower_second=top$lower[2],upper_second=top$upper[2],
    lower_difference=lower,upper_difference=upper,abs_lower_difference=abs(lower),abs_upper_difference=abs(upper),
    largest_difference=max_delta,selected_width=width,
    fraction_of_width=if(is.finite(width)&&width>0)max_delta/width else NA_real_)
}

bc_near_tie_counts <- function(gaps) {
  valid <- is.finite(gaps$gap); near <- valid & gaps$gap<=2
  delta <- gaps$fraction_of_width[near & is.finite(gaps$fraction_of_width)]
  data.frame(attempts=nrow(gaps),finite_gaps=sum(valid),within_2=sum(near),
    share_all=sum(near)/nrow(gaps),share_finite=if(any(valid))sum(near)/sum(valid) else NA_real_,
    near_with_finite_limits=length(delta),median_fraction_of_width=if(length(delta))median(delta) else NA_real_,
    max_fraction_of_width=if(length(delta))max(delta) else NA_real_)
}
