# SPDX-License-Identifier: GPL-3.0-only
# Decisions use first-stage status, so a released fallback cannot cast a pass vote.
h_usable <- function(lower, upper, probs) {
  if (length(probs) == 1L) is.finite(upper) else
    is.finite(lower) & is.finite(upper) & lower < upper
}

h_middle <- function(x) if (length(x)) stats::median(x) else NA_real_

h_median <- function(draws, probs, passing_only = FALSE) {
  keep <- if (passing_only) draws$first_status %in% "reported" else rep(TRUE, nrow(draws))
  usable <- keep & h_usable(draws$raw_lower, draws$raw_upper, probs)
  shares <- keep & is.finite(draws$share) & draws$share >= 0 & draws$share <= 1
  matched <- usable & shares
  data.frame(lower = if (length(probs) == 1L) NA_real_ else h_middle(draws$raw_lower[usable]),
             upper = h_middle(draws$raw_upper[usable]), contributors = sum(usable),
             share = h_middle(draws$share[shares]), share_n = sum(shares),
             matched_share = h_middle(draws$share[matched]), matched_share_n = sum(matched),
             matched_lower = if (length(probs) == 1L) NA_real_ else h_middle(draws$raw_lower[matched]),
             matched_upper = h_middle(draws$raw_upper[matched]))
}

h_votes <- function(first_status) {
  stopifnot(length(first_status) == 50L)
  # Missing or unrecognized decisions remain unresolved attempts.
  declines <- sum(first_status %in% "declined")
  passes <- sum(first_status %in% "reported")
  c(votes_decline = declines, votes_technical = 50L - declines - passes, votes_pass = passes)
}

h_majority <- function(first_status) {
  votes <- h_votes(first_status)
  if (votes[["votes_decline"]] >= 26L) "declined" else
    if (votes[["votes_decline"]] + votes[["votes_technical"]] <= 25L) "reported" else "withheld"
}

h_gate <- function(decision, lower, upper, probs) {
  if (!h_usable(lower, upper, probs) || !decision %in% c("reported", "declined")) "withheld" else decision
}

h_route <- function(majority, lower, upper, full, probs) {
  if (h_gate(majority, lower, upper, probs) == "reported")
    return(data.frame(route = "median", status = "reported", lower = lower, upper = upper, reason = ""))
  route <- "full_nonreturn"
  if (identical(full$status, "reported")) {
    stopifnot(h_usable(full$lower, full$upper, probs))
    route <- if (identical(full$method, "anchoredRI")) "full_first_stage" else "full_fallback"
  }
  data.frame(route = route, status = full$status, lower = full$lower, upper = full$upper, reason = full$reason)
}

h_decision_row <- function(estimator, rule, decision, raw_lower, raw_upper, probs, reason = "") {
  status <- h_gate(decision, raw_lower, raw_upper, probs)
  if (!h_usable(raw_lower, raw_upper, probs)) reason <- "unavailable estimate"
  data.frame(estimator = estimator, rule = rule, raw_lower = raw_lower, raw_upper = raw_upper,
             lower = if (status == "reported") raw_lower else NA_real_,
             upper = if (status == "reported") raw_upper else NA_real_, status = status,
             route = estimator, reason = reason)
}

h_summarize <- function(fits, probs) {
  standard <- subset(fits, fit_type == "full" & cap == 5L)
  stopifnot(nrow(standard) == 1L)
  medians <- decisions <- list()
  for (cap_value in c(3L, 5L)) {
    within_cap <- fits[fits$cap == cap_value, ]
    full <- subset(within_cap, fit_type == "full")
    w <- subset(within_cap, fit_type == "W")
    stopifnot(nrow(full) == 1L, nrow(w) == 1L)
    for (bank in c("A", "B")) {
      draws <- subset(within_cap, fit_type == "subsample" & seed_set == bank)
      stopifnot(identical(sort(draws$draw), 1:50))
      majority <- h_majority(draws$first_status)
      counts <- as.data.frame(as.list(h_votes(draws$first_status)))
      counts[ paste0("k", 1:5) ] <- as.list(tabulate(draws$components, 5L))
      for (definition in c("all", "passing")) {
        key <- data.frame(Index = standard$Index, cap = cap_value, seed_set = bank, definition = definition)
        m <- h_median(draws, probs, definition == "passing")
        medians[[length(medians) + 1L]] <- cbind(key, m, counts)
        rules <- rbind(
          h_decision_row("median", "full", full$first_status, m$lower, m$upper, probs, full$first_reason),
          h_decision_row("median", "majority", majority, m$lower, m$upper, probs),
          h_decision_row("W", "own", w$first_status, w$raw_lower, w$raw_upper, probs, w$first_reason),
          h_decision_row("W", "full", full$first_status, w$raw_lower, w$raw_upper, probs, full$first_reason),
          h_decision_row("W", "majority", majority, w$raw_lower, w$raw_upper, probs))
        returned <- h_route(majority, m$lower, m$upper, standard, probs)
        # This policy's estimates are the actual returned endpoints (or missing).
        combined <- cbind(data.frame(estimator = "median", rule = "majority_or_package",
                                     raw_lower = returned$lower, raw_upper = returned$upper), returned)
        decisions[[length(decisions) + 1L]] <- cbind(key[rep(1L, nrow(rules) + 1L), ], rbind(rules, combined))
      }
    }
  }
  list(medians = do.call(rbind, medians), decisions = do.call(rbind, decisions))
}
