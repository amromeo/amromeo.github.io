library(anchoredRI)

observe <- function(expr) {
  messages <- warnings <- character()
  value <- withCallingHandlers(expr,
    message = function(m) {
      messages <<- c(messages, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, messages = messages, warnings = warnings)
}

point <- function(x, a, ...) fit_anchored_ri(x, rep("one", length(x)), a, "identity", ...)
# Only the screened values distinguish the true unit from a grid twice as coarse.
x <- rep(c(-4, -2, 0, 2, 4), each = 5)
a <- c(-3, -1, 0, 1, 3)
for (step in c(0.5, 5, 10)) {
  seen <- observe(point(x * step, a * step, fallback = FALSE))
  fit <- seen$value
  stopifnot(identical(fit$rounding_unit, step),
            identical(seen$messages, paste0("Inferred rounding unit is ", step, ".")),
            fit$runner_up$transform == "identity",
            fit$runner_up$components == order(fit$bic_original)[2L],
            fit$runner_up$bic_gap == fit$runner_up$bic_original - fit$selected_bic,
            length(fit$optimizer_status) == 1L, length(fit$optimizer_code) == 1L,
            nzchar(fit$optimizer_message))
}
for (unit in c(0, 0.5)) {
  supplied <- observe(point(x, a, rounding_unit = unit, fallback = FALSE))
  stopifnot(identical(supplied$value$rounding_unit, unit), !length(supplied$messages))
}
infer <- getFromNamespace(".measurement_unit", "anchoredRI")
for (step in c(2e-4, 2.5, 5000, 1e-8))
  stopifnot(infer(c(-3, 0, 1, 2) * step, NULL) == step)
large_decimal <- (1e11 + 1:7) / 20
stopifnot(infer(large_decimal, NULL) == 0.05, infer(large_decimal, 0.05) == 0.05)
continuous <- observe(point(x * sqrt(2), a * sqrt(2), fallback = FALSE))
stopifnot(continuous$value$rounding_unit == 0,
          identical(continuous$messages, "Inferred rounding unit is 0 (no shared increment found)."))

# Only routine results count toward the size threshold. Check every reporting status.
size_notice <- "RIbench accuracy between 5,000 and 50,000 routine results has not been tested directly."
fallback_skip <- "The anchorRefine fallback requires at least 40 screened patients."
anchor_warning <- "Task G evaluated both methods with exactly 40 screened patients."
routine <- rep(10:20, length.out = 5001)
at_cutoff <- observe(point(routine[-1L], a, rounding_unit = 1, safeguard = FALSE))
reported <- observe(point(routine, a, rounding_unit = 1, safeguard = FALSE))
declined <- observe(point(routine, a, rounding_unit = 1))
withheld <- observe(point(routine, rep(0, 5), rounding_unit = 1))
stopifnot(!length(at_cutoff$messages),
          identical(reported$messages, size_notice), reported$value$status == "reported",
          identical(declined$messages, c(size_notice, fallback_skip)), declined$value$status == "declined",
          identical(withheld$messages, c(size_notice, fallback_skip)), withheld$value$status == "withheld",
          all(vapply(list(at_cutoff, reported, declined, withheld),
                     function(x) sum(x$warnings == anchor_warning) == 1L, logical(1))),
          identical(reported$value$method, "anchoredRI"),
          is.null(reported$value$diagnostics$anchorRefine),
          is.na(declined$value$method), identical(declined$value$anchoredRI_status, "declined"),
          identical(declined$value$diagnostics$anchorRefine$status, "skipped"),
          identical(declined$value$diagnostics$anchorRefine$reason, fallback_skip),
          is.na(withheld$value$method), identical(withheld$value$anchoredRI_status, "withheld"),
          all(vapply(list(reported, declined, withheld),
                     function(x) x$value$rounding_unit == 1, logical(1))))
disabled <- observe(point(routine, a, rounding_unit = 1, fallback = FALSE))$value
stopifnot(identical(disabled$limits, declined$value$limits),
          identical(disabled$anchoredRI_reason, declined$value$anchoredRI_reason),
          is.null(disabled$diagnostics$anchorRefine))

boot <- observe(bootstrap_anchored_ri(routine, rep("one", length(routine)), a, "identity",
                                     n_boot = 2, safeguard = FALSE))
b <- boot$value
stopifnot(identical(boot$messages, c(size_notice, "Inferred rounding unit is 1.")),
          b$successful == 2L, b$rounding_unit == 1,
          identical(b$optimizer_status, b$fit$optimizer_status),
          identical(b$runner_up, b$fit$runner_up), all(is.finite(b$bootstrap$bic_gap)),
          all(b$bootstrap$runner_up_transform == "identity"),
          all(!is.na(b$bootstrap$optimizer_code)))

# A hit on the iteration limit leaves an otherwise reportable interval available.
status <- getFromNamespace(".ri_optimizer_status", "anchoredRI")
signal <- getFromNamespace(".signal_ri_fit", "anchoredRI")
apply_safeguard <- getFromNamespace(".apply_ri_safeguard", "anchoredRI")
previous <- options(digits = 1)
message_for_step <- observe(signal(list(rounding_unit = 2.5, status = "declined"), TRUE))
options(previous)
stopifnot(identical(message_for_step$messages, "Inferred rounding unit is 2.5."))
raw <- list(limits = setNames(qnorm(c(.025, .975)), c("lower", "upper")),
            transform = "identity", shift = 0, components = 1L,
            optimizer_status = status(1L, "NEW_X"))
anchors <- qnorm(seq(.02, .98, length.out = 40))
accepted <- apply_safeguard(raw, -3:3, anchors)
rejected <- apply_safeguard(raw, -1:3, anchors)
stopifnot(accepted$status == "reported", identical(accepted$limits, raw$limits),
          identical(observe(signal(accepted, FALSE))$warnings,
                    "The optimizer reached its iteration limit for the reported fit."),
          rejected$status == "declined", !length(observe(signal(rejected, FALSE))$warnings),
          identical(status(0L, "CONVERGENCE"), "converged"),
          identical(status(52L, "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"), "line_search_failure"))
accepted$optimizer_status <- "converged"
stopifnot(!length(observe(signal(accepted, FALSE))$warnings))

# Equal BICs use Task E's family order, then the smaller component count.
rank_models <- getFromNamespace(".rank_ri_candidates", "anchoredRI")
choices <- data.frame(transform = c("log", "identity", "identity", "shifted_log"),
                      components = c(1L, 2L, 1L, 1L), bic_original = c(5, 5, 5, Inf))
ranked <- rank_models(choices)
stopifnot(ranked$winner == 3L, ranked$runner_up$transform == "identity",
          ranked$runner_up$components == 2L, ranked$runner_up$bic_gap == 0)
only <- rank_models(choices[3:4, ])
stopifnot(only$winner == 1L, is.na(only$runner_up$transform), is.na(only$runner_up$bic_gap))
