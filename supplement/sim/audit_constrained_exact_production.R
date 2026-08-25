protocol_path <- file.path("sim", "config", "constrained_exact_protocol.yml")
arguments <- commandArgs(trailingOnly = TRUE)
audit_path <- if (length(arguments)) {
  arguments[[1L]]
} else {
  file.path(
    "..", "recomputed",
    "constrained_exact_audit.csv"
  )
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("The production audit requires the yaml package.", call. = FALSE)
}
protocol <- yaml::read_yaml(protocol_path)

output_directory <- file.path("sim", "output", "computational_study")
registered_files <- unlist(lapply(c(12L, 16L), function(p) {
  file.path(
    output_directory,
    paste0(
      "constrained_exact_p", p, "_",
      names(protocol$exact_constrained_study$data_regimes),
      ".csv"
    )
  )
}), use.names = FALSE)
missing_files <- registered_files[!file.exists(registered_files)]
if (length(missing_files) > 0L) {
  stop(
    paste("Missing registered production files:", paste(missing_files, collapse = ", ")),
    call. = FALSE
  )
}

rows <- do.call(rbind, lapply(registered_files, read.csv, stringsAsFactors = FALSE))
if (!"score_oracle" %in% names(rows) ||
    any(rows$score_oracle != "precomputed_exact_table")) {
  stop(
    "Every exact-study row must identify the frozen precomputed exact-table score oracle.",
    call. = FALSE
  )
}
expected_rows <- length(protocol$exact_constrained_study$dimensions) *
  length(protocol$exact_constrained_study$data_regimes) *
  as.integer(protocol$global_design$production_replications_per_cell) *
  length(protocol$global_design$epsilon) *
  length(protocol$exact_constrained_study$interval_modes)
if (nrow(rows) != expected_rows) {
  stop(sprintf("Expected %d production rows but found %d.", expected_rows, nrow(rows)), call. = FALSE)
}

key <- with(rows, paste(p, regime, replication, epsilon, interval_mode, sep = "/"))
if (anyDuplicated(key)) stop("Registered production keys are duplicated.", call. = FALSE)
if (any(!rows$registered_production)) stop("A production row is not registered.", call. = FALSE)
if (any(rows$score_budget != protocol$global_design$production_settings$score_budget)) {
  stop("A production row used the wrong score budget.", call. = FALSE)
}
if (any(rows$budget_overshoot > 0, na.rm = TRUE)) {
  stop("A production run overshot its score budget.", call. = FALSE)
}
if (anyDuplicated(rows$algorithm_seed)) {
  stop("Registered algorithm seeds collide.", call. = FALSE)
}

dataset_key <- with(rows, paste(p, regime, replication, sep = "/"))
data_seed_counts <- tapply(rows$data_seed, dataset_key, function(x) length(unique(x)))
if (any(data_seed_counts != 1L)) stop("A dataset was not reused across methods.", call. = FALSE)

allowed_status <- as.character(unlist(protocol$global_design$status_levels, use.names = FALSE))
unknown_status <- setdiff(unique(rows$status), allowed_status)
if (length(unknown_status) > 0L) {
  stop(paste("Unknown status:", paste(unknown_status, collapse = ", ")), call. = FALSE)
}
if (any(rows$status == "implementation-error")) {
  stop("At least one registered row has implementation-error status.", call. = FALSE)
}

exact_rows <- rows$interval_mode == "exact"
stochastic_rows <- rows$interval_mode == "stochastic"
deterministic_pass <- all(
  rows$cost_sandwich_covered[exact_rows] &
    rows$exact_cost_optimal[exact_rows] &
    rows$status[exact_rows] == "bounds-met"
)
coverage_successes <- sum(rows$simultaneous_interval_coverage[stochastic_rows])
coverage_trials <- sum(stochastic_rows)
coverage_rate <- coverage_successes / coverage_trials

finite_incumbent <- is.finite(rows$certified_cost_upper)
aligned_epsilon <- stochastic_rows & rows$regime == "aligned" & rows$epsilon == 0.10
perturbed_epsilon <- stochastic_rows & rows$regime == "perturbed" & rows$epsilon == 0.10
aligned_rate <- mean(finite_incumbent[aligned_epsilon])
perturbed_rate <- mean(finite_incumbent[perturbed_epsilon])

summary <- data.frame(
  check = c(
    "expected_rows",
    "deterministic_exactness_and_sandwich",
    "stochastic_simultaneous_coverage_rate",
    "aligned_finite_incumbent_rate_epsilon_0.10",
    "perturbed_finite_incumbent_rate_epsilon_0.10",
    "budget_overshoots",
    "implementation_errors"
  ),
  value = c(
    nrow(rows),
    as.numeric(deterministic_pass),
    coverage_rate,
    aligned_rate,
    perturbed_rate,
    sum(rows$budget_overshoot > 0, na.rm = TRUE),
    sum(rows$status == "implementation-error")
  ),
  required = c(
    expected_rows,
    1,
    0.94,
    0.80,
    0.60,
    0,
    0
  ),
  comparison = c("equal", "equal", "at_least", "at_least", "at_least", "equal", "equal"),
  stringsAsFactors = FALSE
)
summary$passed <- with(summary,
  (comparison == "equal" & value == required) |
    (comparison == "at_least" & value >= required) |
    (comparison == "greater_than" & value > required)
)

dir.create(dirname(audit_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, audit_path, row.names = FALSE)
if (!all(summary$passed)) {
  failed <- summary$check[!summary$passed]
  stop(paste("Protocol audit failed:", paste(failed, collapse = ", ")), call. = FALSE)
}
cat(sprintf(
  "Protocol audit passed for %d registered rows; wrote %s.\n",
  nrow(rows), audit_path
))
