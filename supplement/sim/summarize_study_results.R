if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("This summary requires the yaml package.", call. = FALSE)
}
protocol <- yaml::read_yaml(file.path("sim", "config", "constrained_exact_protocol.yml"))
output_directory <- file.path("sim", "output", "computational_study")
arguments <- commandArgs(trailingOnly = TRUE)
summary_path <- if (length(arguments)) {
  arguments[[1L]]
} else {
  file.path("..", "recomputed", "constrained_exact_summary.csv")
}

mean_or_na <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
median_or_na <- function(x) if (length(x) && any(is.finite(x))) stats::median(x[is.finite(x)]) else NA_real_
quantile_or_na <- function(x, probability) if (length(x) && any(is.finite(x))) {
  unname(stats::quantile(x[is.finite(x)], probability, names = FALSE))
} else NA_real_
se_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
}
finite_rate <- function(x) mean(is.finite(x))
write_checked <- function(rows, path) {
  if (!nrow(rows) || anyDuplicated(names(rows))) {
    stop(paste("Invalid summary for", path), call. = FALSE)
  }
  write.csv(rows, path, row.names = FALSE, na = "")
}

exact_paths <- unlist(lapply(c(12L, 16L), function(p) {
  file.path(
    output_directory,
    paste0("constrained_exact_p", p, "_",
           names(protocol$exact_constrained_study$data_regimes), ".csv")
  )
}), use.names = FALSE)
missing_exact_paths <- exact_paths[!file.exists(exact_paths)]
if (length(missing_exact_paths)) {
  stop(
    paste("Missing exact-study result files:",
          paste(missing_exact_paths, collapse = ", ")),
    call. = FALSE
  )
}
if (length(exact_paths)) {
  exact <- do.call(rbind, lapply(exact_paths, read.csv, stringsAsFactors = FALSE))
  cells <- split(exact, list(exact$p, exact$regime, exact$epsilon,
                             exact$interval_mode), drop = TRUE)
  exact_summary <- do.call(rbind, lapply(cells, function(cell) data.frame(
    p = cell$p[[1L]], regime = cell$regime[[1L]], epsilon = cell$epsilon[[1L]],
    interval_mode = cell$interval_mode[[1L]], replications = nrow(cell),
    bounds_met_rate = mean(cell$status == "bounds-met"),
    budget_exhausted_rate = mean(cell$status == "budget-exhausted"),
    tolerance_not_met_rate = mean(cell$status == "tolerance-not-met"),
    uninformative_evidence_rate = mean(cell$status == "uninformative-evidence"),
    finite_incumbent_rate = finite_rate(cell$certified_cost_upper),
    simultaneous_coverage_rate = mean(cell$simultaneous_interval_coverage),
    omitted_mass_coverage_rate = mean(cell$omitted_mass_covered),
    epsilon_feasible_rate = mean(cell$epsilon_feasible),
    exact_cost_optimal_rate = mean(cell$exact_cost_optimal),
    mean_true_omitted_mass = mean_or_na(cell$true_omitted_mass),
    se_true_omitted_mass = se_or_na(cell$true_omitted_mass),
    mean_omitted_mass_bound = mean_or_na(cell$omitted_mass_bound),
    se_omitted_mass_bound = se_or_na(cell$omitted_mass_bound),
    mean_exact_optimal_cost = mean_or_na(cell$exact_optimal_cost),
    se_exact_optimal_cost = se_or_na(cell$exact_optimal_cost),
    mean_certified_cost_lower = mean_or_na(cell$certified_cost_lower),
    se_certified_cost_lower = se_or_na(cell$certified_cost_lower),
    mean_certified_cost_upper = mean_or_na(cell$certified_cost_upper),
    se_certified_cost_upper = se_or_na(cell$certified_cost_upper),
    mean_certified_cost_gap = mean_or_na(cell$certified_cost_gap),
    se_certified_cost_gap = se_or_na(cell$certified_cost_gap),
    median_certified_cost_gap = median_or_na(cell$certified_cost_gap),
    q90_certified_cost_gap = quantile_or_na(cell$certified_cost_gap, 0.90),
    mean_selected_cost = mean_or_na(cell$selected_cost),
    mean_true_cost_suboptimality = mean_or_na(cell$true_cost_suboptimality),
    se_true_cost_suboptimality = se_or_na(cell$true_cost_suboptimality),
    mean_maximum_pip_error = mean_or_na(cell$maximum_pip_error),
    se_maximum_pip_error = se_or_na(cell$maximum_pip_error),
    mean_raw_score_requests = mean_or_na(cell$raw_score_requests),
    se_raw_score_requests = se_or_na(cell$raw_score_requests),
    mean_distinct_score_evaluations = mean_or_na(cell$distinct_score_evaluations),
    se_distinct_score_evaluations = se_or_na(cell$distinct_score_evaluations),
    mean_elapsed_seconds = mean_or_na(cell$elapsed_seconds),
    maximum_budget_overshoot = max(cell$budget_overshoot, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
  exact_summary <- exact_summary[order(
    exact_summary$p, exact_summary$regime, exact_summary$epsilon,
    exact_summary$interval_mode
  ), , drop = FALSE]
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  write_checked(exact_summary, summary_path)
}

cat("Wrote the constrained exact-study summary to", summary_path, "\n")
