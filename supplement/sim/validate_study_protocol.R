arguments <- commandArgs(trailingOnly = TRUE)
protocol_path <- if (length(arguments) > 0L) {
  arguments[[1L]]
} else {
  file.path("sim", "config", "constrained_exact_protocol.yml")
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Protocol validation requires the yaml package.", call. = FALSE)
}

protocol <- yaml::read_yaml(protocol_path)
assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
as_numbers <- function(x) as.numeric(unlist(x, use.names = FALSE))
as_characters <- function(x) as.character(unlist(x, use.names = FALSE))

assert(
  identical(protocol$status, "frozen_before_production"),
  "The exact-study protocol is not marked as frozen before production."
)
assert(
  identical(protocol$protocol_id, "regional-compression-exact-study-2026-08-13"),
  "Unexpected protocol identifier."
)

design <- protocol$global_design
epsilon <- as_numbers(design$epsilon)
assert(
  identical(epsilon, c(0.01, 0.05, 0.10, 0.20)) &&
    all(diff(epsilon) > 0) && all(epsilon >= 0 & epsilon < 1),
  "The epsilon grid is invalid."
)
assert(
  identical(as.numeric(design$delta), 0.05) &&
    identical(as.numeric(design$lattice_spacing), 1),
  "The confidence or cost-lattice setting changed."
)
assert(
  as.integer(design$production_replications_per_cell) == 100L,
  "The production replication count changed."
)

primary_cost <- protocol$cost_contract$tree_parameters
primary_values <- c(
  as.numeric(primary_cost$base_cost),
  as_numbers(primary_cost$state_cost)
)
assert(
  all(primary_values >= 1 & primary_values == round(primary_values)),
  "Every declared regional cost must be a positive integer."
)

exact <- protocol$exact_constrained_study
exact_dimensions <- as_numbers(exact$dimensions)
exact_regimes <- names(exact$data_regimes)
interval_modes <- as_characters(exact$interval_modes)
assert(
  identical(exact_dimensions, c(12, 16)) &&
    identical(exact_regimes, c("aligned", "perturbed", "crossed")) &&
    identical(interval_modes, c("exact", "stochastic")),
  "The exact constrained-study cells changed."
)
assert(
  identical(design$evidence_settings$registered_estimator, "iid"),
  "The stochastic evidence estimator changed."
)
production_settings <- design$production_settings
assert(
  as.integer(production_settings$score_budget) == 25000L &&
    as.integer(production_settings$minimum_paths) == 100L &&
    as.integer(production_settings$action_budget) == 100000L,
  "The exact-study production settings changed."
)
exact_cell_count <- length(exact_dimensions) * length(exact_regimes) *
  length(epsilon) * length(interval_modes)
assert(exact_cell_count == 48L, "Expected 48 exact-study summary cells.")
assert(
  as.integer(exact$declared_tree$max_depth$p12) == 5L &&
    as.integer(exact$declared_tree$max_depth$p16) == 6L,
  "The exact-study tree depths changed."
)

# Confirm data reuse within dimension, regime, and replication, while keeping
# algorithm seeds distinct across tolerance and interval mode.
seed_base <- as.integer(protocol$seed_registry$production_base)
seed_rows <- expand.grid(
  p_index = seq_along(exact_dimensions),
  regime_index = seq_along(exact_regimes),
  replication = seq_len(as.integer(design$production_replications_per_cell)),
  epsilon_index = seq_along(epsilon),
  method_index = seq_along(interval_modes),
  KEEP.OUT.ATTRS = FALSE
)
data_seeds <- with(
  seed_rows,
  seed_base + 100000L * p_index + 10000L * regime_index +
    100L * replication
)
algorithm_seeds <- with(
  seed_rows,
  data_seeds + 10L * epsilon_index + method_index
)
assert(
  length(unique(data_seeds)) == length(exact_dimensions) *
    length(exact_regimes) * as.integer(design$production_replications_per_cell),
  "The data seeds do not have the intended reuse pattern."
)
assert(!anyDuplicated(algorithm_seeds), "The algorithm seeds collide.")
assert(
  identical(protocol$acceptance_checks$status, "passed") &&
    length(protocol$acceptance_checks$requirements) >= 6L,
  "The exact-study acceptance checks are incomplete."
)

cat(
  "Exact-study protocol validation passed:",
  protocol$protocol_id,
  sprintf("(%d summary cells).", exact_cell_count),
  "\n"
)
