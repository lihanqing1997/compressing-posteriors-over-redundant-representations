source(file.path("sim", "src", "region_restricted_smc.R"))
source(file.path("sim", "src", "exact_gaussian_helpers.R"))
source(file.path("sim", "exact_bma.R"))

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("The exact constrained study requires the yaml package.", call. = FALSE)
}

parse_arguments <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (argument in raw) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Arguments must have the form --name=value.", call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    out[[pieces[1L]]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

split_values <- function(value) {
  if (is.null(value) || !nzchar(value)) NULL else strsplit(value, ",", fixed = TRUE)[[1L]]
}

arguments <- parse_arguments()
allowed_arguments <- c(
  "p", "regime", "interval", "epsilon", "replications",
  "replication-ids", "budget", "n-paths", "output"
)
unknown_arguments <- setdiff(names(arguments), allowed_arguments)
if (length(unknown_arguments) > 0L) {
  stop(
    paste0(
      "Unknown argument",
      if (length(unknown_arguments) == 1L) "" else "s",
      ": ",
      paste0("--", unknown_arguments, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}
protocol <- yaml::read_yaml(file.path(
  "sim", "config", "constrained_exact_protocol.yml"
))
if (!identical(protocol$status, "frozen_before_production")) {
  stop("The exact-study protocol must be frozen before any run.", call. = FALSE)
}

registered_dimensions <- as.integer(unlist(
  protocol$exact_constrained_study$dimensions,
  use.names = FALSE
))
registered_regimes <- names(protocol$exact_constrained_study$data_regimes)
registered_intervals <- as.character(unlist(
  protocol$exact_constrained_study$interval_modes,
  use.names = FALSE
))
registered_epsilon <- as.numeric(unlist(
  protocol$global_design$epsilon,
  use.names = FALSE
))

dimensions <- if (is.null(arguments$p)) {
  registered_dimensions
} else {
  as.integer(split_values(arguments$p))
}
regimes <- if (is.null(arguments$regime)) {
  registered_regimes
} else {
  split_values(arguments$regime)
}
interval_modes <- if (is.null(arguments$interval)) {
  registered_intervals
} else {
  split_values(arguments$interval)
}
epsilon_grid <- if (is.null(arguments$epsilon)) {
  registered_epsilon
} else {
  as.numeric(split_values(arguments$epsilon))
}
if (any(!dimensions %in% registered_dimensions) ||
    any(!regimes %in% registered_regimes) ||
    any(!interval_modes %in% registered_intervals) ||
    any(!epsilon_grid %in% registered_epsilon)) {
  stop("A requested cell is outside the frozen protocol.", call. = FALSE)
}
default_replications <- as.integer(
  protocol$global_design$production_replications_per_cell
)
if (!is.null(arguments$replications) && !is.null(arguments$`replication-ids`)) {
  stop("Use only one of --replications or --replication-ids.", call. = FALSE)
}
replications <- if (!is.null(arguments$`replication-ids`)) {
  as.integer(split_values(arguments$`replication-ids`))
} else if (is.null(arguments$replications)) {
  seq_len(default_replications)
} else {
  requested <- as.integer(split_values(arguments$replications))
  if (length(requested) == 1L) seq_len(requested) else requested
}
if (any(!is.finite(replications)) || any(replications < 1L) ||
    any(replications > default_replications)) {
  stop("Invalid replication request.", call. = FALSE)
}

production_settings <- protocol$global_design$production_settings
score_budget <- if (is.null(arguments$budget)) {
  as.numeric(production_settings$score_budget)
} else {
  as.numeric(arguments$budget)
}
n_paths <- if (is.null(arguments$`n-paths`)) {
  as.integer(production_settings$minimum_paths)
} else {
  as.integer(arguments$`n-paths`)
}
if (score_budget != as.numeric(production_settings$score_budget) ||
    n_paths != as.integer(production_settings$minimum_paths)) {
  stop(
    "Registered production must use the protocol score budget and path count.",
    call. = FALSE
  )
}
max_actions <- as.integer(production_settings$action_budget)
if (!is.finite(score_budget) || score_budget < 1 ||
    !is.finite(max_actions) || max_actions < 1L ||
    !n_paths %in% as.integer(unlist(
      protocol$global_design$evidence_settings$production_paths_grid,
      use.names = FALSE
    ))) {
  stop("Invalid score budget or production path count.", call. = FALSE)
}

default_output <- file.path(
  "..", "recomputed",
  "constrained_exact_production.csv"
)
output_path <- if (is.null(arguments$output)) default_output else arguments$output
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

make_declared_groups <- function(p, regime) {
  aligned <- split(seq_len(p), ceiling(seq_len(p) / 2L))
  if (regime == "aligned") {
    return(aligned)
  }
  out <- aligned
  if (regime == "perturbed") {
    return(out)
  }
  first <- vapply(aligned, `[`, integer(1), 1L)
  second <- vapply(aligned, `[`, integer(1), 2L)
  Map(c, first, c(second[-1L], second[1L]))
}

make_region_oracle <- function(supports,
                               log_scores,
                               groups,
                               region_factory) {
  cache <- new.env(hash = TRUE, parent = emptyenv())
  get_summary <- function(states) {
    states <- as.character(states)
    key <- paste(states, collapse = "/")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    region <- region_factory(states, 1L, key)
    log_proposal <- region$log_prob(supports)
    inside <- is.finite(log_proposal)
    if (!any(inside)) {
      stop("A declared region is empty.", call. = FALSE)
    }
    log_g <- log_scores[inside] - log_proposal[inside]
    summary <- list(
      log_evidence = cb_log_sum_exp(log_scores[inside]),
      log_g_lower = min(log_g),
      log_g_upper = max(log_g),
      rho = max(log_g) - min(log_g),
      inside = inside
    )
    assign(key, summary, envir = cache)
    summary
  }
  list(
    summary = get_summary,
    rho = function(states, depth = NULL) get_summary(states)$rho,
    log_g_bounds = function(states, region = NULL) {
      summary <- get_summary(states)
      c(summary$log_g_lower, summary$log_g_upper)
    },
    log_evidence_bounds = function(states, region = NULL) {
      value <- get_summary(states)$log_evidence
      c(value, value)
    }
  )
}

leaf_reference <- function(supports,
                           posterior_weights,
                           groups,
                           max_depth,
                           base_cost,
                           state_cost) {
  split_groups <- groups[seq_len(max_depth)]
  counts <- vapply(split_groups, function(indices) {
    rowSums(supports[, indices, drop = FALSE])
  }, numeric(nrow(supports)))
  if (is.null(dim(counts))) counts <- matrix(counts, ncol = 1L)
  states <- ifelse(counts == 0, "zero", ifelse(counts == 1, "one", "multi"))
  keys <- apply(states, 1L, paste, collapse = "/")
  masses <- rowsum(posterior_weights, keys, reorder = TRUE)[, 1L]
  state_rows <- strsplit(names(masses), "/", fixed = TRUE)
  costs <- vapply(state_rows, function(row) {
    base_cost + sum(state_cost[row])
  }, numeric(1))
  list(
    keys = names(masses),
    log_mass = log(masses),
    costs = costs
  )
}

data_seed_for <- function(p, regime, replication) {
  base <- as.integer(protocol$seed_registry$production_base)
  p_index <- match(p, registered_dimensions)
  regime_index <- match(regime, registered_regimes)
  base + 100000L * p_index + 10000L * regime_index +
    100L * replication
}

algorithm_seed_for <- function(p, regime, replication, epsilon, interval_mode) {
  data_seed_for(p, regime, replication) +
    10L * match(epsilon, registered_epsilon) +
    match(interval_mode, registered_intervals)
}

dataset_cache <- new.env(hash = TRUE, parent = emptyenv())

make_precomputed_support_score <- function(log_scores, p) {
  cb_make_precomputed_support_score(log_scores, p)
}

prepare_dataset <- function(p, regime, replication) {
  cache_key <- paste(p, regime, replication, sep = "/")
  if (exists(cache_key, envir = dataset_cache, inherits = FALSE)) {
    return(get(cache_key, envir = dataset_cache, inherits = FALSE))
  }
  data_seed <- data_seed_for(p, regime, replication)
  set.seed(data_seed)
  model <- protocol$exact_constrained_study$data_model
  n_train <- as.integer(model$n_train[[paste0("p", p)]])
  data <- cb_block_orthogonal_case(
    p = p,
    n_train = n_train,
    n_test = as.integer(model$n_test),
    rho = as.numeric(model$within_pair_correlation),
    signal = as.numeric(model$signal),
    sigma = as.numeric(model$sigma),
    seed = data_seed,
    heterogeneous = isTRUE(model$heterogeneous_effects)
  )
  if (regime == "perturbed") {
    factor_variance <- as.numeric(model$moderate_cross_group_factor_variance)
    set.seed(data_seed + 1L)
    loadings <- rep(c(-1, 1), length.out = p)
    train_factor <- rnorm(nrow(data$X_train))
    test_factor <- rnorm(nrow(data$X_test))
    data$X_train <- sqrt(1 - factor_variance) * data$X_train +
      sqrt(factor_variance) * tcrossprod(train_factor, loadings)
    data$X_test <- sqrt(1 - factor_variance) * data$X_test +
      sqrt(factor_variance) * tcrossprod(test_factor, loadings)
    data$X_train <- scale(data$X_train)
    data$X_test <- scale(data$X_test)
    data$y_train <- as.numeric(
      data$X_train %*% data$beta + rnorm(nrow(data$X_train), sd = data$sigma)
    )
    data$y_test <- as.numeric(
      data$X_test %*% data$beta + rnorm(nrow(data$X_test), sd = data$sigma)
    )
    response_center <- mean(data$y_train)
    data$y_train <- data$y_train - response_center
    data$y_test <- data$y_test - response_center
  }
  inclusion_probability <- as.numeric(model$inclusion_probability)
  reference_score <- make_known_sigma_support_score(
    data$X_train,
    data$y_train,
    sigma = data$sigma,
    inclusion_probability = inclusion_probability,
    tau2 = as.numeric(model$slab_variance_multiplier)
  )
  supports <- enumerate_supports(p)
  sentinel <- unique(as.integer(c(
    1L, 2L, min(12345L, nrow(supports)),
    floor(nrow(supports) / 2), nrow(supports)
  )))
  expected_index <- 1 + as.numeric(
    supports[sentinel, , drop = FALSE] %*% 2^(seq_len(p) - 1L)
  )
  if (!identical(expected_index, as.numeric(sentinel))) {
    stop("The enumerated-support order does not match the lookup oracle.", call. = FALSE)
  }
  log_scores <- rrs_evaluate_log_score(reference_score$score, supports)
  log_normalizer <- cb_log_sum_exp(log_scores)
  posterior_weights <- exp(log_scores - log_normalizer)
  prepared <- list(
    data = data,
    supports = supports,
    log_scores = log_scores,
    posterior_weights = posterior_weights,
    full_pip = colSums(supports * posterior_weights)
  )
  assign(cache_key, prepared, envir = dataset_cache)
  prepared
}

run_cell <- function(p, regime, replication, epsilon, interval_mode) {
  data_seed <- data_seed_for(p, regime, replication)
  algorithm_seed <- algorithm_seed_for(
    p, regime, replication, epsilon, interval_mode
  )
  model <- protocol$exact_constrained_study$data_model
  prepared <- prepare_dataset(p, regime, replication)
  data <- prepared$data
  inclusion_probability <- as.numeric(model$inclusion_probability)
  score <- make_precomputed_support_score(prepared$log_scores, p)
  supports <- prepared$supports
  log_scores <- prepared$log_scores
  posterior_weights <- prepared$posterior_weights
  full_pip <- prepared$full_pip
  groups <- make_declared_groups(p, regime)
  score$reset()
  screened_proposal <- make_screened_group_region_factory(
    log_score_fn = score$score,
    p = p,
    groups = groups
  )
  screening_requests <- unname(score$diagnostics()[["calls"]])
  if (screening_requests >= score_budget) {
    stop("Screening consumed the entire registered score budget.", call. = FALSE)
  }
  oracle <- make_region_oracle(
    supports,
    log_scores,
    groups,
    screened_proposal$region_factory
  )
  primary_cost <- protocol$cost_contract$tree_parameters
  state_cost <- as.numeric(unlist(primary_cost$state_cost))
  names(state_cost) <- names(primary_cost$state_cost)
  reference <- leaf_reference(
    supports,
    posterior_weights,
    groups,
    max_depth = as.integer(
      protocol$exact_constrained_study$declared_tree$max_depth[[paste0("p", p)]]
    ),
    base_cost = as.numeric(primary_cost$base_cost),
    state_cost = state_cost
  )
  exact_optimum <- select_region_union_epsilon(
    log_evidence_lower = reference$log_mass,
    log_evidence_upper = reference$log_mass,
    costs = reference$costs,
    epsilon = epsilon,
    mode = "conservative"
  )
  tree <- new_group_state_tree(
    inclusion_probabilities = rep(inclusion_probability, p),
    groups = groups,
    split_order = seq_along(groups),
    max_depth = as.integer(
      protocol$exact_constrained_study$declared_tree$max_depth[[paste0("p", p)]]
    ),
    rho = oracle$rho,
    log_g_bounds = oracle$log_g_bounds,
    log_evidence_bounds = if (interval_mode == "exact") {
      oracle$log_evidence_bounds
    } else {
      NULL
    },
    target_mutation = make_group_gibbs_target_mutation,
    region_factory = screened_proposal$region_factory,
    state_cost = state_cost,
    base_cost = as.numeric(primary_cost$base_cost)
  )
  set.seed(algorithm_seed)
  elapsed <- system.time({
    search <- if (interval_mode == "exact") {
      select_minimum_cost_region_tree(
        tree = tree,
        log_score_fn = score$score,
        evidence_evaluator = function(node, remaining_score_budget = Inf) {
          stop("An exact-bounded terminal node was unexpectedly reevaluated.")
        },
        epsilon = epsilon,
        cost_tolerance = as.numeric(protocol$global_design$cost_gap_tolerance),
        max_score_evaluations = score_budget - screening_requests,
        max_actions = max_actions
      )
    } else {
      fit_certified_region_selection(
        tree = tree,
        log_score_fn = score$score,
        delta = as.numeric(protocol$global_design$delta),
        n_pilot = as.integer(protocol$global_design$evidence_settings$pilot_paths),
        n_paths = n_paths,
        pilot_mutation_steps = as.integer(
          protocol$global_design$evidence_settings$pilot_mutation_steps
        ),
        production_mutation_steps = as.integer(
          protocol$global_design$evidence_settings$production_mutation_steps
        ),
        maximum_relative_radius = as.numeric(
          protocol$global_design$evidence_settings$maximum_relative_radius
        ),
        maximum_paths = as.integer(
          protocol$global_design$evidence_settings$maximum_production_paths_per_node
        ),
        evidence_method = as.character(
          protocol$global_design$evidence_settings$registered_estimator
        ),
        epsilon = epsilon,
        cost_tolerance = as.numeric(protocol$global_design$cost_gap_tolerance),
        max_score_evaluations = score_budget - screening_requests,
        max_actions = max_actions
      )
    }
  })[["elapsed"]]

  selected_nodes <- search$nodes[search$selected_ids]
  selected_supports <- rep(FALSE, nrow(supports))
  if (length(selected_nodes) > 0L) {
    for (node in selected_nodes) {
      selected_supports <- selected_supports |
        is.finite(node$region$log_prob(supports))
    }
  }
  retained_mass <- sum(posterior_weights[selected_supports])
  true_omitted_mass <- 1 - retained_mass
  selected_cost <- if (length(selected_nodes) == 0L) {
    Inf
  } else {
    sum(vapply(selected_nodes, function(node) node$region$cost, numeric(1)))
  }
  restricted_pip_error <- if (retained_mass > 0) {
    restricted_weights <- posterior_weights * selected_supports / retained_mass
    max(abs(colSums(supports * restricted_weights) - full_pip))
  } else {
    NA_real_
  }
  bounded_nodes <- search$nodes[vapply(search$nodes, function(node) {
    is.finite(node$log_lower) && is.finite(node$log_upper)
  }, logical(1))]
  node_coverage <- vapply(bounded_nodes, function(node) {
    truth <- oracle$summary(node$region$states)$log_evidence
    node$log_lower <= truth && truth <= node$log_upper
  }, logical(1))
  diagnostics <- score$diagnostics()
  preflight_nodes <- search$nodes[vapply(search$nodes, function(node) {
    !is.null(node$evidence_preflight)
  }, logical(1))]
  required_paths <- if (length(preflight_nodes) > 0L) {
    vapply(preflight_nodes, function(node) {
      node$evidence_preflight$required_paths
    }, numeric(1))
  } else {
    numeric(0)
  }
  scheduled_paths <- if (length(preflight_nodes) > 0L) {
    vapply(preflight_nodes, function(node) {
      if (is.null(node$evidence_preflight$scheduled_paths)) {
        NA_real_
      } else {
        node$evidence_preflight$scheduled_paths
      }
    }, numeric(1))
  } else {
    numeric(0)
  }
  scheduled_paths <- scheduled_paths[is.finite(scheduled_paths)]
  accounting <- rrs_common_accounting(
    diagnostics,
    search = search,
    downstream_score_requests = 0,
    elapsed_seconds = elapsed
  )
  data.frame(
    protocol_id = protocol$protocol_id,
    run_mode = "production",
    registered_production = TRUE,
    p = p,
    regime = regime,
    replication = replication,
    epsilon = epsilon,
    delta = as.numeric(protocol$global_design$delta),
    interval_mode = interval_mode,
    score_oracle = "precomputed_exact_table",
    evidence_method = if (interval_mode == "stochastic") {
      as.character(protocol$global_design$evidence_settings$registered_estimator)
    } else {
      "exact"
    },
    data_seed = data_seed,
    algorithm_seed = algorithm_seed,
    n_paths = if (interval_mode == "stochastic") n_paths else NA_integer_,
    score_budget = score_budget,
    raw_status = search$status,
    status = search$registered_status,
    simultaneous_interval_coverage = all(node_coverage),
    covered_nodes = sum(node_coverage),
    bounded_nodes = length(node_coverage),
    exact_optimal_cost = exact_optimum$budget,
    certified_cost_lower = search$cost_lower,
    certified_cost_upper = search$cost_upper,
    certified_cost_gap = search$cost_gap,
    selected_cost = selected_cost,
    true_cost_suboptimality = selected_cost - exact_optimum$budget,
    exact_cost_optimal = is.finite(selected_cost) &&
      selected_cost == exact_optimum$budget,
    retained_mass = retained_mass,
    true_omitted_mass = true_omitted_mass,
    omitted_mass_bound = search$omitted_mass_bound,
    omitted_mass_covered = true_omitted_mass <= search$omitted_mass_bound + 1e-12,
    epsilon_feasible = true_omitted_mass <= epsilon + 1e-12,
    cost_sandwich_covered = search$cost_lower <= exact_optimum$budget &&
      exact_optimum$budget <= search$cost_upper,
    total_variation = true_omitted_mass,
    reverse_kullback_leibler = if (retained_mass > 0) -log(retained_mass) else Inf,
    maximum_pip_error = restricted_pip_error,
    search_actions = search$actions,
    search_expansions = search$expansions,
    screening_score_requests = screening_requests,
    preflight_nodes = length(preflight_nodes),
    preflight_blocked_nodes = sum(vapply(search$nodes, function(node) {
      isTRUE(node$preflight_blocked)
    }, logical(1))),
    minimum_required_paths = if (length(required_paths) > 0L) {
      min(required_paths)
    } else {
      NA_real_
    },
    maximum_required_paths = if (length(required_paths) > 0L) {
      max(required_paths)
    } else {
      NA_real_
    },
    minimum_scheduled_paths = if (length(scheduled_paths) > 0L) {
      min(scheduled_paths)
    } else {
      NA_real_
    },
    maximum_scheduled_paths = if (length(scheduled_paths) > 0L) {
      max(scheduled_paths)
    } else {
      NA_real_
    },
    budget_overshoot = search$budget_overshoot,
    accounting,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

rows <- list()
index <- 0L
bind_output_rows <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- NA
    row[all_names]
  })
  do.call(rbind, rows)
}
for (p in dimensions) {
  for (regime in regimes) {
    for (replication in replications) {
      for (epsilon in epsilon_grid) {
        for (interval_mode in interval_modes) {
          index <- index + 1L
          rows[[index]] <- tryCatch(
            run_cell(p, regime, replication, epsilon, interval_mode),
            error = function(error) data.frame(
              protocol_id = protocol$protocol_id,
              run_mode = "production",
              registered_production = TRUE,
              p = p,
              regime = regime,
              replication = replication,
              epsilon = epsilon,
              delta = as.numeric(protocol$global_design$delta),
              interval_mode = interval_mode,
              score_oracle = "precomputed_exact_table",
              evidence_method = if (interval_mode == "stochastic") {
                as.character(
                  protocol$global_design$evidence_settings$registered_estimator
                )
              } else {
                "exact"
              },
              data_seed = data_seed_for(p, regime, replication),
              algorithm_seed = algorithm_seed_for(
                p, regime, replication, epsilon, interval_mode
              ),
              n_paths = if (interval_mode == "stochastic") n_paths else NA_integer_,
              score_budget = score_budget,
              raw_status = "implementation_error",
              status = "implementation-error",
              error_message = conditionMessage(error),
              stringsAsFactors = FALSE
              )
          )
          message(sprintf(
            "[%d] p=%d regime=%s rep=%d epsilon=%.2f interval=%s status=%s",
            index, p, regime, replication, epsilon, interval_mode,
            rows[[index]]$status[1L]
          ))
          write.csv(
            bind_output_rows(rows),
            paste0(output_path, ".partial"),
            row.names = FALSE,
            na = ""
          )
        }
      }
    }
  }
}

output <- bind_output_rows(rows)
write.csv(output, output_path, row.names = FALSE, na = "")
partial_path <- paste0(output_path, ".partial")
if (file.exists(partial_path)) unlink(partial_path)
cat(sprintf(
  "Wrote %d registered production rows to %s.\n",
  nrow(output), normalizePath(output_path, winslash = "/", mustWork = FALSE)
))
