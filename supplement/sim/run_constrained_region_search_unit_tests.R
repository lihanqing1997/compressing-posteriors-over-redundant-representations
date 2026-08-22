source(file.path("sim", "src", "region_restricted_smc.R"))

fail <- function(message) stop(message, call. = FALSE)
expect_true <- function(value, message) if (!isTRUE(value)) fail(message)
expect_error <- function(expression, message) {
  failed <- inherits(try(force(expression), silent = TRUE), "try-error")
  if (!failed) fail(message)
}

brute_force_minimum_cost <- function(log_lower,
                                     log_upper,
                                     costs,
                                     epsilon,
                                     mode,
                                     eligible = rep(TRUE, length(costs))) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) {
    return(list(feasible = FALSE, budget = Inf, selected = integer(0)))
  }
  best_budget <- Inf
  best_selected <- integer(0)
  for (mask in seq_len(2^length(eligible_ids) - 1L)) {
    chosen_bits <- as.logical(intToBits(mask)[seq_along(eligible_ids)])
    selected <- eligible_ids[chosen_bits]
    summary <- rrs_interval_omitted_mass(log_lower, log_upper, selected)
    if (rrs_interval_epsilon_feasible(summary, epsilon, mode)) {
      budget <- sum(costs[selected])
      if (budget < best_budget) {
        best_budget <- budget
        best_selected <- selected
      }
    }
  }
  list(
    feasible = is.finite(best_budget),
    budget = best_budget,
    selected = best_selected
  )
}

# The older necessary condition U(S) >= (1-epsilon)L(F) accepts S={1}
# here, but the sharp ratio lower bound is 1/(2+1), which exceeds epsilon.
# The optimistic solver must therefore reject the one-cost subset.
sharp_lower <- log(c(1, 1))
sharp_upper <- log(c(2, 1))
sharp_epsilon <- 0.1
weak_condition_holds <- exp(sharp_upper[1L]) >=
  (1 - sharp_epsilon) * sum(exp(sharp_lower))
sharp_selection <- select_region_union_epsilon(
  log_evidence_lower = sharp_lower,
  log_evidence_upper = sharp_upper,
  costs = c(1L, 1L),
  epsilon = sharp_epsilon,
  mode = "optimistic"
)
expect_true(
  weak_condition_holds &&
    sharp_selection$feasible &&
    sharp_selection$budget == 2L &&
    identical(sharp_selection$selected, 1:2),
  "the optimistic solver used the weak total-frontier condition"
)

# Exhaust every subset for randomized small interval frontiers.  This checks
# both the certified incumbent and the sharp optimistic lower-bound problem.
set.seed(20260813)
for (trial in seq_len(500L)) {
  n_regions <- sample.int(8L, 1L)
  truth <- exp(runif(n_regions, -8, 8))
  lower_factor <- exp(runif(n_regions, 0, 2))
  upper_factor <- exp(runif(n_regions, 0, 2))
  log_lower <- log(truth / lower_factor)
  log_upper <- log(truth * upper_factor)
  if (runif(1L) < 0.1) {
    log_lower[sample.int(n_regions, 1L)] <- -Inf
  }
  costs <- sample.int(5L, n_regions, replace = TRUE)
  epsilon <- sample(c(0, 1e-8, 0.01, 0.05, 0.2, 0.5, 0.9), 1L)
  eligible <- runif(n_regions) > 0.25
  for (mode in c("conservative", "optimistic")) {
    exact <- brute_force_minimum_cost(
      log_lower,
      log_upper,
      costs,
      epsilon,
      mode,
      eligible
    )
    dynamic <- select_region_union_epsilon(
      log_evidence_lower = log_lower,
      log_evidence_upper = log_upper,
      costs = costs,
      epsilon = epsilon,
      mode = mode,
      eligible = eligible
    )
    expect_true(
      identical(dynamic$feasible, exact$feasible) &&
        (!exact$feasible || dynamic$budget == exact$budget),
      sprintf(
        paste0(
          "subset optimizer disagreed with exhaustive search at trial %d (%s): ",
          "expected feasible=%s budget=%s; got feasible=%s budget=%s; ",
          "epsilon=%s lower=%s upper=%s costs=%s eligible=%s"
        ),
        trial, mode, exact$feasible, exact$budget,
        dynamic$feasible, dynamic$budget, epsilon,
        paste(log_lower, collapse = ","),
        paste(log_upper, collapse = ","),
        paste(costs, collapse = ","),
        paste(eligible, collapse = ",")
      )
    )
    if (dynamic$feasible) {
      interval <- rrs_interval_omitted_mass(
        log_lower,
        log_upper,
        dynamic$selected
      )
      reported <- if (identical(mode, "conservative")) {
        interval$upper
      } else {
        interval$lower
      }
      expect_true(
        reported <= epsilon + 1e-12,
        "a returned interval subset violated its defining ratio certificate"
      )
    }
  }
}

# Whenever the true evidence lies inside every interval, the optimistic and
# conservative costs must bracket the exact constrained optimum.  A
# conservative incumbent must also satisfy the true omitted-mass tolerance.
set.seed(20260814)
for (trial in seq_len(250L)) {
  n_regions <- sample.int(8L, 1L)
  truth <- exp(runif(n_regions, -5, 5))
  log_truth <- log(truth)
  log_lower <- log_truth - runif(n_regions, 0, 2)
  log_upper <- log_truth + runif(n_regions, 0, 2)
  costs <- sample.int(5L, n_regions, replace = TRUE)
  epsilon <- runif(1L, 0, 0.9)
  optimistic <- select_region_union_epsilon(
    log_lower,
    costs,
    epsilon,
    log_upper,
    "optimistic"
  )
  conservative <- select_region_union_epsilon(
    log_lower,
    costs,
    epsilon,
    log_upper,
    "conservative"
  )
  exact <- select_region_union_epsilon(
    log_truth,
    costs,
    epsilon,
    mode = "conservative"
  )
  expect_true(
    optimistic$budget <= exact$budget &&
      exact$budget <= conservative$budget,
    sprintf("cost bracket missed the exact optimum at trial %d", trial)
  )
  true_omitted <- sum(truth[-conservative$selected]) / sum(truth)
  expect_true(
    true_omitted <= epsilon + 1e-12,
    "a conservative incumbent violated the true mass tolerance"
  )
}

# End-to-end exact star trees check that the selection-only API returns the
# same minimum cost as exhaustive subset optimization and never invokes local
# approximation or a redundant terminal evaluator.
set.seed(20260815)
for (trial in seq_len(100L)) {
  n_leaves <- sample.int(7L, 1L) + 1L
  supports <- diag(n_leaves)
  evidence <- exp(runif(n_leaves, -4, 4))
  costs <- sample.int(5L, n_leaves, replace = TRUE)
  epsilon <- runif(1L, 0, 0.8)
  keys <- rrs_support_keys(supports)
  log_score_fn <- function(gamma) {
    log(evidence[match(paste(gamma, collapse = ","), keys)])
  }
  root_region <- new_enumerated_region(
    supports,
    cost = min(costs),
    label = "root"
  )
  leaf_nodes <- lapply(seq_len(n_leaves), function(j) {
    new_region_tree_node(
      id = paste0("leaf-", j),
      region = new_enumerated_region(
        supports[j, , drop = FALSE],
        cost = costs[j],
        label = paste0("leaf-", j)
      ),
      rho = 0,
      cost_lower = costs[j],
      terminal = TRUE
    )
  })
  tree <- new_region_tree(
    c(
      list(new_region_tree_node(
        id = "root",
        region = root_region,
        children = vapply(leaf_nodes, `[[`, character(1), "id"),
        rho = diff(range(log(evidence))),
        cost_lower = min(costs),
        terminal = FALSE
      )),
      leaf_nodes
    ),
    "root"
  )
  exact <- brute_force_minimum_cost(
    log(evidence),
    log(evidence),
    costs,
    epsilon,
    "conservative"
  )
  selection <- select_minimum_cost_region_tree(
    tree = tree,
    log_score_fn = log_score_fn,
    evidence_evaluator = function(node) {
      stop("point-bounded terminal nodes must not be reevaluated")
    },
    epsilon = epsilon,
    cost_tolerance = 0,
    max_actions = 2L
  )
  selected_cost <- sum(vapply(
    selection$nodes[selection$selected_ids],
    function(node) node$region$cost,
    numeric(1)
  ))
  selected_positions <- match(
    selection$selected_ids,
    paste0("leaf-", seq_len(n_leaves))
  )
  true_omitted <- sum(evidence[-selected_positions]) / sum(evidence)
  expect_true(
    selection$certified &&
      selection$epsilon_feasible &&
      selection$cost_optimality_met &&
      selection$cost_lower == exact$budget &&
      selection$cost_upper == exact$budget &&
      selection$cost_gap == 0 &&
      selected_cost == exact$budget &&
      true_omitted <= epsilon + 1e-12,
    sprintf("exact constrained tree search failed at trial %d", trial)
  )
}

# Selection performs no local sampling.  The optional downstream call is a
# separate operation and accounts for its own score evaluations.
expect_true(
  is.null(selection$local_fits) && is.null(selection$regional_weights),
  "the selection-only API unexpectedly ran a local approximation"
)
set.seed(20260816)
downstream <- approximate_selected_region_tree(
  search = selection,
  log_score_fn = log_score_fn,
  n_local_particles = 8L,
  local_mutation_steps = 0L
)
expect_true(
  length(downstream$local_fits) == length(selection$selected_ids) &&
    downstream$local_score_evaluations > 0 &&
    abs(sum(downstream$regional_weights) - 1) <= 1e-12,
  "the optional downstream approximation did not consume a certified selection"
)

expect_error(
  select_region_union_epsilon(
    log_evidence_lower = log(c(2, 3)),
    log_evidence_upper = log(c(1, 4)),
    costs = c(1L, 1L),
    epsilon = 0.1
  ),
  "the constrained selector accepted an inverted interval"
)
expect_error(
  select_certified_region_tree(
    tree = tree,
    log_score_fn = log_score_fn,
    evidence_evaluator = function(node) NULL,
    selection_mode = "epsilon"
  ),
  "epsilon selection accepted a missing epsilon"
)

expect_true(
  identical(
    rrs_registered_status(c(
      "certified",
      "frontier_exhausted",
      "evidence_interval_conflict",
      "uninformative_evidence",
      "budget_exhausted"
    )),
    c(
      "bounds-met",
      "tolerance-not-met",
      "interval-conflict",
      "uninformative-evidence",
      "budget-exhausted"
    )
  ),
  "raw search statuses do not map to the registered reporting vocabulary"
)

old_dp_limit <- getOption("rrs.max_dynamic_program_states")
options(rrs.max_dynamic_program_states = 3L)
expect_error(
  rrs_log_benefit_frontier_generic(
    log_benefit = log(c(1, 1)),
    costs = c(2L, 2L),
    eligible = c(TRUE, TRUE)
  ),
  "the exact cost dynamic program ignored its audited memory guard"
)
if (is.null(old_dp_limit)) {
  options(rrs.max_dynamic_program_states = NULL)
} else {
  options(rrs.max_dynamic_program_states = old_dp_limit)
}

# Exercise a genuinely adaptive multilevel trajectory with nonzero-width
# stochastic intervals. The target equals the declared Bernoulli proposal, so
# exact node masses are independently available while rho deliberately makes
# every reported stochastic interval nondegenerate.
multilevel_probability <- 0.25
multilevel_groups <- list(1:2, 3:4)
multilevel_tree <- new_group_state_tree(
  inclusion_probabilities = rep(multilevel_probability, 4L),
  groups = multilevel_groups,
  split_order = 1:2,
  max_depth = 2L,
  rho = 0.05,
  target_mutation = make_group_gibbs_target_mutation,
  state_cost = c(zero = 1L, one = 2L, multi = 2L),
  base_cost = 1L
)
multilevel_score <- function(gamma) {
  sum(
    gamma * log(multilevel_probability) +
      (1 - gamma) * log1p(-multilevel_probability)
  )
}
set.seed(20260817)
multilevel_fit <- fit_certified_region_selection(
  tree = multilevel_tree,
  log_score_fn = multilevel_score,
  delta = 0.05,
  n_pilot = 40L,
  n_paths = 200L,
  pilot_mutation_steps = 1L,
  production_mutation_steps = 1L,
  maximum_relative_radius = 0.05,
  selection_mode = "epsilon",
  epsilon = 0.20,
  cost_tolerance = 100,
  max_score_evaluations = 1000000L,
  max_actions = 100L
)
multilevel_supports <- as.matrix(expand.grid(rep(list(c(0L, 1L)), 4L)))
multilevel_probabilities <- exp(apply(
  multilevel_supports,
  1L,
  multilevel_score
))
evaluated_multilevel <- multilevel_fit$nodes[vapply(
  multilevel_fit$nodes,
  function(node) isTRUE(node$evaluated),
  logical(1)
)]
multilevel_interval_valid <- vapply(evaluated_multilevel, function(node) {
  inside <- is.finite(node$region$log_prob(multilevel_supports))
  exact_log_evidence <- log(sum(multilevel_probabilities[inside]))
  node$log_lower <= exact_log_evidence &&
    exact_log_evidence <= node$log_upper &&
    node$log_lower < node$log_upper
}, logical(1))
expect_true(
  multilevel_fit$certified &&
    identical(multilevel_fit$registered_status, "bounds-met") &&
    multilevel_fit$expansions >= 2L &&
    length(evaluated_multilevel) >= 1L &&
    all(multilevel_interval_valid) &&
    multilevel_fit$budget_overshoot == 0,
  "the stochastic multilevel constrained-search trajectory failed validation"
)

cat("constrained region-search unit tests passed\n")
