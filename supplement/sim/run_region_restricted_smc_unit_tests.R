source(file.path("sim", "src", "region_restricted_smc.R"))
source(file.path("sim", "src", "exact_gaussian_helpers.R"))
source(file.path("sim", "exact_bma.R"))

fail <- function(message) stop(message, call. = FALSE)
expect_true <- function(value, message) if (!isTRUE(value)) fail(message)
expect_near <- function(value, target, tolerance, message) {
  if (!is.finite(value) || abs(value - target) > tolerance) {
    fail(sprintf("%s: got %.8f, expected %.8f", message, value, target))
  }
}
expect_error <- function(expression, message) {
  failed <- inherits(try(force(expression), silent = TRUE), "try-error")
  if (!failed) {
    fail(message)
  }
}

supports <- rbind(
  c(0, 0),
  c(1, 0),
  c(0, 1),
  c(1, 1)
)

set.seed(20260725)
exchangeable_region <- new_enumerated_region(
  supports,
  cost = 1L,
  label = "exchangeable"
)
exchangeable_fit <- region_annealed_smc(
  exchangeable_region,
  log_score_fn = function(gamma) log(2),
  n_particles = 32L,
  mutation_steps = 1L
)
expect_near(
  exchangeable_fit$log_evidence,
  log(8),
  1e-12,
  "exact-exchangeability evidence identity failed"
)
expect_true(
  nrow(exchangeable_fit$diagnostics) == 1L,
  "exact exchangeability should use one temperature step"
)

score_values <- c(1, 2, 1.5, 3)
score_keys <- rrs_support_keys(supports)
score_function <- function(gamma) {
  log(score_values[match(paste(gamma, collapse = ","), score_keys)])
}

set.seed(20260725)
one_step_fit <- region_annealed_smc(
  exchangeable_region,
  log_score_fn = score_function,
  n_particles = 20000L,
  ess_fraction = 0.2,
  mutation_steps = 0L,
  phi_schedule = c(0, 1)
)
expect_near(
  one_step_fit$evidence,
  sum(score_values),
  0.08,
  "one-step regional evidence estimate failed"
)

heterogeneous_values <- exp(c(-4, -1, 1, 3))
heterogeneous_function <- function(gamma) {
  log(heterogeneous_values[
    match(paste(gamma, collapse = ","), score_keys)
  ])
}
set.seed(20260725)
annealed_fit <- region_annealed_smc(
  exchangeable_region,
  log_score_fn = heterogeneous_function,
  n_particles = 6000L,
  ess_fraction = 0.8,
  resample_fraction = 0.7,
  mutation_steps = 2L
)
expect_true(
  nrow(annealed_fit$diagnostics) > 1L,
  "heterogeneous region did not trigger annealing"
)
expect_near(
  annealed_fit$evidence,
  sum(heterogeneous_values),
  0.8,
  "annealed regional evidence estimate failed"
)
expect_true(
  all(diff(annealed_fit$diagnostics$phi) > 0) &&
    tail(annealed_fit$diagnostics$phi, 1) == 1,
  "adaptive temperatures are not strictly increasing to one"
)

set.seed(20260725)
ais_fit <- region_annealed_importance_sampling(
  exchangeable_region,
  log_score_fn = heterogeneous_function,
  phi_schedule = c(0, 0.25, 0.5, 0.75, 1),
  n_paths = 20000L,
  mutation_steps = 1L
)
expect_near(
  ais_fit$evidence,
  sum(heterogeneous_values),
  0.5,
  "fixed-schedule annealed importance estimate failed"
)
expect_true(
  identical(ais_fit$phi_schedule, c(0, 0.25, 0.5, 0.75, 1)),
  "production annealing did not preserve the declared schedule"
)

set.seed(20260725)
pilot_fixed <- pilot_fixed_region_evidence(
  exchangeable_region,
  log_score_fn = function(gamma) log(2),
  rho = 0,
  n_declared_nodes = 4L,
  n_pilot = 32L,
  n_paths = 32L
)
expect_near(
  pilot_fixed$evidence,
  8,
  1e-12,
  "pilot/fixed exchangeable evidence calculation failed"
)
expect_true(
  pilot_fixed$relative_radius == 0 &&
    pilot_fixed$log_lower == pilot_fixed$log_upper,
  "exact exchangeability should give a zero-width evidence interval"
)

set.seed(20260725)
bounded_pilot_fixed <- pilot_fixed_region_evidence(
  exchangeable_region,
  log_score_fn = score_function,
  rho = log(3),
  n_declared_nodes = 1L,
  delta = 0.05,
  n_pilot = 64L,
  n_paths = 4000L
)
expect_true(
  bounded_pilot_fixed$finite_interval &&
    bounded_pilot_fixed$log_lower <= log(sum(score_values)) &&
    bounded_pilot_fixed$log_upper >= log(sum(score_values)),
  "pilot/fixed evidence interval missed the exact finite-region evidence"
)

exact_epsilon_selection <- select_region_union_epsilon(
  log_evidence_lower = log(c(60, 30, 10)),
  costs = c(3L, 1L, 1L),
  epsilon = 0.1
)
expect_true(
  exact_epsilon_selection$feasible &&
    exact_epsilon_selection$budget == 4L &&
    identical(exact_epsilon_selection$selected, c(1L, 2L)),
  "minimum-cost epsilon selector did not recover the exact union"
)

brute_force_epsilon_cost <- function(log_lower,
                                     log_upper,
                                     costs,
                                     epsilon,
                                     mode = "conservative",
                                     eligible = rep(TRUE, length(costs))) {
  lower <- exp(log_lower)
  upper <- exp(log_upper)
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) return(Inf)
  best <- Inf
  for (mask in seq_len(2^length(eligible_ids) - 1L)) {
    chosen <- eligible_ids[
      as.logical(intToBits(mask)[seq_along(eligible_ids)])
    ]
    complement <- setdiff(seq_along(costs), chosen)
    omitted <- if (length(complement) == 0L) {
      0
    } else if (identical(mode, "conservative")) {
      sum(upper[complement]) /
        (sum(lower[chosen]) + sum(upper[complement]))
    } else {
      sum(lower[complement]) /
        (sum(upper[chosen]) + sum(lower[complement]))
    }
    if (omitted <= epsilon) best <- min(best, sum(costs[chosen]))
  }
  best
}

set.seed(20260813)
for (trial in seq_len(250L)) {
  n_regions <- sample.int(9L, 1L)
  truth <- exp(runif(n_regions, -4, 4))
  width <- exp(runif(n_regions, 0, 1.5))
  lower <- log(truth / width)
  upper <- log(truth * width)
  costs <- sample.int(5L, n_regions, replace = TRUE)
  epsilon <- runif(1L, 0.001, 0.95)
  eligible <- runif(n_regions) > 0.3
  for (mode in c("optimistic", "conservative")) {
    dynamic <- select_region_union_epsilon(
      log_evidence_lower = lower,
      log_evidence_upper = upper,
      costs = costs,
      epsilon = epsilon,
      mode = mode,
      eligible = if (identical(mode, "conservative")) eligible else NULL
    )
    brute <- brute_force_epsilon_cost(
      lower,
      upper,
      costs,
      epsilon,
      mode = mode,
      eligible = if (identical(mode, "conservative")) eligible else
        rep(TRUE, n_regions)
    )
    expect_true(
      identical(dynamic$feasible, is.finite(brute)) &&
        (is.infinite(brute) || dynamic$budget == brute),
      "epsilon dynamic program disagrees with exhaustive subset search"
    )
    if (identical(mode, "conservative") && dynamic$feasible) {
      exact_omitted <- sum(truth[-dynamic$selected]) / sum(truth)
      expect_true(
        exact_omitted <= epsilon + 1e-13,
        "conservative epsilon selection violated the exact tolerance"
      )
    }
  }
  exact_dynamic <- select_region_union_epsilon(
    log_evidence_lower = log(truth),
    costs = costs,
    epsilon = epsilon
  )
  exact_brute <- brute_force_epsilon_cost(
    log(truth), log(truth), costs, epsilon
  )
  expect_true(
    exact_dynamic$budget == exact_brute,
    "collapsed intervals did not recover the exact constrained optimum"
  )
}

preflight <- regional_evidence_preflight(
  rho = log(2),
  relative_error = 0.1,
  delta_node = 0.05,
  remaining_paths = 100L
)
expect_true(
  preflight$required_paths == ceiling(50 * log(40)) &&
    !preflight$feasible,
  "regional evidence preflight did not expose an inadequate budget"
)

all_supports_four <- as.matrix(expand.grid(
  rep(list(c(0L, 1L)), 4L)
))
group_region <- new_group_state_region(
  inclusion_probabilities = rep(0.25, 4L),
  groups = list(1:2, 3:4),
  states = c("one", "multi"),
  cost = 4L
)
group_probabilities <- exp(group_region$log_prob(all_supports_four))
expect_near(
  sum(group_probabilities),
  1,
  1e-12,
  "group-state proposal does not normalize"
)
set.seed(20260725)
group_samples <- group_region$sample(500L)
expect_true(
  all(rowSums(group_samples[, 1:2, drop = FALSE]) == 1L) &&
    all(rowSums(group_samples[, 3:4, drop = FALSE]) >= 2L),
  "group-state proposal violated its declared constraints"
)
set.seed(20260726)
group_mutation <- group_region$mutate(group_samples)
expect_true(
  identical(group_mutation$type, "group-refresh") &&
    all(rowSums(group_mutation$particles[, 1:2, drop = FALSE]) == 1L) &&
    all(rowSums(group_mutation$particles[, 3:4, drop = FALSE]) >= 2L),
  "group-refresh mutation violated its declared constraints"
)
set.seed(20260726)
group_mutation_fit <- region_annealed_smc(
  group_region,
  log_score_fn = function(gamma) group_region$log_prob(matrix(gamma, nrow = 1L)),
  n_particles = 64L,
  mutation_steps = 2L,
  phi_schedule = c(0, 0.5, 1)
)
expect_near(
  group_mutation_fit$evidence,
  1,
  1e-10,
  "group-refresh mutation failed on its invariant proposal target"
)
expect_true(
  all(group_mutation_fit$diagnostics$mutation_type == "group-refresh"),
  "regional SMC did not report the group-refresh mutation"
)

set.seed(20260726)
score_X <- matrix(rnorm(60), nrow = 20L, ncol = 3L)
score_y <- rnorm(20L)
score_object <- make_gaussian_support_score(
  score_X,
  score_y,
  inclusion_probability = 0.2,
  tau2 = 4,
  a0 = 1,
  b0 = 1
)
score_supports <- enumerate_supports(3L)
direct_scores <- apply(score_supports, 1L, function(gamma) {
  log_marginal_support(gamma, score_X, score_y, tau2 = 4, a0 = 1, b0 = 1) +
    sum(gamma) * log(0.2) +
    (3L - sum(gamma)) * log(0.8)
})
cached_scores <- apply(score_supports, 1L, score_object$score)
expect_true(
  max(abs(direct_scores - cached_scores)) < 1e-8,
  "cached sufficient-statistic support scores disagree with direct scores"
)
invisible(apply(score_supports, 1L, score_object$score))
score_diagnostics <- score_object$diagnostics()
expect_true(
  score_diagnostics[["evaluations"]] == nrow(score_supports) &&
    score_diagnostics[["cache_hits"]] == nrow(score_supports),
  "support-score memoization accounting failed"
)

orthogonal_case <- cb_block_orthogonal_case(
  p = 6L,
  n_train = 12L,
  n_test = 20L,
  rho = 0.95,
  seed = 20260726
)
known_score <- make_known_sigma_support_score(
  orthogonal_case$X_train,
  orthogonal_case$y_train,
  sigma = orthogonal_case$sigma,
  inclusion_probability = 0.2,
  tau2 = 2
)
known_supports <- enumerate_supports(6L)
known_log_scores <- apply(known_supports, 1L, known_score$score)
known_log_evidence <- cb_log_sum_exp(known_log_scores)
factorized_truth <- cb_factorized_group_truth(
  known_score,
  orthogonal_case$groups
)
expect_near(
  factorized_truth$log_evidence,
  known_log_evidence,
  1e-8,
  "block-factorized reference evidence disagrees with enumeration"
)
known_posterior <- exp(known_log_scores - known_log_evidence)
expect_true(
  max(abs(
    factorized_truth$pip -
      cb_pip(known_supports, known_posterior)
  )) < 1e-8,
  "block-factorized reference PIPs disagree with enumeration"
)

group_tree <- new_group_state_tree(
  inclusion_probabilities = rep(0.25, 4L),
  groups = list(1:2, 3:4),
  max_depth = 2L,
  rho = 0,
  materialize = TRUE
)
expect_true(
  group_tree$n_declared_nodes == 13L,
  "group-state tree has the wrong number of declared nodes"
)
lazy_group_tree <- new_group_state_tree(
  inclusion_probabilities = rep(0.25, 4L),
  groups = list(1:2, 3:4),
  max_depth = 2L,
  rho = 0
)
expect_true(
  inherits(lazy_group_tree, "rrs_lazy_region_tree") &&
    length(lazy_group_tree$nodes) == 1L &&
    lazy_group_tree$n_declared_nodes == 13L,
  "default group-state construction materialized the full tree"
)
group_leaf_ids <- names(group_tree$nodes)[vapply(
  group_tree$nodes,
  `[[`,
  logical(1),
  "terminal"
)]
group_memberships <- vapply(group_tree$nodes[group_leaf_ids], function(node) {
  is.finite(node$region$log_prob(all_supports_four))
}, logical(nrow(all_supports_four)))
expect_true(
  all(rowSums(group_memberships) == 1L),
  "terminal group-state regions do not partition the support space"
)
screening_score <- function(gamma) {
  sum(as.numeric(gamma) * c(log(2), log(3), log(5), log(7)))
}
screened_factory <- make_screened_group_region_factory(
  log_score_fn = screening_score,
  p = 4L,
  groups = list(1:2, 3:4)
)
screened_tree <- new_group_state_tree(
  inclusion_probabilities = rep(0.25, 4L),
  groups = list(1:2, 3:4),
  max_depth = 2L,
  rho = 0,
  region_factory = screened_factory$region_factory,
  materialize = TRUE
)
screened_leaf_ids <- names(screened_tree$nodes)[vapply(
  screened_tree$nodes,
  `[[`,
  logical(1),
  "terminal"
)]
screened_memberships <- vapply(
  screened_tree$nodes[screened_leaf_ids],
  function(node) is.finite(node$region$log_prob(all_supports_four)),
  logical(nrow(all_supports_four))
)
expect_true(
  all(rowSums(screened_memberships) == 1L),
  "screened terminal regions do not partition the support space"
)
for (leaf_id in screened_leaf_ids) {
  leaf <- screened_tree$nodes[[leaf_id]]$region
  probabilities <- exp(leaf$log_prob(all_supports_four))
  expect_near(
    sum(probabilities[is.finite(probabilities)]),
    1,
    1e-12,
    "screened region proposal is not normalized"
  )
  draws <- leaf$sample(20L)
  expect_true(
    all(is.finite(leaf$log_prob(draws))),
    "screened region sampler left its declared region"
  )
}
additive_rho <- make_additive_group_state_rho_bound(
  log_lower = matrix(
    c(-4, 0, -1, -4, -4, 0, -1, -4),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("free", "zero", "one", "multi"))
  ),
  log_upper = matrix(
    c(0, 0, -1, -4, 0, 0, -1, -4),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("free", "zero", "one", "multi"))
  )
)
expect_near(
  additive_rho(c("free", "one")),
  4,
  1e-12,
  "additive group-state oscillation helper returned the wrong bound"
)

tree_supports <- rbind(c(0L, 0L), c(1L, 0L), c(0L, 1L))
tree_scores <- c(100, 1, 1)
tree_keys <- rrs_support_keys(tree_supports)
tree_score_function <- function(gamma) {
  log(tree_scores[match(paste(gamma, collapse = ","), tree_keys)])
}
root_region <- new_enumerated_region(
  tree_supports,
  cost = 1L,
  label = "root"
)
leaf_regions <- lapply(seq_len(nrow(tree_supports)), function(i) {
  new_enumerated_region(
    tree_supports[i, , drop = FALSE],
    cost = 1L,
    label = paste0("leaf-", i)
  )
})
tree_nodes <- c(
  list(new_region_tree_node(
    id = "root",
    region = root_region,
    children = paste0("leaf-", 1:3),
    rho = log(100),
    cost_lower = 1L,
    terminal = FALSE
  )),
  lapply(seq_along(leaf_regions), function(i) {
    new_region_tree_node(
      id = paste0("leaf-", i),
      region = leaf_regions[[i]],
      rho = 0,
      cost_lower = 1L,
      terminal = TRUE
    )
  })
)
tree <- new_region_tree(tree_nodes, "root")
expect_error(
  new_region_tree(c(tree_nodes, list(tree_nodes[[2L]])), "root"),
  "region-tree validation did not reject a duplicate node id"
)
exact_tree_evaluator <- function(node) {
  values <- exp(rrs_evaluate_log_score(
    tree_score_function,
    node$region$supports
  ))
  log_evidence <- log(sum(values))
  list(
    log_evidence = log_evidence,
    log_lower = log_evidence,
    log_upper = log_evidence,
    score_evaluations = nrow(node$region$supports),
    fit = list(exact = TRUE)
  )
}
epsilon_tree_fit <- certified_tree_region_search(
  tree = tree,
  log_score_fn = tree_score_function,
  evidence_evaluator = exact_tree_evaluator,
  epsilon = 0.03,
  cost_tolerance = 0
)
expect_true(
  epsilon_tree_fit$certified && epsilon_tree_fit$epsilon_feasible &&
    epsilon_tree_fit$cost_optimality_met &&
    epsilon_tree_fit$cost_lower == epsilon_tree_fit$cost_upper &&
    epsilon_tree_fit$cost_gap == 0 &&
    epsilon_tree_fit$omitted_mass_bound <= 0.03,
  "end-to-end constrained tree search did not certify exact cost"
)

conflict_nodes <- c(
  list(new_region_tree_node(
    id = "root",
    region = root_region,
    children = paste0("leaf-", 1:3),
    rho = log(100),
    cost_lower = 1L,
    terminal = FALSE
  )),
  lapply(seq_along(leaf_regions), function(i) {
    new_region_tree_node(
      id = paste0("leaf-", i),
      region = leaf_regions[[i]],
      rho = 0.1,
      cost_lower = 1L,
      terminal = TRUE
    )
  })
)
conflict_evaluator <- function(node) {
  list(
    log_evidence = 10.5,
    log_lower = 10,
    log_upper = 11,
    score_evaluations = 2L,
    fit = list(conflict_test = TRUE)
  )
}
conflict_fit <- certified_tree_region_search(
  tree = new_region_tree(conflict_nodes, "root"),
  log_score_fn = tree_score_function,
  evidence_evaluator = conflict_evaluator,
  epsilon = 0.03,
  cost_tolerance = 0
)
expect_true(
  identical(conflict_fit$status, "evidence_interval_conflict") &&
    !conflict_fit$certified &&
    !is.na(conflict_fit$conflict_node_id),
  "a conflicting stochastic evidence interval was not reported safely"
)

expect_near(
  gaussian_swap_log_score_bound(
    column_distance = 0,
    column_norm_bound = 2,
    response_norm = 3,
    sigma = 1,
    slab_scale = 1,
    log_prior_ratio_bound = 0.2
  ),
  0.2,
  1e-12,
  "Gaussian swap bound has the wrong zero-distance limit"
)

prior_probabilities <- c(0.25, 0.25)
prior_log_score <- function(gamma) {
  sum(
    gamma * log(prior_probabilities) +
      (1 - gamma) * log1p(-prior_probabilities)
  )
}
small_group_tree <- new_group_state_tree(
  inclusion_probabilities = prior_probabilities,
  groups = list(1:2),
  max_depth = 1L,
  rho = 0
)
set.seed(20260724)
iid_evidence <- iid_bounded_region_evidence(
  region = small_group_tree$nodes[[small_group_tree$root_id]]$region,
  log_score_fn = prior_log_score,
  rho = 0,
  n_declared_nodes = small_group_tree$n_declared_nodes,
  delta = 0.05,
  n_paths = 16L
)
expect_near(
  iid_evidence$log_evidence,
  0,
  1e-12,
  "independent bounded evidence is not exact for target-equals-proposal"
)
expect_true(
  iid_evidence$log_lower == 0 && iid_evidence$log_upper == 0 &&
    iid_evidence$score_evaluations == 16L,
  "independent bounded evidence returned the wrong zero-range certificate"
)
set.seed(20260725)
certified_fit <- fit_certified_region_tree(
  tree = small_group_tree,
  log_score_fn = prior_log_score,
  epsilon = 0.2,
  cost_tolerance = 0,
  n_pilot = 16L,
  n_paths = 16L,
  n_local_particles = 16L,
  max_actions = 10L
)
expect_true(
  certified_fit$certified &&
    length(certified_fit$local_fits) >= 1L &&
    certified_fit$total_score_evaluations ==
      certified_fit$bound_score_evaluations +
      certified_fit$evidence_score_evaluations +
      certified_fit$local_score_evaluations,
  "end-to-end certified tree fit or cost accounting failed"
)

set.seed(20260813)
selection_only_fit <- fit_certified_region_selection(
  tree = small_group_tree,
  log_score_fn = prior_log_score,
  n_pilot = 16L,
  n_paths = 16L,
  max_actions = 10L,
  epsilon = 0.2,
  cost_tolerance = 0
)
expect_true(
  selection_only_fit$certified && selection_only_fit$epsilon_feasible &&
    selection_only_fit$cost_gap == 0 &&
    is.null(selection_only_fit$local_fits),
  "selection-only constrained wrapper did not remain independent of local SMC"
)

set.seed(20260725)
uncertified_fit <- fit_certified_region_tree(
  tree = small_group_tree,
  log_score_fn = prior_log_score,
  epsilon = 0.01,
  cost_tolerance = 0,
  n_pilot = 16L,
  n_paths = 16L,
  n_local_particles = 16L,
  max_score_evaluations = 1L,
  max_actions = 1L
)
expect_true(
  !uncertified_fit$certified &&
    identical(uncertified_fit$selected_ids, character(0)) &&
    length(uncertified_fit$local_fits) == 0L &&
    uncertified_fit$local_score_evaluations == 0,
  "local SMC ran before the tree search earned a certificate"
)

test_schedule <- c(0, 0.02, 0.12, 1)
test_diagnostics <- data.frame(
  phi_previous = head(test_schedule, -1L),
  phi = tail(test_schedule, -1L),
  ess = c(90, 20, 80)
)
refined_schedule <- rrs_refine_phi_schedule(
  test_schedule,
  test_diagnostics,
  n_particles = 100L,
  ess_fraction = 0.5,
  max_increment = 0.2
)
expect_true(
  any(abs(refined_schedule - 0.07) < 1e-12) &&
    max(diff(refined_schedule)) <= 0.2 + 1e-12 &&
    !any(abs(refined_schedule - 0.01) < 1e-12),
  "fixed-schedule refinement did not split exactly the failing intervals"
)

pilot_particles <- rbind(
  c(0, 0, 0, 0),
  c(0, 0, 1, 0),
  c(0, 0, 0, 1),
  c(1, 0, 1, 0)
)
pilot_weights <- c(0.45, 0.25, 0.20, 0.10)
mass_particles <- rbind(
  c(0, 0, 0, 0),
  c(0, 0, 1, 0),
  c(1, 0, 0, 0),
  c(1, 0, 1, 0)
)
mass_weights <- c(0.40, 0.30, 0.20, 0.10)
pilot_cells <- cb_select_pilot_cells(
  pilot_particles,
  pilot_weights,
  groups = list(1:2, 3:4),
  mass_particles = mass_particles,
  mass_weights = mass_weights,
  n_split_groups = 2L,
  n_cells = 4L,
  target_mass = 0.85,
  minimum_confidence = 0.8,
  split_rule = "confidence"
)
expect_true(
  identical(pilot_cells$split_groups, 1L) &&
    pilot_cells$selected_pilot_mass >= 0.85 &&
    length(pilot_cells$states) == 2L,
  "independent mass-targeted pilot-cell selection failed"
)

unrestricted_cells <- cb_select_pilot_cells(
  pilot_particles,
  pilot_weights,
  groups = list(1:2, 3:4),
  n_split_groups = 2L,
  n_cells = 4L,
  target_mass = 0.85,
  minimum_confidence = 0.99,
  split_rule = "confidence"
)
expect_true(
  length(unrestricted_cells$split_groups) == 0L &&
    length(unrestricted_cells$states) == 1L &&
    unrestricted_cells$selected_pilot_mass == 1,
  "pilot-cell fallback to the unrestricted region failed"
)
unrestricted_keys <- cb_group_state_keys(
  supports,
  unrestricted_cells$groups
)
unrestricted_selected_key <- vapply(
  unrestricted_cells$states,
  paste,
  collapse = "/",
  character(1)
)
expect_true(
  all(unrestricted_keys %in% unrestricted_selected_key),
  "the fallback cell does not retain every support"
)

cell_regions <- cb_make_cell_regions(
  inclusion_probabilities = c(0.25, 0.25),
  all_groups = list(1:2),
  split_groups = 1L,
  states = list("zero", "one"),
  new_group_state_region_fn = new_group_state_region
)
cell_membership <- vapply(cell_regions, function(region) {
  is.finite(region$log_prob(supports))
}, logical(nrow(supports)))
expect_true(
  all(rowSums(cell_membership) <= 1L) &&
    all(rowSums(cell_membership)[rowSums(supports) <= 1L] == 1L),
  "selected regional cells are not a disjoint partition of their union"
)

# Independently verify the unknown-variance Gaussian BMA score. This
# calculation works in observation space and therefore does not reuse the
# implementation's Woodbury formula.
unknown_x <- matrix(c(
  -1.2, 0.3, 1.1,
  -0.4, 0.8, -0.5,
  0.1, -1.0, 0.7,
  0.8, 0.4, -0.9,
  1.3, -0.5, 0.2
), nrow = 5L, byrow = TRUE)
unknown_y <- c(-0.7, 0.2, -0.1, 0.9, -0.3)
unknown_q <- 0.2
unknown_tau2 <- 3
unknown_a0 <- 1.5
unknown_b0 <- 0.8
unknown_score <- make_gaussian_support_score(
  unknown_x,
  unknown_y,
  inclusion_probability = unknown_q,
  tau2 = unknown_tau2,
  a0 = unknown_a0,
  b0 = unknown_b0,
  memoize = FALSE
)
unknown_supports <- as.matrix(expand.grid(rep(list(c(0L, 1L)), 3L)))
unknown_direct <- apply(unknown_supports, 1L, function(gamma) {
  active <- which(gamma == 1L)
  covariance <- diag(nrow(unknown_x))
  if (length(active) > 0L) {
    active_x <- unknown_x[, active, drop = FALSE]
    covariance <- covariance + unknown_tau2 * tcrossprod(active_x)
  }
  factor <- chol(covariance)
  log_determinant <- 2 * sum(log(diag(factor)))
  quadratic <- sum(unknown_y * backsolve(
    factor,
    forwardsolve(t(factor), unknown_y)
  ))
  a_n <- unknown_a0 + length(unknown_y) / 2
  log_marginal <- lgamma(a_n) - lgamma(unknown_a0) -
    length(unknown_y) / 2 * log(2 * pi) +
    unknown_a0 * log(unknown_b0) -
    0.5 * log_determinant -
    a_n * log(unknown_b0 + 0.5 * quadratic)
  size <- sum(gamma)
  log_marginal + size * log(unknown_q) +
    (length(gamma) - size) * log1p(-unknown_q)
})
unknown_implemented <- apply(
  unknown_supports,
  1L,
  unknown_score$score
)
expect_true(
  max(abs(unknown_direct - unknown_implemented)) < 1e-9,
  "the unknown-variance Gaussian BMA score failed an independent formula check"
)

cell_prior_score <- function(gamma) {
  sum(gamma * log(0.25) + (1 - gamma) * log(0.75))
}
set.seed(20260726)
cell_fits <- lapply(cell_regions, function(region) {
  region_annealed_smc(
    region,
    log_score_fn = cell_prior_score,
    n_particles = 32L,
    mutation_steps = 1L,
    phi_schedule = c(0, 1)
  )
})
cell_combined <- cb_combine_region_fits(cell_fits)
expected_cell_weights <- c(0.75^2, 2 * 0.25 * 0.75)
expected_cell_weights <- expected_cell_weights / sum(expected_cell_weights)
expect_true(
  max(abs(
    cell_combined$region_weights - expected_cell_weights
  )) < 1e-12,
  "regional evidence weights do not recover the posterior conditional on the union"
)

# A finite search ceiling is enforced before every score batch, including
# batches inside a partially completed stochastic evidence action.
budget_region <- new_enumerated_region(
  supports = rbind(c(0L), c(1L)),
  cost = 1L,
  label = "budget-root"
)
budget_node <- new_region_tree_node(
  id = "budget-root",
  region = budget_region,
  rho = 2,
  terminal = TRUE
)
budget_tree <- new_region_tree(
  stats::setNames(list(budget_node), "budget-root"),
  "budget-root"
)
budget_score <- function(gamma) -0.25 * sum(gamma)
budget_evaluator <- make_pilot_fixed_evidence_evaluator(
  log_score_fn = budget_score,
  n_declared_nodes = 1L,
  n_pilot = 2L,
  n_paths = 2L,
  pilot_mutation_steps = 0L,
  production_mutation_steps = 0L
)
budget_fit <- certified_tree_region_search(
  tree = budget_tree,
  log_score_fn = budget_score,
  evidence_evaluator = budget_evaluator,
  epsilon = 0.5,
  cost_tolerance = 0,
  max_score_evaluations = 2L
)
expect_true(
  identical(budget_fit$status, "budget_exhausted") &&
    budget_fit$total_score_evaluations == 2L &&
    budget_fit$budget_overshoot == 0 &&
    isTRUE(budget_fit$budget_checked_before_score_batches),
  "the evidence action exceeded or misreported its atomic score budget"
)

# The production-path preflight is now an enforced search outcome rather than
# a callable diagnostic that applications can accidentally bypass.
preflight_evaluator <- make_pilot_fixed_evidence_evaluator(
  log_score_fn = budget_score,
  n_declared_nodes = 1L,
  delta = 0.05,
  n_pilot = 2L,
  n_paths = 2L,
  pilot_mutation_steps = 0L,
  production_mutation_steps = 0L,
  maximum_relative_radius = 0.01
)
preflight_fit <- certified_tree_region_search(
  tree = budget_tree,
  log_score_fn = budget_score,
  evidence_evaluator = preflight_evaluator,
  epsilon = 0.5,
  cost_tolerance = 0,
  max_score_evaluations = 100L
)
expect_true(
  identical(preflight_fit$status, "certified") &&
    identical(preflight_fit$registered_status, "bounds-met") &&
    preflight_fit$total_score_evaluations == 1L &&
    isTRUE(preflight_fit$nodes[["budget-root"]]$preflight_blocked) &&
    identical(
      preflight_fit$nodes[["budget-root"]]$evidence_fit$method,
      "deterministic_structural_fallback"
    ),
  "an infeasible evidence preflight did not use the valid structural fallback"
)

accounting <- rrs_common_accounting(
  score_diagnostics = c(calls = 12, evaluations = 9, cache_hits = 3),
  search = list(
    bound_score_evaluations = 2,
    evidence_score_evaluations = 7,
    total_score_evaluations = 9
  ),
  downstream_score_requests = 3,
  elapsed_seconds = 0.25,
  peak_memory_bytes = 1024
)
expect_true(
  accounting$raw_score_requests == 12 &&
    accounting$distinct_score_evaluations == 9 &&
    accounting$selection_score_requests == 9 &&
    accounting$downstream_score_requests == 3 &&
    accounting$peak_memory_bytes == 1024,
  "the common computation-accounting schema changed or lost a stage"
)

message("region-restricted SMC unit tests passed")
