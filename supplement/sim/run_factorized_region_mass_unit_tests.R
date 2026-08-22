source(file.path("sim", "src", "region_restricted_smc.R"))
source(file.path("sim", "src", "exact_gaussian_helpers.R"))
source(file.path("sim", "exact_bma.R"))

expect <- function(condition, message) if (!isTRUE(condition)) stop(message, call. = FALSE)
lookup <- cb_make_precomputed_support_score(seq_len(8L), 3L)
lookup_values <- c(
  lookup$score(c(0L, 0L, 0L)),
  lookup$score(c(1L, 1L, 1L)),
  lookup$score(c(0L, 0L, 0L))
)
lookup_diagnostics <- lookup$diagnostics()
expect(
  identical(lookup_values, c(1, 8, 1)) &&
    unname(lookup_diagnostics[["calls"]]) == 3L &&
    unname(lookup_diagnostics[["evaluations"]]) == 2L &&
    unname(lookup_diagnostics[["cache_hits"]]) == 1L,
  "The precomputed support oracle misindexed scores or accounting."
)
lookup$reset(clear_cache = TRUE)
lookup_batch <- rbind(
  c(0L, 0L, 0L), c(0L, 0L, 0L),
  c(1L, 0L, 1L), c(1L, 0L, 1L), c(1L, 1L, 1L)
)
lookup_batch_values <- attr(lookup$score, "matrix_score")(lookup_batch)
lookup_batch_diagnostics <- lookup$diagnostics()
expect(
  identical(lookup_batch_values, c(1, 1, 6, 6, 8)) &&
    unname(lookup_batch_diagnostics[["calls"]]) == 5L &&
    unname(lookup_batch_diagnostics[["evaluations"]]) == 3L &&
    unname(lookup_batch_diagnostics[["cache_hits"]]) == 2L,
  "The precomputed support oracle batch path misindexed scores or accounting."
)
lookup$reset(clear_cache = FALSE)
lookup$score(c(1L, 1L, 1L))
expect(
  unname(lookup$diagnostics()[["cache_hits"]]) == 1L,
  "The precomputed support oracle did not preserve its cache on request."
)
lookup$reset(clear_cache = TRUE)
lookup$score(c(1L, 1L, 1L))
expect(
  unname(lookup$diagnostics()[["evaluations"]]) == 1L,
  "The precomputed support oracle did not clear its cache."
)
set.seed(20260813L)
data <- cb_block_orthogonal_case(
  p = 6L, n_train = 30L, n_test = 10L, rho = 0.85,
  signal = 1.4, sigma = 1, seed = 20260813L, heterogeneous = TRUE
)
base <- make_known_sigma_support_score(
  data$X_train, data$y_train, sigma = 1,
  inclusion_probability = 0.1, tau2 = 2
)
score <- make_factorized_known_sigma_score(base, data$groups)
truth <- cb_factorized_group_truth(score, data$groups)
score$reset()
duplicate_batch <- rbind(integer(6L), integer(6L), c(1L, 0L, 0L, 0L, 0L, 0L))
batch_value <- attr(score$score, "matrix_score")(duplicate_batch)
batch_diagnostics <- score$diagnostics()
expect(
  length(batch_value) == 3L && batch_value[[1L]] == batch_value[[2L]] &&
    unname(batch_diagnostics[["calls"]]) == 3L &&
    unname(batch_diagnostics[["evaluations"]]) == 2L &&
    unname(batch_diagnostics[["cache_hits"]]) == 1L,
  "The factorized matrix score does not honor memoized distinct accounting."
)
supports <- enumerate_supports(6L)
log_scores <- rrs_evaluate_log_score(score$score, supports)
weights <- exp(log_scores - cb_log_sum_exp(log_scores))
aligned <- data$groups
crossed <- list(c(1L, 4L), c(3L, 6L), c(5L, 2L))
for (declared in list(aligned, crossed)) {
  score$reset()
  screened <- make_screened_group_region_factory(score$score, 6L, declared)
  for (trial in seq_len(100L)) {
    states <- sample(c("free", "zero", "one", "multi"), 3L, replace = TRUE)
    exact_inside <- rep(TRUE, nrow(supports))
    for (group in seq_along(declared)) {
      count <- rowSums(supports[, declared[[group]], drop = FALSE])
      exact_inside <- exact_inside & switch(
        states[[group]], free = TRUE, zero = count == 0L,
        one = count == 1L, multi = count >= 2L
      )
    }
    exact <- sum(weights[exact_inside])
    factorized <- cb_factorized_region_mass(
      truth, aligned, declared, states
    )
    expect(
      abs(exact - factorized) < 1e-12,
      sprintf("Factorized regional mass failed at trial %d.", trial)
    )
    fixed_variable <- sample.int(6L, 1L)
    fixed_value <- sample(0:1, 1L)
    fixed <- setNames(fixed_value, fixed_variable)
    exact_fixed <- sum(weights[exact_inside & supports[, fixed_variable] == fixed_value])
    factorized_fixed <- cb_factorized_region_mass(
      truth, aligned, declared, states, fixed_bits = fixed
    )
    expect(
      abs(exact_fixed - factorized_fixed) < 1e-12,
      sprintf("Factorized fixed-bit regional mass failed at trial %d.", trial)
    )
    region <- screened$region_factory(states, 1L, "test")
    inside <- is.finite(region$log_prob(supports))
    exact_log_g <- log_scores[inside] - region$log_prob(supports[inside, , drop = FALSE])
    bounds <- cb_factorized_region_log_g_bounds(
      score, truth, aligned, screened$configuration_tables, declared, states
    )
    expect(
      abs(bounds[[1L]] - min(exact_log_g)) < 1e-10 &&
        abs(bounds[[2L]] - max(exact_log_g)) < 1e-10,
      sprintf("Factorized log-g bounds failed at trial %d.", trial)
    )
  }
}
cat("factorized regional-mass unit tests passed\n")
