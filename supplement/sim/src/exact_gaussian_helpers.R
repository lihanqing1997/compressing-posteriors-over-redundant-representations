# Gaussian and exact-table helpers for the constrained exact study.
#
# This focused release module contains only the score, synthetic-design, and
# factorized-reference utilities required by the released production runner and
# its unit tests. The broader internal benchmark helper is intentionally not a
# release dependency.

cb_log_sum_exp <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(-Inf)
  }
  maximum <- max(x[finite])
  maximum + log(sum(exp(x[finite] - maximum)))
}

cb_support_key <- function(gamma) {
  paste0(as.integer(gamma), collapse = "")
}

cb_make_precomputed_support_score <- function(log_scores, p, metadata = list()) {
  log_scores <- as.numeric(log_scores)
  p <- as.integer(p)
  expected <- 2^p
  if (length(log_scores) != expected || any(!is.finite(log_scores))) {
    stop("The precomputed support-score table is invalid.", call. = FALSE)
  }
  powers <- 2^(seq_len(p) - 1L)
  seen <- rep(FALSE, expected)
  counters <- new.env(parent = emptyenv())
  counters$calls <- 0L
  counters$evaluations <- 0L
  counters$cache_hits <- 0L
  score <- function(gamma) {
    gamma <- as.integer(gamma)
    if (length(gamma) != p || any(!gamma %in% c(0L, 1L))) {
      stop("A support score received an invalid support.", call. = FALSE)
    }
    counters$calls <- counters$calls + 1L
    index <- 1L + sum(gamma * powers)
    if (seen[[index]]) {
      counters$cache_hits <- counters$cache_hits + 1L
    } else {
      seen[[index]] <<- TRUE
      counters$evaluations <- counters$evaluations + 1L
    }
    log_scores[[index]]
  }
  attr(score, "matrix_score") <- function(particles) {
    particles <- as.matrix(particles)
    if (ncol(particles) != p || any(!particles %in% c(0L, 1L))) {
      stop("A support-score batch contained an invalid support.", call. = FALSE)
    }
    indices <- 1L + as.numeric(particles %*% powers)
    newly_seen <- unique(indices[!seen[indices]])
    counters$calls <- counters$calls + nrow(particles)
    counters$evaluations <- counters$evaluations + length(newly_seen)
    counters$cache_hits <- counters$cache_hits +
      nrow(particles) - length(newly_seen)
    if (length(newly_seen)) seen[newly_seen] <<- TRUE
    log_scores[indices]
  }
  diagnostics <- function() c(
    calls = counters$calls,
    evaluations = counters$evaluations,
    cache_hits = counters$cache_hits,
    cached_supports = sum(seen)
  )
  reset <- function(clear_cache = TRUE) {
    counters$calls <- 0L
    counters$evaluations <- 0L
    counters$cache_hits <- 0L
    if (clear_cache) seen[] <<- FALSE
    invisible(NULL)
  }
  c(
    list(
      score = score,
      diagnostics = diagnostics,
      reset = reset,
      p = p,
      score_oracle = "precomputed_exact_table"
    ),
    metadata
  )
}

cb_safe_chol <- function(A, jitter = 1e-10, max_tries = 8L) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2
  for (attempt in seq_len(max_tries)) {
    ridge <- jitter * 10^(attempt - 1L)
    factor <- tryCatch(
      chol(A + diag(ridge, nrow(A))),
      error = function(error) NULL
    )
    if (!is.null(factor)) {
      return(factor)
    }
  }
  eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  ridge <- max(jitter, -min(eigenvalues) + jitter)
  chol(A + diag(ridge, nrow(A)))
}

make_gaussian_support_score <- function(X,
                                        y,
                                        inclusion_probability = 0.05,
                                        tau2 = 4,
                                        a0 = 1,
                                        b0 = 1,
                                        memoize = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n ||
      !(inclusion_probability > 0 && inclusion_probability < 1) ||
      tau2 <= 0 ||
      a0 <= 0 ||
      b0 <= 0) {
    stop("Invalid Gaussian support-score inputs.", call. = FALSE)
  }
  gram <- crossprod(X)
  xty <- as.numeric(crossprod(X, y))
  yty <- sum(y^2)
  a_n <- a0 + n / 2
  log_constant <- lgamma(a_n) - lgamma(a0) -
    (n / 2) * log(2 * pi) + a0 * log(b0)
  log_prior_zero <- log1p(-inclusion_probability)
  log_prior_one <- log(inclusion_probability)
  cache <- new.env(hash = TRUE, parent = emptyenv())
  counters <- new.env(parent = emptyenv())
  counters$calls <- 0L
  counters$evaluations <- 0L
  counters$cache_hits <- 0L

  log_marginal <- function(gamma) {
    gamma <- as.integer(gamma)
    if (length(gamma) != p || any(!gamma %in% c(0L, 1L))) {
      stop("A support score received an invalid support.", call. = FALSE)
    }
    active <- which(gamma == 1L)
    q <- length(active)
    if (q == 0L) {
      log_determinant <- 0
      quadratic <- yty
    } else {
      system <- diag(q) + tau2 * gram[active, active, drop = FALSE]
      factor <- cb_safe_chol(system)
      log_determinant <- 2 * sum(log(diag(factor)))
      active_xty <- xty[active]
      solution <- backsolve(
        factor,
        forwardsolve(t(factor), active_xty)
      )
      quadratic <- max(yty - tau2 * sum(active_xty * solution), 0)
    }
    log_constant -
      0.5 * log_determinant -
      a_n * log(b0 + 0.5 * quadratic)
  }

  score <- function(gamma) {
    counters$calls <- counters$calls + 1L
    key <- cb_support_key(gamma)
    if (memoize && exists(key, envir = cache, inherits = FALSE)) {
      counters$cache_hits <- counters$cache_hits + 1L
      return(get(key, envir = cache, inherits = FALSE))
    }
    gamma <- as.integer(gamma)
    size <- sum(gamma)
    value <- log_marginal(gamma) +
      size * log_prior_one +
      (p - size) * log_prior_zero
    counters$evaluations <- counters$evaluations + 1L
    if (memoize) {
      assign(key, value, envir = cache)
    }
    value
  }

  diagnostics <- function() {
    c(
      calls = counters$calls,
      evaluations = counters$evaluations,
      cache_hits = counters$cache_hits,
      cached_supports = length(ls(cache, all.names = TRUE))
    )
  }
  reset <- function(clear_cache = TRUE) {
    counters$calls <- 0L
    counters$evaluations <- 0L
    counters$cache_hits <- 0L
    if (clear_cache) {
      rm(list = ls(cache, all.names = TRUE), envir = cache)
    }
    invisible(NULL)
  }
  log_marginal_upper <- log_constant - a_n * log(b0)

  structure(
    list(
      score = score,
      log_marginal = log_marginal,
      diagnostics = diagnostics,
      reset = reset,
      gram = gram,
      xty = xty,
      yty = yty,
      n = n,
      p = p,
      inclusion_probability = inclusion_probability,
      tau2 = tau2,
      a0 = a0,
      b0 = b0,
      a_n = a_n,
      log_constant = log_constant,
      log_marginal_upper = log_marginal_upper
    ),
    class = "cb_gaussian_support_score"
  )
}

make_known_sigma_support_score <- function(X,
                                           y,
                                           sigma = 1,
                                           inclusion_probability = 0.05,
                                           tau2 = 4,
                                           memoize = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n ||
      sigma <= 0 ||
      !(inclusion_probability > 0 && inclusion_probability < 1) ||
      tau2 <= 0) {
    stop("Invalid known-variance support-score inputs.", call. = FALSE)
  }
  gram <- crossprod(X)
  xty <- as.numeric(crossprod(X, y))
  yty <- sum(y^2)
  log_constant <- -(n / 2) * log(2 * pi * sigma^2) -
    yty / (2 * sigma^2)
  log_prior_zero <- log1p(-inclusion_probability)
  log_prior_one <- log(inclusion_probability)
  cache <- new.env(hash = TRUE, parent = emptyenv())
  counters <- new.env(parent = emptyenv())
  counters$calls <- 0L
  counters$evaluations <- 0L
  counters$cache_hits <- 0L

  log_likelihood_increment <- function(gamma) {
    gamma <- as.integer(gamma)
    active <- which(gamma == 1L)
    q <- length(active)
    if (q == 0L) {
      return(0)
    }
    system <- diag(q) + tau2 * gram[active, active, drop = FALSE]
    factor <- cb_safe_chol(system)
    log_determinant <- 2 * sum(log(diag(factor)))
    active_xty <- xty[active]
    solution <- backsolve(
      factor,
      forwardsolve(t(factor), active_xty)
    )
    -0.5 * log_determinant +
      (tau2 / (2 * sigma^2)) * sum(active_xty * solution)
  }

  score <- function(gamma) {
    counters$calls <- counters$calls + 1L
    key <- cb_support_key(gamma)
    if (memoize && exists(key, envir = cache, inherits = FALSE)) {
      counters$cache_hits <- counters$cache_hits + 1L
      return(get(key, envir = cache, inherits = FALSE))
    }
    gamma <- as.integer(gamma)
    if (length(gamma) != p || any(!gamma %in% c(0L, 1L))) {
      stop("A support score received an invalid support.", call. = FALSE)
    }
    size <- sum(gamma)
    value <- log_constant +
      log_likelihood_increment(gamma) +
      size * log_prior_one +
      (p - size) * log_prior_zero
    counters$evaluations <- counters$evaluations + 1L
    if (memoize) {
      assign(key, value, envir = cache)
    }
    value
  }
  diagnostics <- function() {
    c(
      calls = counters$calls,
      evaluations = counters$evaluations,
      cache_hits = counters$cache_hits,
      cached_supports = length(ls(cache, all.names = TRUE))
    )
  }
  reset <- function(clear_cache = TRUE) {
    counters$calls <- 0L
    counters$evaluations <- 0L
    counters$cache_hits <- 0L
    if (clear_cache) {
      rm(list = ls(cache, all.names = TRUE), envir = cache)
    }
    invisible(NULL)
  }

  structure(
    list(
      score = score,
      log_likelihood_increment = log_likelihood_increment,
      diagnostics = diagnostics,
      reset = reset,
      gram = gram,
      xty = xty,
      yty = yty,
      n = n,
      p = p,
      sigma = sigma,
      inclusion_probability = inclusion_probability,
      tau2 = tau2,
      log_constant = log_constant
    ),
    class = "cb_known_sigma_support_score"
  )
}

make_factorized_known_sigma_score <- function(base_score, groups) {
  if (!inherits(base_score, "cb_known_sigma_support_score")) {
    stop("base_score must be a known-variance Gaussian score.", call. = FALSE)
  }
  p <- base_score$p
  if (!identical(sort(unlist(groups, use.names = FALSE)), seq_len(p))) {
    stop("groups must partition the support coordinates.", call. = FALSE)
  }
  local_tables <- lapply(groups, function(indices) {
    supports <- as.matrix(expand.grid(
      rep(list(c(0L, 1L)), length(indices))
    ))
    values <- apply(supports, 1L, function(local) {
      support <- integer(p)
      support[indices] <- local
      base_score$log_likelihood_increment(support)
    })
    keys <- apply(supports, 1L, paste0, collapse = "")
    stats::setNames(values, keys)
  })
  cache <- new.env(hash = TRUE, parent = emptyenv())
  counters <- new.env(parent = emptyenv())
  counters$calls <- 0L
  counters$evaluations <- 0L
  counters$cache_hits <- 0L
  theta <- base_score$inclusion_probability

  log_likelihood_increment <- function(gamma) {
    gamma <- as.integer(gamma)
    sum(vapply(seq_along(groups), function(group) {
      key <- paste0(gamma[groups[[group]]], collapse = "")
      unname(local_tables[[group]][key])
    }, numeric(1)))
  }
  score <- function(gamma) {
    counters$calls <- counters$calls + 1L
    gamma <- as.integer(gamma)
    key <- cb_support_key(gamma)
    if (exists(key, envir = cache, inherits = FALSE)) {
      counters$cache_hits <- counters$cache_hits + 1L
      return(get(key, envir = cache, inherits = FALSE))
    }
    size <- sum(gamma)
    value <- base_score$log_constant +
      log_likelihood_increment(gamma) +
      size * log(theta) +
      (p - size) * log1p(-theta)
    counters$evaluations <- counters$evaluations + 1L
    assign(key, value, envir = cache)
    value
  }
  matrix_score <- function(gamma) {
    gamma <- as.matrix(gamma)
    if (ncol(gamma) != p || any(!gamma %in% c(0L, 1L))) {
      stop("A support score received invalid supports.", call. = FALSE)
    }
    n_supports <- nrow(gamma)
    keys <- apply(gamma, 1L, paste0, collapse = "")
    counters$calls <- counters$calls + n_supports
    out <- numeric(n_supports)
    cached <- vapply(keys, exists, logical(1), envir = cache, inherits = FALSE)
    if (any(cached)) {
      out[cached] <- vapply(keys[cached], get, numeric(1), envir = cache,
                            inherits = FALSE)
      counters$cache_hits <- counters$cache_hits + sum(cached)
    }
    missing_rows <- which(!cached)
    if (length(missing_rows)) {
      unique_missing_key <- unique(keys[missing_rows])
      first <- missing_rows[match(unique_missing_key, keys[missing_rows])]
      supports <- gamma[first, , drop = FALSE]
      increment <- numeric(nrow(supports))
      for (group in seq_along(groups)) {
        indices <- groups[[group]]
        code <- 1L + as.integer(
          supports[, indices, drop = FALSE] %*%
            2^(seq_along(indices) - 1L)
        )
        increment <- increment + unname(local_tables[[group]][code])
      }
      size <- rowSums(supports)
      values <- base_score$log_constant + increment +
        size * log(theta) + (p - size) * log1p(-theta)
      for (index in seq_along(unique_missing_key)) {
        assign(unique_missing_key[[index]], values[[index]], envir = cache)
      }
      out[missing_rows] <- values[match(keys[missing_rows], unique_missing_key)]
      counters$evaluations <- counters$evaluations + length(unique_missing_key)
      counters$cache_hits <- counters$cache_hits +
        length(missing_rows) - length(unique_missing_key)
    }
    out
  }
  attr(score, "matrix_score") <- matrix_score
  diagnostics <- function() {
    c(
      calls = counters$calls,
      evaluations = counters$evaluations,
      cache_hits = counters$cache_hits,
      cached_supports = length(ls(cache, all.names = TRUE))
    )
  }
  reset <- function(clear_cache = TRUE) {
    counters$calls <- 0L
    counters$evaluations <- 0L
    counters$cache_hits <- 0L
    if (clear_cache) {
      rm(list = ls(cache, all.names = TRUE), envir = cache)
    }
    invisible(NULL)
  }
  base_score$score <- score
  base_score$log_likelihood_increment <- log_likelihood_increment
  base_score$diagnostics <- diagnostics
  base_score$reset <- reset
  base_score$factorized_groups <- groups
  base_score$factorized_tables <- local_tables
  class(base_score) <- c(
    "cb_factorized_known_sigma_score",
    class(base_score)
  )
  base_score
}

cb_block_orthogonal_case <- function(p,
                                     n_train = ceiling(1.6 * p),
                                     n_test = 500L,
                                     rho = 0.98,
                                     signal = 1.4,
                                     sigma = 1,
                                     seed = 1L,
                                     heterogeneous = FALSE) {
  if (p %% 2L != 0L || n_train < 3L * p / 2L) {
    stop(
      "The block-orthogonal design needs even p and n_train at least 3p/2.",
      call. = FALSE
    )
  }
  set.seed(seed)
  n_groups <- p / 2L
  raw_basis <- matrix(rnorm(n_train * 3L * n_groups), nrow = n_train)
  basis <- qr.Q(qr(raw_basis))[, seq_len(3L * n_groups), drop = FALSE]
  X_train <- matrix(0, nrow = n_train, ncol = p)
  for (group in seq_len(n_groups)) {
    basis_indices <- (3L * group - 2L):(3L * group)
    predictor_indices <- (2L * group - 1L):(2L * group)
    local_rho <- if (heterogeneous && group %% 3L == 0L) {
      max(0.15, rho - 0.65)
    } else {
      rho
    }
    X_train[, predictor_indices[1L]] <-
      sqrt(local_rho) * basis[, basis_indices[1L]] +
      sqrt(1 - local_rho) * basis[, basis_indices[2L]]
    X_train[, predictor_indices[2L]] <-
      sqrt(local_rho) * basis[, basis_indices[1L]] +
      sqrt(1 - local_rho) * basis[, basis_indices[3L]]
  }
  X_train <- sweep(X_train, 2L, sqrt(colSums(X_train^2)), "/") *
    sqrt(n_train)

  latent_test <- matrix(rnorm(n_test * n_groups), nrow = n_test)
  X_test <- matrix(0, nrow = n_test, ncol = p)
  for (group in seq_len(n_groups)) {
    predictor_indices <- (2L * group - 1L):(2L * group)
    local_rho <- if (heterogeneous && group %% 3L == 0L) {
      max(0.15, rho - 0.65)
    } else {
      rho
    }
    X_test[, predictor_indices] <-
      sqrt(local_rho) * latent_test[, group] +
      sqrt(1 - local_rho) *
        matrix(rnorm(n_test * 2L), nrow = n_test)
  }
  X_test <- scale(X_test)

  beta <- numeric(p)
  active_groups <- unique(pmin(c(1L, 3L, 5L), n_groups))
  beta[2L * active_groups - 1L] <- signal * c(1, -0.8, 0.7)[
    seq_along(active_groups)
  ]
  y_train <- as.numeric(X_train %*% beta + rnorm(n_train, sd = sigma))
  y_test <- as.numeric(X_test %*% beta + rnorm(n_test, sd = sigma))
  y_center <- mean(y_train)
  list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train - y_center,
    y_test = y_test - y_center,
    beta = beta,
    true_support = as.integer(beta != 0),
    groups = split(seq_len(p), ceiling(seq_len(p) / 2)),
    rho = rho,
    sigma = sigma,
    seed = seed,
    heterogeneous = heterogeneous
  )
}

cb_factorized_group_truth <- function(score_object, groups) {
  inclusion_probability <- score_object$inclusion_probability
  group_tables <- lapply(groups, function(indices) {
    local_supports <- as.matrix(expand.grid(
      rep(list(c(0L, 1L)), length(indices))
    ))
    values <- apply(local_supports, 1L, function(local) {
      support <- integer(score_object$p)
      support[indices] <- local
      score_object$score(support) - score_object$log_constant -
        (score_object$p - length(indices)) *
          log1p(-inclusion_probability)
    })
    normalizer <- cb_log_sum_exp(values)
    list(
      supports = local_supports,
      log_weights = values - normalizer,
      weights = exp(values - normalizer),
      log_evidence = normalizer
    )
  })
  pip <- numeric(score_object$p)
  group_state_probabilities <- matrix(
    0,
    nrow = length(groups),
    ncol = 3L,
    dimnames = list(NULL, c("zero", "one", "multi"))
  )
  for (group in seq_along(groups)) {
    table <- group_tables[[group]]
    pip[groups[[group]]] <- cb_pip(table$supports, table$weights)
    counts <- rowSums(table$supports)
    group_state_probabilities[group, ] <- c(
      sum(table$weights[counts == 0L]),
      sum(table$weights[counts == 1L]),
      sum(table$weights[counts >= 2L])
    )
  }
  list(
    pip = pip,
    group_state_probabilities = group_state_probabilities,
    group_tables = group_tables,
    log_evidence = score_object$log_constant +
      sum(vapply(group_tables, `[[`, numeric(1), "log_evidence"))
  )
}

cb_factorized_region_mass <- function(truth,
                                      target_groups,
                                      declared_groups,
                                      states,
                                      fixed_bits = integer(0)) {
  target_groups <- lapply(target_groups, as.integer)
  declared_groups <- lapply(declared_groups, as.integer)
  states <- as.character(states)
  p <- max(unlist(target_groups, use.names = FALSE))
  if (!identical(sort(unlist(target_groups, use.names = FALSE)), seq_len(p)) ||
      length(declared_groups) != length(states) ||
      any(!states %in% c("free", "zero", "one", "multi")) ||
      !setequal(unlist(declared_groups, use.names = FALSE), seq_len(p)) ||
      length(truth$group_tables) != length(target_groups)) {
    stop("Invalid factorized regional-mass inputs.", call. = FALSE)
  }
  fixed_bit_names <- names(fixed_bits)
  fixed_bits <- stats::setNames(as.integer(fixed_bits), fixed_bit_names)
  if (length(fixed_bits) &&
      (is.null(names(fixed_bits)) || anyNA(as.integer(names(fixed_bits))) ||
       any(!fixed_bits %in% c(0L, 1L)) ||
       any(as.integer(names(fixed_bits)) < 1L | as.integer(names(fixed_bits)) > p) ||
       anyDuplicated(as.integer(names(fixed_bits))))) {
    stop("fixed_bits must be a named binary vector indexed by predictor.", call. = FALSE)
  }
  make_factor <- function(vars, values) {
    vars <- as.integer(vars)
    values <- as.numeric(values)
    if (length(values) != 2^length(vars) || any(!is.finite(values)) ||
        any(values < 0)) {
      stop("A binary factor has invalid values.", call. = FALSE)
    }
    list(vars = vars, values = values)
  }
  factors <- lapply(seq_along(target_groups), function(group) {
    table <- truth$group_tables[[group]]
    supports <- as.matrix(table$supports)
    code <- 1L + as.integer(supports %*% 2^(seq_len(ncol(supports)) - 1L))
    values <- numeric(nrow(supports))
    values[code] <- table$weights
    make_factor(target_groups[[group]], values)
  })
  for (group in seq_along(declared_groups)) {
    if (states[[group]] == "free") next
    size <- length(declared_groups[[group]])
    supports <- as.matrix(expand.grid(rep(list(c(0L, 1L)), size)))
    count <- rowSums(supports)
    allowed <- switch(
      states[[group]], zero = count == 0L, one = count == 1L,
      multi = count >= 2L
    )
    factors[[length(factors) + 1L]] <- make_factor(
      declared_groups[[group]], as.numeric(allowed)
    )
  }
  if (length(fixed_bits)) {
    for (index in seq_along(fixed_bits)) {
      variable <- as.integer(names(fixed_bits)[[index]])
      value <- fixed_bits[[index]]
      factors[[length(factors) + 1L]] <- make_factor(
        variable, if (value == 0L) c(1, 0) else c(0, 1)
      )
    }
  }
  evaluate_factor <- function(factor, assignments, union_vars) {
    position <- match(factor$vars, union_vars)
    code <- 1L + as.integer(
      assignments[, position, drop = FALSE] %*%
        2^(seq_along(position) - 1L)
    )
    factor$values[code]
  }
  for (variable in seq_len(p)) {
    selected <- which(vapply(factors, function(factor) {
      variable %in% factor$vars
    }, logical(1)))
    if (!length(selected)) next
    active <- factors[selected]
    union_vars <- sort(unique(unlist(lapply(active, `[[`, "vars"))))
    assignments <- as.matrix(expand.grid(
      rep(list(c(0L, 1L)), length(union_vars)), KEEP.OUT.ATTRS = FALSE
    ))
    product <- Reduce(`*`, lapply(active, evaluate_factor,
      assignments = assignments, union_vars = union_vars
    ))
    remaining_vars <- setdiff(union_vars, variable)
    if (!length(remaining_vars)) {
      replacement <- make_factor(integer(0), sum(product))
    } else {
      remaining_code <- 1L + as.integer(
        assignments[, match(remaining_vars, union_vars), drop = FALSE] %*%
          2^(seq_along(remaining_vars) - 1L)
      )
      values <- rowsum(product, remaining_code, reorder = TRUE)[, 1L]
      replacement <- make_factor(remaining_vars, values)
    }
    factors <- c(factors[-selected], list(replacement))
  }
  mass <- prod(vapply(factors, function(factor) factor$values[[1L]], numeric(1)))
  if (!is.finite(mass) || mass < -1e-12 || mass > 1 + 1e-12) {
    stop("Factorized regional mass is outside [0,1].", call. = FALSE)
  }
  min(1, max(0, mass))
}

cb_factorized_region_log_g_bounds <- function(score_object,
                                              truth,
                                              target_groups,
                                              configuration_tables,
                                              declared_groups,
                                              states) {
  target_groups <- lapply(target_groups, as.integer)
  declared_groups <- lapply(declared_groups, as.integer)
  states <- as.character(states)
  p <- score_object$p
  if (length(truth$group_tables) != length(target_groups) ||
      length(configuration_tables) != length(declared_groups) ||
      length(states) != length(declared_groups) ||
      any(!states %in% c("free", "zero", "one", "multi")) ||
      !setequal(unlist(target_groups, use.names = FALSE), seq_len(p)) ||
      !setequal(unlist(declared_groups, use.names = FALSE), seq_len(p))) {
    stop("Invalid factorized log-g-bound inputs.", call. = FALSE)
  }
  make_factor <- function(vars, values) list(
    vars = as.integer(vars), values = as.numeric(values)
  )
  factors <- lapply(seq_along(target_groups), function(group) {
    table <- truth$group_tables[[group]]
    supports <- as.matrix(table$supports)
    code <- 1L + as.integer(supports %*% 2^(seq_len(ncol(supports)) - 1L))
    values <- rep(NA_real_, nrow(supports))
    values[code] <- table$log_weights + table$log_evidence
    make_factor(target_groups[[group]], values)
  })
  for (group in seq_along(declared_groups)) {
    table <- configuration_tables[[group]]
    supports <- as.matrix(table$supports)
    count <- rowSums(supports)
    inside <- switch(
      states[[group]], free = rep(TRUE, nrow(supports)), zero = count == 0L,
      one = count == 1L, multi = count >= 2L
    )
    probability <- as.numeric(table$probabilities)
    conditional <- rep(NA_real_, length(probability))
    conditional[inside] <- probability[inside] / sum(probability[inside])
    code <- 1L + as.integer(supports %*% 2^(seq_len(ncol(supports)) - 1L))
    values <- rep(NA_real_, nrow(supports))
    values[code[inside]] <- -log(conditional[inside])
    factors[[length(factors) + 1L]] <- make_factor(
      declared_groups[[group]], values
    )
  }
  factor_values <- function(factor, assignments, union_vars) {
    position <- match(factor$vars, union_vars)
    code <- 1L + as.integer(
      assignments[, position, drop = FALSE] %*%
        2^(seq_along(position) - 1L)
    )
    factor$values[code]
  }
  eliminate <- function(operation) {
    active_factors <- factors
    for (variable in seq_len(p)) {
      selected <- which(vapply(active_factors, function(factor) {
        variable %in% factor$vars
      }, logical(1)))
      if (!length(selected)) next
      active <- active_factors[selected]
      union_vars <- sort(unique(unlist(lapply(active, `[[`, "vars"))))
      assignments <- as.matrix(expand.grid(
        rep(list(c(0L, 1L)), length(union_vars)), KEEP.OUT.ATTRS = FALSE
      ))
      values <- Reduce(`+`, lapply(active, factor_values,
        assignments = assignments, union_vars = union_vars
      ))
      remaining <- setdiff(union_vars, variable)
      if (!length(remaining)) {
        finite <- values[is.finite(values)]
        replacement <- make_factor(integer(0), operation(finite))
      } else {
        code <- 1L + as.integer(
          assignments[, match(remaining, union_vars), drop = FALSE] %*%
            2^(seq_along(remaining) - 1L)
        )
        split_values <- split(values, code)
        reduced <- vapply(split_values, function(value) {
          finite <- value[is.finite(value)]
          if (length(finite)) operation(finite) else NA_real_
        }, numeric(1))
        replacement <- make_factor(remaining, reduced)
      }
      active_factors <- c(active_factors[-selected], list(replacement))
    }
    constants <- vapply(active_factors, function(factor) factor$values[[1L]], numeric(1))
    if (any(!is.finite(constants))) {
      stop("A declared region has no finite configuration.", call. = FALSE)
    }
    score_object$log_constant + sum(constants)
  }
  bounds <- c(eliminate(min), eliminate(max))
  if (any(!is.finite(bounds)) || bounds[[1L]] > bounds[[2L]]) {
    stop("Factorized log-g bounds are invalid.", call. = FALSE)
  }
  bounds
}

cb_collapse_particles <- function(particles, weights) {
  particles <- as.matrix(particles)
  weights <- as.numeric(weights)
  if (nrow(particles) != length(weights) ||
      any(!is.finite(weights)) ||
      sum(weights) <= 0) {
    stop("Invalid particles or weights.", call. = FALSE)
  }
  keys <- apply(particles, 1L, paste0, collapse = "")
  unique_keys <- unique(keys)
  collapsed_weights <- vapply(
    unique_keys,
    function(key) sum(weights[keys == key]),
    numeric(1)
  )
  first <- match(unique_keys, keys)
  list(
    supports = particles[first, , drop = FALSE],
    weights = collapsed_weights / sum(collapsed_weights),
    keys = unique_keys
  )
}

cb_combine_region_fits <- function(fits, log_evidence = NULL) {
  if (length(fits) < 1L) {
    stop("At least one regional fit is required.", call. = FALSE)
  }
  if (is.null(log_evidence)) {
    log_evidence <- vapply(fits, `[[`, numeric(1), "log_evidence")
  }
  region_weights <- exp(log_evidence - cb_log_sum_exp(log_evidence))
  particles <- do.call(rbind, lapply(fits, `[[`, "particles"))
  weights <- unlist(Map(
    function(fit, region_weight) fit$weights * region_weight,
    fits,
    region_weights
  ))
  collapsed <- cb_collapse_particles(particles, weights)
  collapsed$region_weights <- region_weights
  collapsed$log_evidence <- cb_log_sum_exp(log_evidence)
  collapsed
}

cb_pip <- function(supports, weights) {
  supports <- as.matrix(supports)
  weights <- as.numeric(weights) / sum(weights)
  as.numeric(crossprod(weights, supports))
}

cb_group_state_keys <- function(supports, groups) {
  supports <- as.matrix(supports)
  state <- vapply(groups, function(indices) {
    count <- rowSums(supports[, indices, drop = FALSE])
    ifelse(count == 0L, "zero", ifelse(count == 1L, "one", "multi"))
  }, character(nrow(supports)))
  if (is.vector(state)) {
    state <- matrix(state, ncol = length(groups))
  }
  apply(state, 1L, paste, collapse = "/")
}

cb_select_pilot_cells <- function(particles,
                                  weights,
                                  groups,
                                  mass_particles = NULL,
                                  mass_weights = NULL,
                                  n_split_groups = 5L,
                                  n_cells = 8L,
                                  target_mass = 0.9,
                                  minimum_confidence = 0,
                                  split_rule = c(
                                    "confidence",
                                    "uncertainty"
                                  )) {
  split_rule <- match.arg(split_rule)
  particles <- as.matrix(particles)
  weights <- as.numeric(weights) / sum(weights)
  if (is.null(mass_particles)) {
    mass_particles <- particles
  }
  if (is.null(mass_weights)) {
    mass_weights <- weights
  }
  mass_particles <- as.matrix(mass_particles)
  mass_weights <- as.numeric(mass_weights)
  if (ncol(mass_particles) != ncol(particles) ||
      length(mass_weights) != nrow(mass_particles) ||
      any(!is.finite(mass_weights)) ||
      any(mass_weights < 0) ||
      sum(mass_weights) <= 0) {
    stop(
      "The independent cell-mass sample is incompatible with the pilot sample.",
      call. = FALSE
    )
  }
  mass_weights <- mass_weights / sum(mass_weights)
  group_probabilities <- vapply(groups, function(indices) {
    sum(weights * (rowSums(particles[, indices, drop = FALSE]) > 0))
  }, numeric(1))
  state_probabilities <- t(vapply(groups, function(indices) {
    count <- rowSums(particles[, indices, drop = FALSE])
    c(
      zero = sum(weights[count == 0L]),
      one = sum(weights[count == 1L]),
      multi = sum(weights[count >= 2L])
    )
  }, numeric(3)))
  split_score <- if (split_rule == "confidence") {
    apply(state_probabilities, 1L, max)
  } else {
    1 - apply(state_probabilities, 1L, max)
  }
  if (!is.finite(target_mass) ||
      target_mass <= 0 ||
      target_mass > 1 ||
      !is.finite(minimum_confidence) ||
      minimum_confidence < 0 ||
      minimum_confidence > 1) {
    stop(
      "target_mass and minimum_confidence must lie in (0,1] and [0,1].",
      call. = FALSE
    )
  }
  order_groups <- order(split_score, decreasing = TRUE)
  eligible <- if (split_rule == "confidence") {
    order_groups[split_score[order_groups] >= minimum_confidence]
  } else {
    order_groups
  }
  split_groups <- eligible[
    seq_len(min(n_split_groups, length(eligible)))
  ]
  if (length(split_groups) == 0L) {
    return(list(
      split_groups = integer(0),
      groups = list(),
      states = list(character(0)),
      pilot_mass = 1,
      selected_pilot_mass = 1,
      group_probabilities = group_probabilities,
      state_probabilities = state_probabilities,
      split_rule = split_rule,
      target_mass = target_mass,
      minimum_confidence = minimum_confidence
    ))
  }
  restricted_groups <- groups[split_groups]
  keys <- cb_group_state_keys(mass_particles, restricted_groups)
  cell_mass <- sort(
    vapply(
      unique(keys),
      function(key) sum(mass_weights[keys == key]),
      numeric(1)
    ),
    decreasing = TRUE
  )
  required_cells <- which(cumsum(cell_mass) >= target_mass)[1L]
  if (is.na(required_cells)) {
    required_cells <- length(cell_mass)
  }
  selected_count <- min(n_cells, required_cells, length(cell_mass))
  selected_keys <- names(cell_mass)[seq_len(selected_count)]
  states <- strsplit(selected_keys, "/", fixed = TRUE)
  list(
    split_groups = split_groups,
    groups = restricted_groups,
    states = states,
    pilot_mass = unname(cell_mass[selected_keys]),
    selected_pilot_mass = sum(cell_mass[selected_keys]),
    group_probabilities = group_probabilities,
    state_probabilities = state_probabilities,
    split_rule = split_rule,
    target_mass = target_mass,
    minimum_confidence = minimum_confidence
  )
}

cb_make_cell_regions <- function(inclusion_probabilities,
                                 all_groups,
                                 split_groups,
                                 states,
                                 new_group_state_region_fn) {
  lapply(seq_along(states), function(index) {
    full_states <- rep("free", length(all_groups))
    full_states[split_groups] <- states[[index]]
    new_group_state_region_fn(
      inclusion_probabilities = inclusion_probabilities,
      groups = all_groups,
      states = full_states,
      cost = 1L + sum(full_states != "free"),
      label = paste0("pilot-cell-", index)
    )
  })
}
