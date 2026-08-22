rrs_log_sum_exp <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(-Inf)
  }
  m <- max(x[finite])
  m + log(sum(exp(x[finite] - m)))
}

rrs_normalize_log_weights <- function(log_weights) {
  normalizer <- rrs_log_sum_exp(log_weights)
  if (!is.finite(normalizer)) {
    stop("Particle weights have no finite normalizer.", call. = FALSE)
  }
  log_weights - normalizer
}

rrs_effective_sample_size <- function(log_weights) {
  weights <- exp(rrs_normalize_log_weights(log_weights))
  1 / sum(weights^2)
}

rrs_registered_status <- function(status) {
  mapping <- c(
    certified = "bounds-met",
    frontier_exhausted = "tolerance-not-met",
    evidence_interval_conflict = "interval-conflict",
    uninformative_evidence = "uninformative-evidence",
    budget_exhausted = "budget-exhausted",
    search_incomplete = "implementation-error"
  )
  status <- as.character(status)
  out <- unname(mapping[status])
  out[is.na(out)] <- "implementation-error"
  out
}

rrs_common_accounting <- function(score_diagnostics,
                                  search = NULL,
                                  downstream_score_requests = 0,
                                  elapsed_seconds,
                                  peak_memory_bytes = NA_real_) {
  diagnostic_names <- names(score_diagnostics)
  score_diagnostics <- as.numeric(score_diagnostics)
  names(score_diagnostics) <- diagnostic_names
  required <- c("calls", "evaluations")
  if (is.null(names(score_diagnostics)) ||
      any(!required %in% names(score_diagnostics)) ||
      length(elapsed_seconds) != 1L || !is.finite(elapsed_seconds) ||
      elapsed_seconds < 0 ||
      length(downstream_score_requests) != 1L ||
      !is.finite(downstream_score_requests) ||
      downstream_score_requests < 0 ||
      length(peak_memory_bytes) != 1L ||
      (!is.na(peak_memory_bytes) &&
       (!is.finite(peak_memory_bytes) || peak_memory_bytes < 0))) {
    stop("Invalid common computation-accounting inputs.", call. = FALSE)
  }
  raw_requests <- unname(score_diagnostics[["calls"]])
  distinct_evaluations <- unname(score_diagnostics[["evaluations"]])
  if (raw_requests < distinct_evaluations) {
    stop("Distinct score evaluations cannot exceed raw requests.", call. = FALSE)
  }
  bound_requests <- if (is.null(search)) 0 else search$bound_score_evaluations
  evidence_requests <- if (is.null(search)) 0 else search$evidence_score_evaluations
  selection_requests <- bound_requests + evidence_requests
  if (!is.null(search) &&
      selection_requests != search$total_score_evaluations) {
    stop("Selection-stage score accounting is internally inconsistent.", call. = FALSE)
  }
  data.frame(
    raw_score_requests = raw_requests,
    distinct_score_evaluations = distinct_evaluations,
    bound_score_requests = bound_requests,
    evidence_score_requests = evidence_requests,
    downstream_score_requests = downstream_score_requests,
    selection_score_requests = selection_requests,
    elapsed_seconds = elapsed_seconds,
    peak_memory_bytes = peak_memory_bytes,
    stringsAsFactors = FALSE
  )
}

rrs_conditional_ess <- function(log_weights, log_increment) {
  log_weights <- rrs_normalize_log_weights(log_weights)
  n_particles <- length(log_weights)
  log_first <- rrs_log_sum_exp(log_weights + log_increment)
  log_second <- rrs_log_sum_exp(log_weights + 2 * log_increment)
  exp(log(n_particles) + 2 * log_first - log_second)
}

rrs_systematic_resample <- function(weights) {
  n_particles <- length(weights)
  cumulative <- cumsum(weights / sum(weights))
  positions <- (runif(1) + 0:(n_particles - 1L)) / n_particles
  pmin(findInterval(positions, cumulative) + 1L, n_particles)
}

rrs_support_keys <- function(supports) {
  supports <- as.matrix(supports)
  apply(supports, 1L, paste, collapse = ",")
}

new_enumerated_region <- function(supports,
                                  proposal_prob = NULL,
                                  cost = 1,
                                  label = NULL) {
  supports <- as.matrix(supports)
  if (nrow(supports) < 1L) {
    stop("An enumerated region must contain at least one support.", call. = FALSE)
  }
  if (is.null(proposal_prob)) {
    proposal_prob <- rep(1 / nrow(supports), nrow(supports))
  }
  if (length(proposal_prob) != nrow(supports) ||
      any(!is.finite(proposal_prob)) ||
      any(proposal_prob <= 0)) {
    stop("proposal_prob must be positive and match the region size.", call. = FALSE)
  }
  proposal_prob <- proposal_prob / sum(proposal_prob)
  keys <- rrs_support_keys(supports)
  if (anyDuplicated(keys)) {
    stop("The supports in an enumerated region must be distinct.", call. = FALSE)
  }

  sample_region <- function(n) {
    supports[
      sample.int(nrow(supports), size = n, replace = TRUE, prob = proposal_prob),
      ,
      drop = FALSE
    ]
  }
  log_prob <- function(gamma) {
    gamma <- as.matrix(gamma)
    index <- match(rrs_support_keys(gamma), keys)
    out <- rep(-Inf, nrow(gamma))
    inside <- !is.na(index)
    out[inside] <- log(proposal_prob[index[inside]])
    out
  }

  structure(
    list(
      sample = sample_region,
      log_prob = log_prob,
      cost = cost,
      label = if (is.null(label)) "region" else as.character(label),
      supports = supports,
      proposal_prob = proposal_prob
    ),
    class = "rrs_region"
  )
}

new_group_table_region <- function(configuration_tables,
                                   states,
                                   groups = NULL,
                                   cost = 1L,
                                   label = NULL) {
  states <- as.character(states)
  if (length(configuration_tables) < 1L ||
      length(states) != length(configuration_tables) ||
      any(!states %in% c("free", "zero", "one", "multi"))) {
    stop("Invalid group-table region inputs.", call. = FALSE)
  }
  sizes <- vapply(configuration_tables, function(table) {
    ncol(as.matrix(table$supports))
  }, integer(1))
  if (is.null(groups)) {
    groups <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
  }
  if (length(groups) != length(sizes) ||
      !identical(vapply(groups, length, integer(1)), sizes) ||
      !setequal(unlist(groups, use.names = FALSE), seq_len(sum(sizes)))) {
    stop("The screened configuration tables do not match their groups.", call. = FALSE)
  }
  p <- sum(sizes)
  tables <- lapply(seq_along(configuration_tables), function(group) {
    table <- configuration_tables[[group]]
    supports <- as.matrix(table$supports)
    probabilities <- as.numeric(table$probabilities)
    if (nrow(supports) != 2^ncol(supports) ||
        length(probabilities) != nrow(supports) ||
        any(!supports %in% c(0L, 1L)) ||
        any(!is.finite(probabilities)) || any(probabilities <= 0)) {
      stop("Every screened group table must have positive full support.", call. = FALSE)
    }
    counts <- rowSums(supports)
    allowed <- switch(
      states[group],
      free = seq_len(nrow(supports)),
      zero = which(counts == 0L),
      one = which(counts == 1L),
      multi = which(counts >= 2L)
    )
    if (length(allowed) == 0L) {
      stop("A declared group state has no configurations.", call. = FALSE)
    }
    conditional <- numeric(nrow(supports))
    conditional[allowed] <- probabilities[allowed] / sum(probabilities[allowed])
    list(
      supports = supports,
      conditional = conditional,
      allowed = allowed,
      anchor = allowed[which.max(conditional[allowed])]
    )
  })
  sample_region <- function(n) {
    n <- as.integer(n)
    out <- matrix(0L, nrow = n, ncol = p)
    for (group in seq_along(groups)) {
      table <- tables[[group]]
      selected <- table$allowed[sample.int(
        length(table$allowed),
        size = n,
        replace = TRUE,
        prob = table$conditional[table$allowed]
      )]
      out[, groups[[group]]] <- table$supports[selected, , drop = FALSE]
    }
    out
  }
  log_prob_inside <- function(gamma) {
    gamma <- as.matrix(gamma)
    if (ncol(gamma) != p || any(!gamma %in% c(0L, 1L))) {
      stop("A group-table proposal received invalid supports.", call. = FALSE)
    }
    out <- numeric(nrow(gamma))
    valid <- rep(TRUE, nrow(gamma))
    for (group in seq_along(groups)) {
      local <- gamma[, groups[[group]], drop = FALSE]
      code <- 1L + as.integer(local %*% 2^(seq_len(ncol(local)) - 1L))
      probabilities <- tables[[group]]$conditional[code]
      valid <- valid & probabilities > 0
      out <- out + ifelse(probabilities > 0, log(probabilities), 0)
    }
    out[!valid] <- -Inf
    out
  }
  anchor <- function() {
    out <- matrix(0L, nrow = 1L, ncol = p)
    for (group in seq_along(groups)) {
      table <- tables[[group]]
      out[1L, groups[[group]]] <- table$supports[table$anchor, ]
    }
    out
  }
  structure(
    list(
      sample = sample_region,
      log_prob = log_prob_inside,
      log_prob_inside = log_prob_inside,
      anchor = anchor,
      cost = as.integer(cost),
      label = if (is.null(label)) paste(states, collapse = "/") else label,
      states = states,
      groups = groups,
      inclusion_probabilities = rep(0.5, p),
      configuration_tables = configuration_tables
    ),
    class = "rrs_region"
  )
}

make_screened_group_region_factory <- function(log_score_fn, p, groups) {
  if (!is.function(log_score_fn) || length(p) != 1L || p < 1L ||
      length(groups) < 1L ||
      !setequal(unlist(groups, use.names = FALSE), seq_len(p))) {
    stop("Invalid screened group-proposal inputs.", call. = FALSE)
  }
  configuration_tables <- lapply(groups, function(indices) {
    supports <- as.matrix(expand.grid(rep(list(c(0L, 1L)), length(indices))))
    full_supports <- matrix(0L, nrow = nrow(supports), ncol = p)
    full_supports[, indices] <- supports
    log_values <- rrs_evaluate_log_score(log_score_fn, full_supports)
    list(
      supports = supports,
      probabilities = exp(log_values - rrs_log_sum_exp(log_values))
    )
  })
  list(
    configuration_tables = configuration_tables,
    region_factory = function(states, cost, label) {
      region <- new_group_table_region(
        configuration_tables = configuration_tables,
        states = states,
        groups = groups,
        cost = cost,
        label = label
      )
      region
    }
  )
}

rrs_validate_region <- function(region) {
  required <- c("sample", "log_prob", "cost", "label")
  missing <- setdiff(required, names(region))
  if (length(missing) > 0L) {
    stop(
      sprintf("Region is missing fields: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.function(region$sample) || !is.function(region$log_prob)) {
    stop("Region sampling and log-probability fields must be functions.", call. = FALSE)
  }
  if (length(region$cost) != 1L || !is.finite(region$cost) || region$cost <= 0) {
    stop("Region cost must be one positive finite number.", call. = FALSE)
  }
  invisible(TRUE)
}

rrs_evaluate_log_score <- function(log_score_fn, particles) {
  particles <- as.matrix(particles)
  matrix_score <- attr(log_score_fn, "matrix_score", exact = TRUE)
  out <- if (is.function(matrix_score)) {
    as.numeric(matrix_score(particles))
  } else {
    vapply(
      seq_len(nrow(particles)),
      function(i) log_score_fn(particles[i, ]),
      numeric(1)
    )
  }
  if (length(out) != nrow(particles)) {
    stop("A vectorized log-score returned the wrong length.", call. = FALSE)
  }
  if (any(!is.finite(out))) {
    stop("log_score_fn returned a non-finite value.", call. = FALSE)
  }
  out
}

rrs_budgeted_log_score <- function(log_score_fn, budget) {
  if (!is.function(log_score_fn) || length(budget) != 1L ||
      is.na(budget) || budget < 0) {
    stop("Invalid score function or score budget.", call. = FALSE)
  }
  used <- 0
  reserve <- function(requested) {
    requested <- as.numeric(requested)
    if (length(requested) != 1L || !is.finite(requested) || requested < 0) {
      stop("A score batch declared an invalid request size.", call. = FALSE)
    }
    if (used + requested > budget) {
      condition <- structure(
        list(
          message = "The next score batch would exceed the declared budget.",
          call = NULL,
          used = used,
          requested = requested,
          budget = budget
        ),
        class = c("rrs_score_budget_exhausted", "error", "condition")
      )
      stop(condition)
    }
    used <<- used + requested
  }
  wrapped <- function(gamma) {
    reserve(1)
    log_score_fn(gamma)
  }
  matrix_score <- attr(log_score_fn, "matrix_score", exact = TRUE)
  if (is.function(matrix_score)) {
    attr(wrapped, "matrix_score") <- function(particles) {
      particles <- as.matrix(particles)
      reserve(nrow(particles))
      matrix_score(particles)
    }
  }
  attr(wrapped, "score_budget_usage") <- function() used
  wrapped
}

rrs_choose_temperature <- function(phi,
                                   log_weights,
                                   log_g,
                                   ess_fraction,
                                   tolerance = 1e-8) {
  n_particles <- length(log_weights)
  target <- ess_fraction * n_particles
  final_cess <- rrs_conditional_ess(
    log_weights,
    (1 - phi) * log_g
  )
  if (final_cess >= target) {
    return(1)
  }

  lower <- phi
  upper <- 1
  for (iteration in seq_len(80L)) {
    middle <- (lower + upper) / 2
    middle_cess <- rrs_conditional_ess(
      log_weights,
      (middle - phi) * log_g
    )
    if (middle_cess >= target) {
      lower <- middle
    } else {
      upper <- middle
    }
    if (upper - lower <= tolerance) {
      break
    }
  }
  max(lower, phi + tolerance)
}

rrs_region_mutation_proposal <- function(region, particles) {
  particles <- as.matrix(particles)
  if (is.function(region$mutate)) {
    proposal <- region$mutate(particles)
    if (is.matrix(proposal)) {
      proposal <- list(
        particles = proposal,
        proposed = rep(TRUE, nrow(particles)),
        type = "region-local"
      )
    }
    required <- c("particles", "proposed")
    if (!is.list(proposal) || any(!required %in% names(proposal))) {
      stop(
        "A regional mutation must return particles and proposed indicators.",
        call. = FALSE
      )
    }
    proposal$particles <- as.matrix(proposal$particles)
    proposal$proposed <- as.logical(proposal$proposed)
    if (!identical(dim(proposal$particles), dim(particles)) ||
        length(proposal$proposed) != nrow(particles)) {
      stop("A regional mutation returned incompatible dimensions.", call. = FALSE)
    }
    if (is.null(proposal$type)) {
      proposal$type <- "region-local"
    }
    return(proposal)
  }
  list(
    particles = as.matrix(region$sample(nrow(particles))),
    proposed = rep(TRUE, nrow(particles)),
    type = "independence"
  )
}

rrs_reversible_mutation <- function(particles,
                                    log_g,
                                    phi,
                                    region,
                                    log_score_fn,
                                    n_steps) {
  if (n_steps <= 0L) {
    return(list(
      particles = particles,
      log_g = log_g,
      accepted = 0L,
      proposed = 0L,
      score_evaluations = 0L,
      type = if (is.function(region$mutate)) {
        "region-local"
      } else {
        "independence"
      }
    ))
  }

  n_particles <- nrow(particles)
  if (is.function(region$target_mutate)) {
    accepted <- 0L
    proposed_count <- 0L
    score_evaluations <- 0L
    mutation_type <- "target-invariant"
    for (step in seq_len(n_steps)) {
      mutation <- region$target_mutate(
        particles = particles,
        phi = phi,
        log_score_fn = log_score_fn
      )
      required <- c("particles", "log_g", "score_evaluations")
      if (!is.list(mutation) || any(!required %in% names(mutation))) {
        stop(
          "A target-invariant mutation returned an incomplete result.",
          call. = FALSE
        )
      }
      particles <- as.matrix(mutation$particles)
      log_g <- as.numeric(mutation$log_g)
      if (nrow(particles) != n_particles ||
          length(log_g) != n_particles ||
          any(!is.finite(log_g))) {
        stop(
          "A target-invariant mutation returned invalid particles.",
          call. = FALSE
        )
      }
      accepted <- accepted + if (is.null(mutation$accepted)) {
        n_particles
      } else {
        mutation$accepted
      }
      proposed_count <- proposed_count + if (is.null(mutation$proposed)) {
        n_particles
      } else {
        mutation$proposed
      }
      score_evaluations <- score_evaluations + mutation$score_evaluations
      if (!is.null(mutation$type)) {
        mutation_type <- mutation$type
      }
    }
    return(list(
      particles = particles,
      log_g = log_g,
      accepted = accepted,
      proposed = proposed_count,
      score_evaluations = score_evaluations,
      type = mutation_type
    ))
  }

  accepted <- 0L
  proposed_count <- 0L
  score_evaluations <- 0L
  mutation_type <- if (is.function(region$mutate)) {
    "region-local"
  } else {
    "independence"
  }
  for (step in seq_len(n_steps)) {
    proposal <- rrs_region_mutation_proposal(region, particles)
    proposed <- proposal$particles
    proposed_rows <- which(proposal$proposed)
    mutation_type <- proposal$type
    if (length(proposed_rows) == 0L) {
      next
    }
    proposed_log_nu <- region$log_prob(proposed[proposed_rows, , drop = FALSE])
    if (any(!is.finite(proposed_log_nu))) {
      stop("A regional mutation left the proposal support.", call. = FALSE)
    }
    proposed_log_score <- rrs_evaluate_log_score(
      log_score_fn,
      proposed[proposed_rows, , drop = FALSE]
    )
    proposed_log_g <- proposed_log_score - proposed_log_nu
    accept_rows <- proposed_rows[
      log(runif(length(proposed_rows))) <=
        phi * (proposed_log_g - log_g[proposed_rows])
    ]
    accept <- rep(FALSE, n_particles)
    accept[accept_rows] <- TRUE
    if (any(accept)) {
      particles[accept, ] <- proposed[accept, , drop = FALSE]
      accepted_positions <- match(accept_rows, proposed_rows)
      log_g[accept_rows] <- proposed_log_g[accepted_positions]
    }
    accepted <- accepted + sum(accept)
    proposed_count <- proposed_count + length(proposed_rows)
    score_evaluations <- score_evaluations + length(proposed_rows)
  }

  list(
    particles = particles,
    log_g = log_g,
    accepted = accepted,
    proposed = proposed_count,
    score_evaluations = score_evaluations,
    type = mutation_type
  )
}

rrs_independence_mutation <- rrs_reversible_mutation

region_annealed_smc <- function(region,
                                log_score_fn,
                                n_particles = 1000L,
                                ess_fraction = 0.8,
                                resample_fraction = 0.5,
                                mutation_steps = 2L,
                                phi_schedule = NULL,
                                temperature_tolerance = 1e-8,
                                max_stages = 1000L) {
  rrs_validate_region(region)
  n_particles <- as.integer(n_particles)
  mutation_steps <- as.integer(mutation_steps)
  if (n_particles < 2L) {
    stop("n_particles must be at least two.", call. = FALSE)
  }
  if (!(ess_fraction > 0 && ess_fraction < 1)) {
    stop("ess_fraction must lie strictly between zero and one.", call. = FALSE)
  }
  if (!(resample_fraction > 0 && resample_fraction <= 1)) {
    stop("resample_fraction must lie in (0, 1].", call. = FALSE)
  }
  if (mutation_steps < 0L) {
    stop("mutation_steps must be nonnegative.", call. = FALSE)
  }

  if (!is.null(phi_schedule)) {
    phi_schedule <- as.numeric(phi_schedule)
    if (length(phi_schedule) < 2L ||
        phi_schedule[1] != 0 ||
        tail(phi_schedule, 1) != 1 ||
        any(diff(phi_schedule) <= 0)) {
      stop("phi_schedule must increase strictly from zero to one.", call. = FALSE)
    }
  }

  particles <- as.matrix(region$sample(n_particles))
  if (nrow(particles) != n_particles) {
    stop("Region sampler returned the wrong number of particles.", call. = FALSE)
  }
  log_nu <- region$log_prob(particles)
  if (any(!is.finite(log_nu))) {
    stop("Region sampler generated a point outside its proposal support.", call. = FALSE)
  }
  log_score <- rrs_evaluate_log_score(log_score_fn, particles)
  log_g <- log_score - log_nu
  log_weights <- rep(-log(n_particles), n_particles)
  log_evidence <- 0
  phi <- 0
  stage <- 0L
  schedule_index <- 1L
  score_evaluations <- n_particles
  diagnostic_rows <- list()

  while (phi < 1 - temperature_tolerance) {
    stage <- stage + 1L
    if (stage > max_stages) {
      stop("The tempering schedule exceeded max_stages.", call. = FALSE)
    }
    if (is.null(phi_schedule)) {
      next_phi <- rrs_choose_temperature(
        phi = phi,
        log_weights = log_weights,
        log_g = log_g,
        ess_fraction = ess_fraction,
        tolerance = temperature_tolerance
      )
    } else {
      schedule_index <- schedule_index + 1L
      next_phi <- phi_schedule[schedule_index]
    }
    next_phi <- min(1, next_phi)
    delta <- next_phi - phi
    if (delta <= 0) {
      stop("The adaptive temperature update failed to advance.", call. = FALSE)
    }

    log_increment <- delta * log_g
    log_ratio <- rrs_log_sum_exp(log_weights + log_increment)
    log_evidence <- log_evidence + log_ratio
    log_weights <- rrs_normalize_log_weights(
      log_weights + log_increment
    )
    ess <- rrs_effective_sample_size(log_weights)

    resampled <- ess <= resample_fraction * n_particles
    if (resampled) {
      ancestors <- rrs_systematic_resample(exp(log_weights))
      particles <- particles[ancestors, , drop = FALSE]
      log_g <- log_g[ancestors]
      log_weights <- rep(-log(n_particles), n_particles)
    }

    mutation <- rrs_reversible_mutation(
      particles = particles,
      log_g = log_g,
      phi = next_phi,
      region = region,
      log_score_fn = log_score_fn,
      n_steps = mutation_steps
    )
    particles <- mutation$particles
    log_g <- mutation$log_g
    score_evaluations <- score_evaluations + mutation$score_evaluations
    acceptance_rate <- if (mutation$proposed > 0L) {
      mutation$accepted / mutation$proposed
    } else {
      NA_real_
    }

    diagnostic_rows[[stage]] <- data.frame(
      stage = stage,
      phi_previous = phi,
      phi = next_phi,
      delta = delta,
      ess = ess,
      resampled = resampled,
      mutation_type = mutation$type,
      mutation_acceptance = acceptance_rate,
      log_evidence = log_evidence
    )
    phi <- next_phi
  }

  diagnostics <- do.call(rbind, diagnostic_rows)
  weights <- exp(rrs_normalize_log_weights(log_weights))
  list(
    label = region$label,
    log_evidence = log_evidence,
    evidence = exp(log_evidence),
    particles = particles,
    weights = weights,
    diagnostics = diagnostics,
    score_evaluations = score_evaluations,
    cost = region$cost,
    final_ess = 1 / sum(weights^2)
  )
}

regional_particle_requirement <- function(rho, relative_error, n_regions, delta) {
  if (rho < 0 ||
      !(relative_error > 0 && relative_error < 1) ||
      n_regions < 1 ||
      !(delta > 0 && delta < 1)) {
    stop("Invalid concentration-bound inputs.", call. = FALSE)
  }
  if (rho == 0) {
    return(1L)
  }
  ceiling(
    ((exp(rho) - 1)^2 / (2 * relative_error^2)) *
      log(2 * n_regions / delta)
  )
}

# Budget preflight for one node with its own predictable failure allocation.
# This is the direct inversion of the adaptive evidence radius used in the
# manuscript.  `regional_particle_requirement()` above is retained for the
# older equal-allocation interface used by the released penalized experiments.
regional_evidence_preflight <- function(rho,
                                        relative_error,
                                        delta_node,
                                        remaining_paths = Inf) {
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 ||
      length(relative_error) != 1L || !is.finite(relative_error) ||
      !(relative_error > 0 && relative_error < 1) ||
      length(delta_node) != 1L || !is.finite(delta_node) ||
      !(delta_node > 0 && delta_node < 1) ||
      length(remaining_paths) != 1L || is.na(remaining_paths) ||
      remaining_paths < 0) {
    stop("Invalid regional evidence preflight inputs.", call. = FALSE)
  }
  required_paths <- if (rho == 0) {
    1L
  } else {
    raw_required_paths <- ceiling(
      ((exp(rho) - 1)^2 / (2 * relative_error^2)) *
        log(2 / delta_node)
    )
    if (!is.finite(raw_required_paths) ||
        raw_required_paths > .Machine$integer.max) {
      Inf
    } else {
      as.integer(raw_required_paths)
    }
  }
  list(
    required_paths = required_paths,
    remaining_paths = remaining_paths,
    feasible = is.infinite(remaining_paths) ||
      required_paths <= remaining_paths,
    relative_error = relative_error,
    delta_node = delta_node,
    rho = rho
  )
}

region_evidence_frontier <- function(log_evidence, costs) {
  log_evidence <- as.numeric(log_evidence)
  raw_costs <- as.numeric(costs)
  if (length(log_evidence) != length(raw_costs) ||
      length(raw_costs) < 1L ||
      any(!is.finite(log_evidence)) ||
      any(!is.finite(raw_costs)) ||
      any(raw_costs < 1) ||
      any(abs(raw_costs - round(raw_costs)) > sqrt(.Machine$double.eps))) {
    stop(
      "log_evidence and positive integer costs must have the same nonzero length.",
      call. = FALSE
    )
  }
  costs <- as.integer(round(raw_costs))

  offset <- max(log_evidence)
  evidence <- exp(log_evidence - offset)
  if (length(unique(costs)) == 1L) {
    order_by_evidence <- order(evidence, decreasing = TRUE)
    cumulative <- cumsum(evidence[order_by_evidence])
    common_cost <- costs[1L]
    return(data.frame(
      budget = seq_along(cumulative) * common_cost,
      log_evidence = log(cumulative) + offset,
      selected = I(lapply(
        seq_along(cumulative),
        function(k) order_by_evidence[seq_len(k)]
      ))
    ))
  }
  max_budget <- sum(costs)
  best <- rep(-Inf, max_budget + 1L)
  selected <- vector("list", max_budget + 1L)
  best[1] <- 0
  selected[[1]] <- integer(0)

  for (j in seq_along(costs)) {
    for (budget in seq.int(max_budget, costs[j], by = -1L)) {
      previous <- budget - costs[j]
      if (!is.finite(best[previous + 1L])) {
        next
      }
      candidate <- best[previous + 1L] + evidence[j]
      if (candidate > best[budget + 1L]) {
        best[budget + 1L] <- candidate
        selected[[budget + 1L]] <- c(selected[[previous + 1L]], j)
      }
    }
  }

  feasible <- which(is.finite(best) & best > 0) - 1L
  data.frame(
    budget = feasible,
    log_evidence = log(best[feasible + 1L]) + offset,
    selected = I(selected[feasible + 1L])
  )
}

rrs_validate_interval_frontier_inputs <- function(log_evidence_lower,
                                                   log_evidence_upper,
                                                   costs,
                                                   epsilon,
                                                   eligible) {
  log_evidence_lower <- as.numeric(log_evidence_lower)
  log_evidence_upper <- as.numeric(log_evidence_upper)
  raw_costs <- as.numeric(costs)
  n_regions <- length(log_evidence_lower)
  inversion <- log_evidence_lower - log_evidence_upper
  rounding_inversion <- is.finite(inversion) & inversion > 0 & inversion <= 1e-10
  if (any(rounding_inversion)) {
    midpoint <- (
      log_evidence_lower[rounding_inversion] +
        log_evidence_upper[rounding_inversion]
    ) / 2
    log_evidence_lower[rounding_inversion] <- midpoint
    log_evidence_upper[rounding_inversion] <- midpoint
  }
  if (n_regions < 1L ||
      length(log_evidence_upper) != n_regions ||
      length(raw_costs) != n_regions ||
      any(is.na(log_evidence_lower)) ||
      any(log_evidence_lower == Inf) ||
      any(!is.finite(log_evidence_upper)) ||
      any(log_evidence_lower > log_evidence_upper) ||
      any(!is.finite(raw_costs)) ||
      any(raw_costs < 1) ||
      any(abs(raw_costs - round(raw_costs)) >
          sqrt(.Machine$double.eps)) ||
      length(epsilon) != 1L ||
      !is.finite(epsilon) ||
      epsilon < 0 ||
      epsilon >= 1) {
    stop(
      paste0(
        "Evidence intervals must be ordered, upper endpoints finite, ",
        "costs positive integers, and epsilon in [0, 1). Diagnostics: ",
        "n=", n_regions,
        ", lower_na=", sum(is.na(log_evidence_lower)),
        ", lower_posinf=", sum(log_evidence_lower == Inf, na.rm = TRUE),
        ", upper_nonfinite=", sum(!is.finite(log_evidence_upper)),
        ", inversions=", sum(
          log_evidence_lower > log_evidence_upper, na.rm = TRUE
        ),
        ", invalid_costs=", sum(
          !is.finite(raw_costs) | raw_costs < 1 |
            abs(raw_costs - round(raw_costs)) > sqrt(.Machine$double.eps),
          na.rm = TRUE
        ),
        ", epsilon_length=", length(epsilon),
        ", epsilon_value=", if (length(epsilon) == 1L) epsilon else NA_real_,
        "."
      ),
      call. = FALSE
    )
  }
  if (is.null(eligible)) {
    eligible <- rep(TRUE, n_regions)
  }
  if (length(eligible) != n_regions ||
      any(is.na(eligible)) ||
      !is.logical(eligible)) {
    stop("eligible must be one logical value per region.", call. = FALSE)
  }
  costs <- as.integer(round(raw_costs))
  max_budget <- sum(as.double(costs[eligible]))
  if (!is.finite(max_budget) || max_budget > .Machine$integer.max - 1) {
    stop("The eligible integer-cost budget is too large.", call. = FALSE)
  }
  list(
    log_lower = log_evidence_lower,
    log_upper = log_evidence_upper,
    costs = costs,
    epsilon = epsilon,
    eligible = eligible,
    max_budget = as.integer(max_budget)
  )
}

rrs_log_weighted_pair <- function(log_x, weight_x, log_y, weight_y) {
  term_x <- if (weight_x == 0 || is.infinite(log_x) && log_x < 0) {
    -Inf
  } else {
    log(weight_x) + log_x
  }
  term_y <- if (weight_y == 0 || is.infinite(log_y) && log_y < 0) {
    -Inf
  } else {
    log(weight_y) + log_y
  }
  rrs_log_sum_exp(c(term_x, term_y))
}

rrs_interval_omitted_mass <- function(log_evidence_lower,
                                      log_evidence_upper,
                                      selected) {
  n_regions <- length(log_evidence_lower)
  selected <- sort(unique(as.integer(selected)))
  if (length(selected) < 1L ||
      any(selected < 1L) ||
      any(selected > n_regions)) {
    stop("selected must identify a nonempty region subset.", call. = FALSE)
  }
  complement <- setdiff(seq_len(n_regions), selected)
  log_lower_selected <- rrs_log_sum_exp(log_evidence_lower[selected])
  log_upper_selected <- rrs_log_sum_exp(log_evidence_upper[selected])
  if (length(complement) == 0L) {
    return(list(
      lower = 0,
      upper = 0,
      log_lower_selected = log_lower_selected,
      log_upper_selected = log_upper_selected,
      log_lower_complement = -Inf,
      log_upper_complement = -Inf
    ))
  }
  log_lower_complement <- rrs_log_sum_exp(log_evidence_lower[complement])
  log_upper_complement <- rrs_log_sum_exp(log_evidence_upper[complement])
  lower <- if (is.infinite(log_lower_complement) &&
      log_lower_complement < 0) {
    0
  } else {
    exp(
      log_lower_complement -
        rrs_log_sum_exp(c(log_upper_selected, log_lower_complement))
    )
  }
  upper <- if (is.infinite(log_lower_selected) &&
      log_lower_selected < 0) {
    1
  } else {
    exp(
      log_upper_complement -
        rrs_log_sum_exp(c(log_lower_selected, log_upper_complement))
    )
  }
  list(
    lower = lower,
    upper = upper,
    log_lower_selected = log_lower_selected,
    log_upper_selected = log_upper_selected,
    log_lower_complement = log_lower_complement,
    log_upper_complement = log_upper_complement
  )
}

rrs_interval_epsilon_feasible <- function(interval_summary,
                                          epsilon,
                                          mode) {
  log_epsilon <- if (epsilon == 0) -Inf else log(epsilon)
  log_one_minus_epsilon <- log1p(-epsilon)
  if (identical(mode, "conservative")) {
    # Sharp sufficient condition:
    # U(S^c) / {L(S) + U(S^c)} <= epsilon.
    log_left <- log_one_minus_epsilon +
      interval_summary$log_upper_complement
    log_right <- log_epsilon + interval_summary$log_lower_selected
  } else {
    # Sharp necessary condition:
    # L(S^c) / {U(S) + L(S^c)} <= epsilon.
    # This is stronger than U(S) >= (1-epsilon)L(F).
    log_left <- log_one_minus_epsilon +
      interval_summary$log_lower_complement
    log_right <- log_epsilon + interval_summary$log_upper_selected
  }
  if (is.nan(log_left)) {
    log_left <- -Inf
  }
  if (is.nan(log_right)) {
    log_right <- -Inf
  }
  log_left <= log_right
}

rrs_cumulative_log_sum_exp <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) {
    return(numeric(0))
  }
  out <- numeric(length(x))
  running <- -Inf
  for (j in seq_along(x)) {
    running <- rrs_log_sum_exp(c(running, x[j]))
    out[j] <- running
  }
  out
}

rrs_strict_log_improvement <- function(candidate, incumbent) {
  if (is.infinite(incumbent) && incumbent < 0) {
    return(is.finite(candidate))
  }
  tolerance <- 16 * .Machine$double.eps * max(
    1,
    abs(candidate),
    abs(incumbent)
  )
  candidate > incumbent + tolerance
}

rrs_log_benefit_frontier_generic <- function(log_benefit,
                                              costs,
                                              eligible) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) {
    return(list(
      budget = integer(0),
      log_benefit = numeric(0),
      selected = list()
    ))
  }
  max_budget <- sum(costs[eligible_ids])
  state_limit <- getOption("rrs.max_dynamic_program_states", 5000000L)
  if (length(state_limit) != 1L || !is.finite(state_limit) ||
      state_limit < 2 || state_limit != round(state_limit)) {
    stop("rrs.max_dynamic_program_states must be an integer of at least two.", call. = FALSE)
  }
  if (max_budget + 1 > state_limit) {
    stop(
      paste0(
        "The exact cost dynamic program needs ", max_budget + 1,
        " states, exceeding rrs.max_dynamic_program_states=", state_limit,
        ". Rescale the declared integer costs or raise the audited limit."
      ),
      call. = FALSE
    )
  }
  reachable <- rep(FALSE, max_budget + 1L)
  best <- rep(-Inf, max_budget + 1L)
  selected <- vector("list", max_budget + 1L)
  reachable[1L] <- TRUE
  selected[[1L]] <- integer(0)

  for (j in eligible_ids) {
    for (budget in seq.int(max_budget, costs[j], by = -1L)) {
      previous <- budget - costs[j]
      if (!reachable[previous + 1L]) {
        next
      }
      candidate <- rrs_log_sum_exp(c(
        best[previous + 1L],
        log_benefit[j]
      ))
      if (!reachable[budget + 1L] ||
          rrs_strict_log_improvement(candidate, best[budget + 1L])) {
        reachable[budget + 1L] <- TRUE
        best[budget + 1L] <- candidate
        selected[[budget + 1L]] <- c(selected[[previous + 1L]], j)
      }
    }
  }
  reachable_budget <- which(reachable)[-1L] - 1L
  list(
    budget = reachable_budget,
    log_benefit = best[reachable_budget + 1L],
    selected = selected[reachable_budget + 1L]
  )
}

rrs_log_benefit_frontier_equal_cost <- function(log_benefit,
                                                 costs,
                                                 eligible,
                                                 ranked_ids = NULL) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) {
    return(list(
      budget = integer(0),
      log_benefit = numeric(0),
      selected = list()
    ))
  }
  common_cost <- costs[eligible_ids[1L]]
  if (any(costs[eligible_ids] != common_cost)) {
    stop("The equal-cost frontier received unequal eligible costs.", call. = FALSE)
  }
  if (is.null(ranked_ids)) {
    ranked_ids <- eligible_ids[order(
      -log_benefit[eligible_ids],
      eligible_ids,
      method = "radix"
    )]
  }
  selected <- vector("list", length(ranked_ids))
  running_ids <- integer(0)
  for (k in seq_along(ranked_ids)) {
    insertion <- findInterval(ranked_ids[k], running_ids)
    running_ids <- append(running_ids, ranked_ids[k], after = insertion)
    selected[[k]] <- running_ids
  }
  list(
    budget = as.integer(seq_along(ranked_ids) * common_cost),
    log_benefit = rrs_cumulative_log_sum_exp(log_benefit[ranked_ids]),
    selected = selected
  )
}

rrs_log_benefit_frontier <- function(log_benefit, costs, eligible) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) > 0L &&
      length(unique(costs[eligible_ids])) == 1L) {
    return(rrs_log_benefit_frontier_equal_cost(
      log_benefit = log_benefit,
      costs = costs,
      eligible = eligible
    ))
  }
  rrs_log_benefit_frontier_generic(
    log_benefit = log_benefit,
    costs = costs,
    eligible = eligible
  )
}

region_interval_evidence_frontier <- function(
    log_evidence_lower,
    costs,
    epsilon,
    log_evidence_upper = log_evidence_lower,
    mode = c("conservative", "optimistic"),
    eligible = NULL) {
  mode <- match.arg(mode)
  inputs <- rrs_validate_interval_frontier_inputs(
    log_evidence_lower = log_evidence_lower,
    log_evidence_upper = log_evidence_upper,
    costs = costs,
    epsilon = epsilon,
    eligible = eligible
  )
  log_lower <- inputs$log_lower
  log_upper <- inputs$log_upper
  epsilon <- inputs$epsilon
  retain_weight <- 1 - epsilon
  log_benefit <- vapply(seq_along(log_lower), function(j) {
    if (identical(mode, "conservative")) {
      rrs_log_weighted_pair(
        log_upper[j], retain_weight,
        log_lower[j], epsilon
      )
    } else {
      rrs_log_weighted_pair(
        log_lower[j], retain_weight,
        log_upper[j], epsilon
      )
    }
  }, numeric(1))
  benefit_frontier <- rrs_log_benefit_frontier(
    log_benefit = log_benefit,
    costs = inputs$costs,
    eligible = inputs$eligible
  )
  if (length(benefit_frontier$budget) == 0L) {
    return(data.frame(
      budget = integer(0),
      log_benefit = numeric(0),
      log_evidence_lower = numeric(0),
      log_evidence_upper = numeric(0),
      omitted_mass_lower = numeric(0),
      omitted_mass_upper = numeric(0),
      feasible = logical(0),
      selected = I(list())
    ))
  }

  interval_summaries <- lapply(
    benefit_frontier$selected,
    function(selected) {
      rrs_interval_omitted_mass(log_lower, log_upper, selected)
    }
  )
  omitted_lower <- vapply(interval_summaries, `[[`, numeric(1), "lower")
  omitted_upper <- vapply(interval_summaries, `[[`, numeric(1), "upper")
  feasible <- vapply(interval_summaries, function(summary) {
    rrs_interval_epsilon_feasible(summary, epsilon, mode)
  }, logical(1))
  data.frame(
    budget = benefit_frontier$budget,
    log_benefit = benefit_frontier$log_benefit,
    log_evidence_lower = vapply(
      interval_summaries,
      `[[`,
      numeric(1),
      "log_lower_selected"
    ),
    log_evidence_upper = vapply(
      interval_summaries,
      `[[`,
      numeric(1),
      "log_upper_selected"
    ),
    omitted_mass_lower = omitted_lower,
    omitted_mass_upper = omitted_upper,
    feasible = feasible,
    selected = I(benefit_frontier$selected)
  )
}

rrs_minimum_cost_log_benefit <- function(log_benefit,
                                         costs,
                                         eligible,
                                         log_threshold,
                                         feasibility_check = NULL) {
  eligible_ids <- which(eligible)
  if (length(eligible_ids) == 0L) {
    return(list(feasible = FALSE, selected = integer(0), budget = Inf))
  }
  if (is.infinite(log_threshold) && log_threshold < 0) {
    selected <- eligible_ids[which.min(costs[eligible_ids])]
    return(list(
      feasible = TRUE,
      selected = selected,
      budget = costs[selected],
      log_benefit = log_benefit[selected]
    ))
  }
  max_budget <- sum(costs[eligible_ids])
  state_limit <- getOption("rrs.max_dynamic_program_states", 5000000L)
  if (max_budget + 1 > state_limit) {
    stop(
      paste0(
        "The exact cost dynamic program needs ", max_budget + 1,
        " states, exceeding rrs.max_dynamic_program_states=", state_limit,
        ". Rescale the declared integer costs or raise the audited limit."
      ),
      call. = FALSE
    )
  }
  best <- rep(-Inf, max_budget + 1L)
  best[1L] <- -Inf
  take <- matrix(
    FALSE,
    nrow = length(eligible_ids),
    ncol = max_budget + 1L
  )
  for (position in seq_along(eligible_ids)) {
    item <- eligible_ids[position]
    cost <- costs[item]
    previous <- best
    source <- seq_len(max_budget + 1L - cost)
    target <- source + cost
    left <- previous[source]
    right <- rep(log_benefit[item], length(source))
    candidate <- pmax(left, right) + log1p(exp(-abs(left - right)))
    candidate[is.infinite(left) & left < 0] <- right[is.infinite(left) & left < 0]
    tolerance <- 16 * .Machine$double.eps * pmax(
      1,
      abs(candidate),
      abs(previous[target])
    )
    improve <- candidate > previous[target] + tolerance
    improve[is.infinite(previous[target]) & previous[target] < 0 &
      is.finite(candidate)] <- TRUE
    best[target[improve]] <- candidate[improve]
    take[position, target[improve]] <- TRUE
  }
  feasible_budgets <- which(best >= log_threshold)
  if (length(feasible_budgets) == 0L) {
    if (all(eligible)) {
      selected <- eligible_ids
      if (is.null(feasibility_check) || isTRUE(feasibility_check(selected))) {
        return(list(
          feasible = TRUE,
          selected = selected,
          budget = sum(costs[selected]),
          log_benefit = rrs_log_sum_exp(log_benefit[selected])
        ))
      }
    }
    return(list(feasible = FALSE, selected = integer(0), budget = Inf))
  }
  for (feasible_budget in feasible_budgets) {
    budget_index <- feasible_budget
    selected <- integer(0)
    for (position in rev(seq_along(eligible_ids))) {
      if (take[position, budget_index]) {
        item <- eligible_ids[position]
        selected <- c(item, selected)
        budget_index <- budget_index - costs[item]
      }
    }
    if (is.null(feasibility_check) || isTRUE(feasibility_check(selected))) {
      return(list(
        feasible = TRUE,
        selected = selected,
        budget = feasible_budget - 1L,
        log_benefit = best[feasible_budget]
      ))
    }
  }
  list(feasible = FALSE, selected = integer(0), budget = Inf)
}

select_region_union_epsilon <- function(
    log_evidence_lower,
    costs,
    epsilon,
    log_evidence_upper = log_evidence_lower,
    mode = c("conservative", "optimistic"),
    eligible = NULL) {
  mode <- match.arg(mode)
  inputs <- rrs_validate_interval_frontier_inputs(
    log_evidence_lower = log_evidence_lower,
    log_evidence_upper = log_evidence_upper,
    costs = costs,
    epsilon = epsilon,
    eligible = eligible
  )
  if (inputs$epsilon == 0) {
    required <- if (identical(mode, "conservative")) {
      which(is.finite(inputs$log_upper))
    } else {
      which(is.finite(inputs$log_lower))
    }
    if (any(!inputs$eligible[required])) {
      return(list(
        feasible = FALSE,
        selected = integer(0),
        budget = Inf,
        log_evidence_lower = -Inf,
        log_evidence_upper = -Inf,
        omitted_mass_lower = 0,
        omitted_mass_upper = 1,
        epsilon = inputs$epsilon,
        mode = mode,
        frontier = NULL
      ))
    }
    if (length(required) == 0L) {
      eligible_ids <- which(inputs$eligible)
      if (length(eligible_ids) == 0L) {
        return(list(
          feasible = FALSE,
          selected = integer(0),
          budget = Inf,
          log_evidence_lower = -Inf,
          log_evidence_upper = -Inf,
          omitted_mass_lower = 0,
          omitted_mass_upper = 1,
          epsilon = inputs$epsilon,
          mode = mode,
          frontier = NULL
        ))
      }
      required <- eligible_ids[which.min(inputs$costs[eligible_ids])]
    }
    interval <- rrs_interval_omitted_mass(
      inputs$log_lower,
      inputs$log_upper,
      required
    )
    return(list(
      feasible = TRUE,
      selected = required,
      budget = sum(inputs$costs[required]),
      log_evidence_lower = interval$log_lower_selected,
      log_evidence_upper = interval$log_upper_selected,
      omitted_mass_lower = interval$lower,
      omitted_mass_upper = interval$upper,
      epsilon = inputs$epsilon,
      mode = mode,
      frontier = NULL
    ))
  }
  retain_weight <- 1 - inputs$epsilon
  log_benefit <- vapply(seq_along(inputs$log_lower), function(j) {
    if (identical(mode, "conservative")) {
      rrs_log_weighted_pair(
        inputs$log_upper[j], retain_weight,
        inputs$log_lower[j], inputs$epsilon
      )
    } else {
      rrs_log_weighted_pair(
        inputs$log_lower[j], retain_weight,
        inputs$log_upper[j], inputs$epsilon
      )
    }
  }, numeric(1))
  total_log_evidence <- if (identical(mode, "conservative")) {
    rrs_log_sum_exp(inputs$log_upper)
  } else {
    rrs_log_sum_exp(inputs$log_lower)
  }
  log_threshold <- if (inputs$epsilon == 1) {
    -Inf
  } else {
    log1p(-inputs$epsilon) + total_log_evidence
  }
  optimum <- rrs_minimum_cost_log_benefit(
    log_benefit = log_benefit,
    costs = inputs$costs,
    eligible = inputs$eligible,
    log_threshold = log_threshold,
    feasibility_check = function(selected) {
      interval <- rrs_interval_omitted_mass(
        inputs$log_lower,
        inputs$log_upper,
        selected
      )
      rrs_interval_epsilon_feasible(interval, inputs$epsilon, mode)
    }
  )
  if (!optimum$feasible) {
    return(list(
      feasible = FALSE,
      selected = integer(0),
      budget = Inf,
      log_evidence_lower = -Inf,
      log_evidence_upper = -Inf,
      omitted_mass_lower = 0,
      omitted_mass_upper = 1,
      epsilon = epsilon,
      mode = mode,
      frontier = NULL
    ))
  }
  interval <- rrs_interval_omitted_mass(
    inputs$log_lower,
    inputs$log_upper,
    optimum$selected
  )
  if (!rrs_interval_epsilon_feasible(interval, inputs$epsilon, mode)) {
    stop("The validated exact benefit program returned an infeasible subset.", call. = FALSE)
  }
  list(
    feasible = TRUE,
    selected = optimum$selected,
    budget = optimum$budget,
    log_evidence_lower = interval$log_lower_selected,
    log_evidence_upper = interval$log_upper_selected,
    omitted_mass_lower = interval$lower,
    omitted_mass_upper = interval$upper,
    epsilon = inputs$epsilon,
    mode = mode,
    frontier = NULL
  )
}

select_region_union <- function(log_evidence, costs, lambda) {
  if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be one positive finite number.", call. = FALSE)
  }
  frontier <- region_evidence_frontier(log_evidence, costs)
  frontier$objective <- -frontier$log_evidence + lambda * frontier$budget
  optimum <- which.min(frontier$objective)
  list(
    selected = frontier$selected[[optimum]],
    budget = frontier$budget[optimum],
    log_evidence = frontier$log_evidence[optimum],
    objective = frontier$objective[optimum],
    frontier = frontier
  )
}

regional_omitted_mass_bound <- function(log_evidence,
                                        selected,
                                        log_error_bound) {
  log_evidence <- as.numeric(log_evidence)
  selected <- sort(unique(as.integer(selected)))
  if (length(selected) < 1L ||
      any(selected < 1L) ||
      any(selected > length(log_evidence)) ||
      any(!is.finite(log_evidence)) ||
      length(log_error_bound) != 1L ||
      !is.finite(log_error_bound) ||
      log_error_bound < 0) {
    stop("Invalid omitted-mass-bound inputs.", call. = FALSE)
  }
  complement <- setdiff(seq_along(log_evidence), selected)
  if (length(complement) == 0L) {
    return(0)
  }
  log_selected <- rrs_log_sum_exp(log_evidence[selected])
  log_complement <- rrs_log_sum_exp(log_evidence[complement])
  log_numerator <- log_error_bound + log_complement
  log_denominator <- rrs_log_sum_exp(c(
    -log_error_bound + log_selected,
    log_error_bound + log_complement
  ))
  exp(log_numerator - log_denominator)
}

fit_region_restricted_smc <- function(regions,
                                      log_score_fn,
                                      lambda,
                                      log_error_bound = NULL,
                                      ...) {
  if (length(regions) < 1L) {
    stop("At least one region is required.", call. = FALSE)
  }
  fits <- lapply(
    regions,
    region_annealed_smc,
    log_score_fn = log_score_fn,
    ...
  )
  log_evidence <- vapply(fits, `[[`, numeric(1), "log_evidence")
  costs <- vapply(regions, `[[`, numeric(1), "cost")
  selection <- select_region_union(log_evidence, costs, lambda)
  regional_weights <- exp(
    log_evidence[selection$selected] -
      rrs_log_sum_exp(log_evidence[selection$selected])
  )
  omitted_bound <- if (is.null(log_error_bound)) {
    NA_real_
  } else {
    regional_omitted_mass_bound(
      log_evidence,
      selection$selected,
      log_error_bound
    )
  }

  list(
    selected = selection$selected,
    selected_labels = vapply(
      regions[selection$selected],
      `[[`,
      character(1),
      "label"
    ),
    regional_weights = regional_weights,
    omitted_mass_bound = omitted_bound,
    selection = selection,
    region_fits = fits,
    total_score_evaluations = sum(
      vapply(fits, `[[`, numeric(1), "score_evaluations")
    )
  )
}

rrs_validate_phi_schedule <- function(phi_schedule) {
  phi_schedule <- as.numeric(phi_schedule)
  if (length(phi_schedule) < 2L ||
      phi_schedule[1L] != 0 ||
      tail(phi_schedule, 1L) != 1 ||
      any(!is.finite(phi_schedule)) ||
      any(diff(phi_schedule) <= 0)) {
    stop(
      "phi_schedule must increase strictly from zero to one.",
      call. = FALSE
    )
  }
  phi_schedule
}

rrs_refine_phi_schedule <- function(phi_schedule,
                                    diagnostics,
                                    n_particles,
                                    ess_fraction = 0.5,
                                    max_increment = 0.05) {
  phi_schedule <- rrs_validate_phi_schedule(phi_schedule)
  diagnostics <- as.data.frame(diagnostics)
  n_particles <- as.integer(n_particles)
  if (nrow(diagnostics) != length(phi_schedule) - 1L ||
      !all(c("phi_previous", "phi", "ess") %in% names(diagnostics)) ||
      n_particles < 2L ||
      !(ess_fraction > 0 && ess_fraction < 1) ||
      !(max_increment > 0 && max_increment <= 1)) {
    stop("Invalid fixed-schedule refinement inputs.", call. = FALSE)
  }
  if (any(abs(diagnostics$phi_previous - head(phi_schedule, -1L)) >
          sqrt(.Machine$double.eps)) ||
      any(abs(diagnostics$phi - tail(phi_schedule, -1L)) >
          sqrt(.Machine$double.eps))) {
    stop("Diagnostics do not match the declared temperature schedule.", call. = FALSE)
  }
  increments <- diff(phi_schedule)
  pieces <- pmax(1L, ceiling(increments / max_increment))
  pieces[diagnostics$ess < ess_fraction * n_particles] <- pmax(
    pieces[diagnostics$ess < ess_fraction * n_particles],
    2L
  )
  refined <- phi_schedule[1L]
  for (stage in seq_along(increments)) {
    local <- seq(
      phi_schedule[stage],
      phi_schedule[stage + 1L],
      length.out = pieces[stage] + 1L
    )
    refined <- c(refined, local[-1L])
  }
  sort(unique(refined))
}

region_annealed_importance_sampling <- function(region,
                                                 log_score_fn,
                                                 phi_schedule,
                                                 n_paths = 1000L,
                                                 mutation_steps = 2L) {
  rrs_validate_region(region)
  phi_schedule <- rrs_validate_phi_schedule(phi_schedule)
  n_paths <- as.integer(n_paths)
  mutation_steps <- as.integer(mutation_steps)
  if (n_paths < 2L) {
    stop("n_paths must be at least two.", call. = FALSE)
  }
  if (mutation_steps < 0L) {
    stop("mutation_steps must be nonnegative.", call. = FALSE)
  }

  particles <- as.matrix(region$sample(n_paths))
  log_nu <- region$log_prob(particles)
  if (nrow(particles) != n_paths || any(!is.finite(log_nu))) {
    stop("The regional proposal returned invalid paths.", call. = FALSE)
  }
  log_score <- rrs_evaluate_log_score(log_score_fn, particles)
  log_g <- log_score - log_nu
  log_path_weights <- rep(0, n_paths)
  score_evaluations <- n_paths
  diagnostic_rows <- vector("list", length(phi_schedule) - 1L)

  for (stage in seq_len(length(phi_schedule) - 1L)) {
    phi_previous <- phi_schedule[stage]
    phi <- phi_schedule[stage + 1L]
    delta <- phi - phi_previous
    log_path_weights <- log_path_weights + delta * log_g
    mutation <- rrs_reversible_mutation(
      particles = particles,
      log_g = log_g,
      phi = phi,
      region = region,
      log_score_fn = log_score_fn,
      n_steps = mutation_steps
    )
    particles <- mutation$particles
    log_g <- mutation$log_g
    score_evaluations <- score_evaluations + mutation$score_evaluations
    diagnostic_rows[[stage]] <- data.frame(
      stage = stage,
      phi_previous = phi_previous,
      phi = phi,
      delta = delta,
      path_ess = rrs_effective_sample_size(log_path_weights),
      mutation_type = mutation$type,
      mutation_acceptance = if (mutation$proposed > 0L) {
        mutation$accepted / mutation$proposed
      } else {
        NA_real_
      }
    )
  }

  log_evidence <- rrs_log_sum_exp(log_path_weights) - log(n_paths)
  list(
    label = region$label,
    log_evidence = log_evidence,
    evidence = exp(log_evidence),
    log_path_weights = log_path_weights,
    particles = particles,
    phi_schedule = phi_schedule,
    diagnostics = do.call(rbind, diagnostic_rows),
    score_evaluations = score_evaluations,
    path_ess = rrs_effective_sample_size(log_path_weights)
  )
}

rrs_evidence_interval <- function(log_estimate,
                                  rho,
                                  n_paths,
                                  n_declared_nodes,
                                  delta) {
  if (length(log_estimate) != 1L ||
      !is.finite(log_estimate) ||
      length(rho) != 1L ||
      !is.finite(rho) ||
      rho < 0 ||
      n_paths < 1L ||
      n_declared_nodes < 1L ||
      !(delta > 0 && delta < 1)) {
    stop("Invalid evidence-interval inputs.", call. = FALSE)
  }
  relative_radius <- if (rho == 0) {
    0
  } else {
    (exp(rho) - 1) *
      sqrt(log(2 * n_declared_nodes / delta) / (2 * n_paths))
  }
  if (relative_radius >= 1) {
    return(list(
      log_lower = -Inf,
      log_upper = Inf,
      relative_radius = relative_radius,
      finite = FALSE
    ))
  }
  list(
    log_lower = log_estimate - log1p(relative_radius),
    log_upper = log_estimate - log1p(-relative_radius),
    relative_radius = relative_radius,
    finite = TRUE
  )
}

pilot_fixed_region_evidence <- function(region,
                                        log_score_fn,
                                        rho,
                                        n_declared_nodes,
                                        delta = 0.05,
                                        n_pilot = 250L,
                                        n_paths = 1000L,
                                        ess_fraction = 0.8,
                                        resample_fraction = 0.5,
                                        pilot_mutation_steps = 1L,
                                        production_mutation_steps = 2L) {
  pilot <- region_annealed_smc(
    region = region,
    log_score_fn = log_score_fn,
    n_particles = n_pilot,
    ess_fraction = ess_fraction,
    resample_fraction = resample_fraction,
    mutation_steps = pilot_mutation_steps
  )
  phi_schedule <- c(0, pilot$diagnostics$phi)
  production <- region_annealed_importance_sampling(
    region = region,
    log_score_fn = log_score_fn,
    phi_schedule = phi_schedule,
    n_paths = n_paths,
    mutation_steps = production_mutation_steps
  )
  interval <- rrs_evidence_interval(
    log_estimate = production$log_evidence,
    rho = rho,
    n_paths = n_paths,
    n_declared_nodes = n_declared_nodes,
    delta = delta
  )
  list(
    log_evidence = production$log_evidence,
    evidence = production$evidence,
    log_lower = interval$log_lower,
    log_upper = interval$log_upper,
    relative_radius = interval$relative_radius,
    finite_interval = interval$finite,
    phi_schedule = phi_schedule,
    pilot = pilot,
    production = production,
    score_evaluations =
      pilot$score_evaluations + production$score_evaluations
  )
}

iid_bounded_region_evidence <- function(region,
                                        log_score_fn,
                                        rho,
                                        n_declared_nodes,
                                        delta = 0.05,
                                        n_paths = 1000L) {
  rrs_validate_region(region)
  if (length(n_paths) != 1L || !is.finite(n_paths) || n_paths < 1L) {
    stop("n_paths must be one positive finite integer.", call. = FALSE)
  }
  n_paths <- as.integer(n_paths)
  particles <- region$sample(n_paths)
  log_score <- rrs_evaluate_log_score(log_score_fn, particles)
  log_proposal <- region$log_prob_inside(particles)
  if (length(log_proposal) != n_paths || any(!is.finite(log_proposal))) {
    stop("The exact region proposal returned invalid log probabilities.", call. = FALSE)
  }
  log_weights <- log_score - log_proposal
  log_evidence <- rrs_log_sum_exp(log_weights) - log(n_paths)
  interval <- rrs_evidence_interval(
    log_estimate = log_evidence,
    rho = rho,
    n_paths = n_paths,
    n_declared_nodes = n_declared_nodes,
    delta = delta
  )
  list(
    log_evidence = log_evidence,
    evidence = exp(log_evidence),
    log_lower = interval$log_lower,
    log_upper = interval$log_upper,
    relative_radius = interval$relative_radius,
    finite_interval = interval$finite,
    particles = particles,
    log_weights = log_weights,
    score_evaluations = n_paths
  )
}

new_region_tree_node <- function(id,
                                 region,
                                 children = character(0),
                                 rho,
                                 cost_lower = NULL,
                                 terminal = length(children) == 0L) {
  rrs_validate_region(region)
  if (length(id) != 1L || !nzchar(id)) {
    stop("A tree node needs one nonempty id.", call. = FALSE)
  }
  if (length(rho) != 1L || !is.finite(rho) || rho < 0) {
    stop("A tree node needs a finite nonnegative rho.", call. = FALSE)
  }
  if (abs(region$cost - round(region$cost)) >
      sqrt(.Machine$double.eps)) {
    stop("Region-tree costs must be positive integers.", call. = FALSE)
  }
  region$cost <- as.integer(round(region$cost))
  if (is.null(cost_lower)) {
    cost_lower <- region$cost
  }
  if (length(cost_lower) != 1L ||
      !is.finite(cost_lower) ||
      cost_lower < 1 ||
      abs(cost_lower - round(cost_lower)) > sqrt(.Machine$double.eps)) {
    stop("cost_lower must be a positive integer.", call. = FALSE)
  }
  if (isTRUE(terminal) && cost_lower > region$cost) {
    stop(
      "A terminal cost lower bound cannot exceed the terminal region cost.",
      call. = FALSE
    )
  }
  list(
    id = as.character(id),
    region = region,
    children = as.character(children),
    rho = rho,
    cost_lower = as.integer(round(cost_lower)),
    terminal = isTRUE(terminal),
    evaluated = FALSE,
    log_lower = -Inf,
    log_upper = Inf,
    bound_score_evaluations = 0L,
    evidence_fit = NULL
  )
}

new_region_tree <- function(nodes, root_id) {
  if (length(nodes) < 1L) {
    stop("A region tree must contain at least one node.", call. = FALSE)
  }
  ids <- vapply(nodes, `[[`, character(1), "id")
  if (anyDuplicated(ids)) {
    stop("Region-tree node ids must be unique.", call. = FALSE)
  }
  names(nodes) <- ids
  if (!root_id %in% ids) {
    stop("root_id is not a node in the tree.", call. = FALSE)
  }
  for (node in nodes) {
    if (any(!node$children %in% ids)) {
      stop("Every child id must name a declared tree node.", call. = FALSE)
    }
    if (node$terminal && length(node$children) > 0L) {
      stop("A terminal node cannot have children.", call. = FALSE)
    }
    if (!node$terminal && length(node$children) == 0L) {
      stop("A nonterminal node must have children.", call. = FALSE)
    }
  }
  child_ids <- unlist(lapply(nodes, `[[`, "children"), use.names = FALSE)
  if (anyDuplicated(child_ids)) {
    stop("Every nonroot tree node must have exactly one parent.", call. = FALSE)
  }
  if (root_id %in% child_ids) {
    stop("The tree root cannot also be a child.", call. = FALSE)
  }
  nonroot_ids <- setdiff(ids, root_id)
  if (!setequal(child_ids, nonroot_ids)) {
    stop(
      "Every declared nonroot node must be reachable from exactly one parent.",
      call. = FALSE
    )
  }
  visited <- character(0)
  queue <- root_id
  while (length(queue) > 0L) {
    current <- queue[1L]
    queue <- queue[-1L]
    if (current %in% visited) {
      stop("The declared region tree contains a cycle.", call. = FALSE)
    }
    visited <- c(visited, current)
    queue <- c(queue, nodes[[current]]$children)
  }
  if (!setequal(visited, ids)) {
    stop("Every declared node must be reachable from the root.", call. = FALSE)
  }
  structure(
    list(nodes = nodes, root_id = root_id, n_declared_nodes = length(nodes)),
    class = "rrs_region_tree"
  )
}

new_lazy_region_tree <- function(root_node,
                                 n_declared_nodes,
                                 expand_node) {
  if (!is.function(expand_node) ||
      length(n_declared_nodes) != 1L ||
      !is.finite(n_declared_nodes) ||
      n_declared_nodes < 1 ||
      abs(n_declared_nodes - round(n_declared_nodes)) >
        sqrt(.Machine$double.eps)) {
    stop("A lazy tree needs a finite declared size and an expansion function.", call. = FALSE)
  }
  if (root_node$terminal && n_declared_nodes != 1) {
    stop("A terminal lazy root can only declare one node.", call. = FALSE)
  }
  nodes <- list(root_node)
  names(nodes) <- root_node$id
  structure(
    list(
      nodes = nodes,
      root_id = root_node$id,
      n_declared_nodes = as.numeric(n_declared_nodes),
      expand_node = expand_node
    ),
    class = c("rrs_lazy_region_tree", "rrs_region_tree")
  )
}

rrs_anchor_evidence_bounds <- function(node, log_score_fn) {
  region <- node$region
  anchor <- if (is.function(region$anchor)) {
    as.matrix(region$anchor())
  } else {
    as.matrix(region$sample(1L))
  }
  if (nrow(anchor) != 1L) {
    stop("A regional anchor must contain exactly one support.", call. = FALSE)
  }
  log_nu <- region$log_prob(anchor)
  if (length(log_nu) != 1L || !is.finite(log_nu)) {
    stop("The regional anchor has invalid proposal probability.", call. = FALSE)
  }
  anchor_log_score <- rrs_evaluate_log_score(log_score_fn, anchor)
  log_g <- anchor_log_score - log_nu
  log_lower <- log_g - node$rho
  log_upper <- log_g + node$rho
  if (is.function(region$log_g_bounds)) {
    direct <- as.numeric(region$log_g_bounds())
    if (length(direct) != 2L ||
        any(!is.finite(direct)) ||
        direct[1L] > direct[2L]) {
      stop("Regional log-g bounds must be two finite ordered values.", call. = FALSE)
    }
    log_lower <- max(log_lower, direct[1L], anchor_log_score)
    log_upper <- min(log_upper, direct[2L])
    if (log_lower > log_upper + 1e-10) {
      stop("Regional structural evidence bounds are inconsistent.", call. = FALSE)
    }
  }
  if (is.function(region$log_evidence_bounds)) {
    direct_evidence <- as.numeric(region$log_evidence_bounds())
    if (length(direct_evidence) != 2L ||
        any(!is.finite(direct_evidence)) ||
        direct_evidence[1L] > direct_evidence[2L]) {
      stop(
        "Regional log-evidence bounds must be two finite ordered values.",
        call. = FALSE
      )
    }
    log_lower <- max(log_lower, direct_evidence[1L], anchor_log_score)
    log_upper <- min(log_upper, direct_evidence[2L])
    if (log_lower > log_upper + 1e-10) {
      stop("Regional evidence bounds are inconsistent.", call. = FALSE)
    }
  }
  list(
    log_lower = log_lower,
    log_upper = log_upper,
    anchor = anchor,
    score_evaluations = 1L
  )
}

rrs_bound_tree_node <- function(node, log_score_fn) {
  bounds <- rrs_anchor_evidence_bounds(node, log_score_fn)
  node$log_lower <- bounds$log_lower
  node$log_upper <- bounds$log_upper
  node$bound_score_evaluations <- bounds$score_evaluations
  node$anchor <- bounds$anchor
  exact_from_bounds <- isTRUE(all.equal(
    bounds$log_lower,
    bounds$log_upper,
    tolerance = 1e-12
  ))
  if (node$terminal && exact_from_bounds) {
    node$evaluated <- TRUE
    node$evidence_fit <- list(
      log_evidence = bounds$log_lower,
      evidence = exp(bounds$log_lower),
      log_lower = bounds$log_lower,
      log_upper = bounds$log_upper,
      relative_radius = 0,
      finite_interval = TRUE,
      phi_schedule = c(0, 1),
      exact_from_anchor = node$rho == 0,
      exact_from_structural_bounds = TRUE,
      score_evaluations = 0L
    )
  }
  node
}

rrs_log_omitted_mass_bound <- function(log_selected_lower,
                                       log_omitted_upper) {
  if (length(log_selected_lower) != 1L ||
      length(log_omitted_upper) != 1L) {
    stop("Omitted-mass inputs must be scalar.", call. = FALSE)
  }
  if (is.infinite(log_omitted_upper) && log_omitted_upper < 0) {
    return(0)
  }
  if (is.infinite(log_selected_lower) && log_selected_lower < 0) {
    return(1)
  }
  if (!is.finite(log_selected_lower) || !is.finite(log_omitted_upper)) {
    stop("Omitted-mass log bounds must be finite or negative infinity.", call. = FALSE)
  }
  exp(
    log_omitted_upper -
      rrs_log_sum_exp(c(log_selected_lower, log_omitted_upper))
  )
}

rrs_tree_search_state <- function(nodes, frontier_ids, lambda) {
  frontier <- nodes[frontier_ids]
  log_upper <- vapply(frontier, `[[`, numeric(1), "log_upper")
  costs_lower <- vapply(frontier, `[[`, numeric(1), "cost_lower")
  if (any(!is.finite(log_upper))) {
    stop("Every frontier node needs a finite evidence upper bound.", call. = FALSE)
  }
  optimistic <- select_region_union(log_upper, costs_lower, lambda)
  optimistic_ids <- frontier_ids[optimistic$selected]

  evaluated <- vapply(frontier, `[[`, logical(1), "evaluated")
  evaluated_ids <- frontier_ids[evaluated]
  incumbent <- NULL
  selected_ids <- character(0)
  if (length(evaluated_ids) > 0L) {
    evaluated_nodes <- nodes[evaluated_ids]
    log_lower <- vapply(evaluated_nodes, `[[`, numeric(1), "log_lower")
    finite_lower <- is.finite(log_lower)
    evaluated_ids <- evaluated_ids[finite_lower]
    evaluated_nodes <- evaluated_nodes[finite_lower]
    log_lower <- log_lower[finite_lower]
    if (length(evaluated_ids) > 0L) {
      incumbent <- select_region_union(
        log_lower,
        vapply(evaluated_nodes, function(node) node$region$cost, numeric(1)),
        lambda
      )
      selected_ids <- evaluated_ids[incumbent$selected]
    }
  }

  objective_gap <- if (is.null(incumbent)) {
    Inf
  } else {
    max(0, incumbent$objective - optimistic$objective)
  }
  omitted_ids <- setdiff(frontier_ids, selected_ids)
  log_selected_lower <- if (length(selected_ids) > 0L) {
    rrs_log_sum_exp(vapply(
      nodes[selected_ids],
      `[[`,
      numeric(1),
      "log_lower"
    ))
  } else {
    -Inf
  }
  log_omitted_upper <- if (length(omitted_ids) > 0L) {
    rrs_log_sum_exp(vapply(
      nodes[omitted_ids],
      `[[`,
      numeric(1),
      "log_upper"
    ))
  } else {
    -Inf
  }
  list(
    optimistic = optimistic,
    optimistic_ids = optimistic_ids,
    incumbent = incumbent,
    selected_ids = selected_ids,
    objective_gap = objective_gap,
    omitted_mass_bound = rrs_log_omitted_mass_bound(
      log_selected_lower,
      log_omitted_upper
    )
  )
}

rrs_tree_epsilon_search_state <- function(nodes, frontier_ids, epsilon) {
  if (!is.list(nodes) ||
      is.null(names(nodes)) ||
      length(frontier_ids) < 1L ||
      any(!frontier_ids %in% names(nodes))) {
    stop("Invalid nodes or frontier ids for epsilon search.", call. = FALSE)
  }
  frontier <- nodes[frontier_ids]
  log_lower <- vapply(frontier, `[[`, numeric(1), "log_lower")
  log_upper <- vapply(frontier, `[[`, numeric(1), "log_upper")
  costs_lower <- vapply(frontier, function(node) {
    if (isTRUE(node$terminal)) node$region$cost else node$cost_lower
  }, numeric(1))

  # The optimistic problem uses the sharp necessary ratio condition
  # L(S^c)/{U(S)+L(S^c)} <= epsilon.  Its optimum is therefore a valid
  # lower bound on the minimum cost of any feasible terminal union.
  optimistic <- select_region_union_epsilon(
    log_evidence_lower = log_lower,
    log_evidence_upper = log_upper,
    costs = costs_lower,
    epsilon = epsilon,
    mode = "optimistic"
  )
  if (!optimistic$feasible) {
    stop(
      "The full frontier must be optimistically epsilon-feasible.",
      call. = FALSE
    )
  }

  evaluated_terminal <- vapply(frontier, function(node) {
    isTRUE(node$evaluated) && isTRUE(node$terminal)
  }, logical(1))
  incumbent_costs <- vapply(frontier, function(node) {
    if (isTRUE(node$evaluated) && isTRUE(node$terminal)) {
      node$region$cost
    } else {
      node$cost_lower
    }
  }, numeric(1))
  incumbent <- select_region_union_epsilon(
    log_evidence_lower = log_lower,
    log_evidence_upper = log_upper,
    costs = incumbent_costs,
    epsilon = epsilon,
    mode = "conservative",
    eligible = evaluated_terminal
  )
  cost_lower <- optimistic$budget
  cost_upper <- if (incumbent$feasible) incumbent$budget else Inf
  bracket_consistent <- !is.finite(cost_upper) || cost_lower <= cost_upper
  cost_gap <- if (is.finite(cost_upper) && bracket_consistent) {
    cost_upper - cost_lower
  } else {
    Inf
  }
  selected_ids <- if (incumbent$feasible) {
    frontier_ids[incumbent$selected]
  } else {
    character(0)
  }
  list(
    optimistic = optimistic,
    optimistic_ids = frontier_ids[optimistic$selected],
    incumbent = incumbent,
    selected_ids = selected_ids,
    cost_lower = cost_lower,
    cost_upper = cost_upper,
    cost_gap = cost_gap,
    bracket_consistent = bracket_consistent,
    omitted_mass_bound = if (incumbent$feasible) {
      incumbent$omitted_mass_upper
    } else {
      1
    },
    optimistic_omitted_mass_lower = optimistic$omitted_mass_lower
  )
}

make_pilot_fixed_evidence_evaluator <- function(log_score_fn,
                                                n_declared_nodes,
                                                delta = 0.05,
                                                n_pilot = 250L,
                                                n_paths = 1000L,
                                                ess_fraction = 0.8,
                                                resample_fraction = 0.5,
                                                pilot_mutation_steps = 1L,
                                                production_mutation_steps = 2L,
                                                maximum_relative_radius = NULL,
                                                maximum_paths = NULL) {
  force(log_score_fn)
  function(node, remaining_score_budget = Inf) {
    scheduled_paths <- n_paths
    if (!is.null(maximum_relative_radius)) {
      path_cap <- if (is.null(maximum_paths)) n_paths else maximum_paths
      preflight <- regional_evidence_preflight(
        rho = node$rho,
        relative_error = maximum_relative_radius,
        delta_node = delta / n_declared_nodes,
        remaining_paths = path_cap
      )
      if (!isTRUE(preflight$feasible)) {
        return(list(
          preflight_failed = TRUE,
          preflight = preflight,
          score_evaluations = 0
        ))
      }
      scheduled_paths <- max(n_paths, preflight$required_paths)
      preflight$scheduled_paths <- scheduled_paths
    } else {
      preflight <- NULL
    }
    budgeted_log_score <- rrs_budgeted_log_score(
      log_score_fn,
      remaining_score_budget
    )
    fit <- tryCatch(pilot_fixed_region_evidence(
      region = node$region,
      log_score_fn = budgeted_log_score,
      rho = node$rho,
      n_declared_nodes = n_declared_nodes,
      delta = delta,
      n_pilot = n_pilot,
      n_paths = scheduled_paths,
      ess_fraction = ess_fraction,
      resample_fraction = resample_fraction,
      pilot_mutation_steps = pilot_mutation_steps,
      production_mutation_steps = production_mutation_steps
    ), rrs_score_budget_exhausted = function(condition) condition)
    used <- attr(budgeted_log_score, "score_budget_usage")()
    if (inherits(fit, "rrs_score_budget_exhausted")) {
      return(list(
        budget_exhausted = TRUE,
        score_evaluations = used,
        budget_condition = fit,
        preflight = preflight
      ))
    }
    if (fit$score_evaluations != used) {
      stop("Evidence score accounting disagrees with the budget guard.", call. = FALSE)
    }
    list(
      log_evidence = fit$log_evidence,
      log_lower = fit$log_lower,
      log_upper = fit$log_upper,
      score_evaluations = fit$score_evaluations,
      fit = fit,
      preflight = preflight
    )
  }
}

make_iid_bounded_evidence_evaluator <- function(log_score_fn,
                                                n_declared_nodes,
                                                delta = 0.05,
                                                n_paths = 1000L,
                                                maximum_relative_radius = NULL,
                                                maximum_paths = NULL) {
  force(log_score_fn)
  function(node, remaining_score_budget = Inf) {
    scheduled_paths <- n_paths
    if (!is.null(maximum_relative_radius)) {
      path_cap <- if (is.null(maximum_paths)) n_paths else maximum_paths
      preflight <- regional_evidence_preflight(
        rho = node$rho,
        relative_error = maximum_relative_radius,
        delta_node = delta / n_declared_nodes,
        remaining_paths = path_cap
      )
      if (!isTRUE(preflight$feasible)) {
        return(list(
          preflight_failed = TRUE,
          preflight = preflight,
          score_evaluations = 0
        ))
      }
      scheduled_paths <- max(n_paths, preflight$required_paths)
      preflight$scheduled_paths <- scheduled_paths
    } else {
      preflight <- NULL
    }
    budgeted_log_score <- rrs_budgeted_log_score(
      log_score_fn,
      remaining_score_budget
    )
    fit <- tryCatch(iid_bounded_region_evidence(
      region = node$region,
      log_score_fn = budgeted_log_score,
      rho = node$rho,
      n_declared_nodes = n_declared_nodes,
      delta = delta,
      n_paths = scheduled_paths
    ), rrs_score_budget_exhausted = function(condition) condition)
    used <- attr(budgeted_log_score, "score_budget_usage")()
    if (inherits(fit, "rrs_score_budget_exhausted")) {
      return(list(
        budget_exhausted = TRUE,
        score_evaluations = used,
        budget_condition = fit,
        preflight = preflight
      ))
    }
    if (fit$score_evaluations != used) {
      stop("Evidence score accounting disagrees with the budget guard.", call. = FALSE)
    }
    list(
      log_evidence = fit$log_evidence,
      log_lower = fit$log_lower,
      log_upper = fit$log_upper,
      score_evaluations = fit$score_evaluations,
      fit = fit,
      preflight = preflight
    )
  }
}

certified_tree_region_search <- function(tree,
                                         log_score_fn,
                                         evidence_evaluator,
                                         lambda = NULL,
                                         objective_tolerance = 0.05,
                                         omitted_mass_tolerance = 0.05,
                                         max_score_evaluations = Inf,
                                         max_actions = Inf,
                                         selection_mode = c(
                                           "penalized",
                                           "epsilon"
                                         ),
                                         epsilon = NULL,
                                         cost_tolerance = 0) {
  selection_mode <- match.arg(selection_mode)
  if (!inherits(tree, "rrs_region_tree")) {
    stop("tree must be created by new_region_tree.", call. = FALSE)
  }
  if (!is.function(evidence_evaluator) || !is.function(log_score_fn)) {
    stop("Score and evidence evaluators must be functions.", call. = FALSE)
  }
  if (!(objective_tolerance >= 0) ||
      !(omitted_mass_tolerance >= 0 && omitted_mass_tolerance < 1) ||
      !(max_score_evaluations > 0) ||
      !(max_actions > 0)) {
    stop("Invalid certified-search tolerances or budgets.", call. = FALSE)
  }
  if (identical(selection_mode, "penalized")) {
    if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
      stop("lambda must be one positive finite number.", call. = FALSE)
    }
  } else {
    if (length(epsilon) != 1L ||
        !is.finite(epsilon) ||
        epsilon < 0 ||
        epsilon >= 1 ||
        length(cost_tolerance) != 1L ||
        !is.finite(cost_tolerance) ||
        cost_tolerance < 0) {
      stop(
        paste0(
          "Epsilon search requires epsilon in [0, 1) and ",
          "nonnegative cost_tolerance."
        ),
        call. = FALSE
      )
    }
  }

  nodes <- tree$nodes
  root_id <- tree$root_id
  nodes[[root_id]] <- rrs_bound_tree_node(nodes[[root_id]], log_score_fn)
  frontier_ids <- root_id
  score_evaluations <- nodes[[root_id]]$bound_score_evaluations
  bound_score_evaluations <- nodes[[root_id]]$bound_score_evaluations
  evidence_score_evaluations <- 0L
  actions <- 0L
  expansions <- 0L
  trace <- list()
  status <- "search_incomplete"
  conflict_node_id <- NA_character_

  repeat {
    state <- if (identical(selection_mode, "penalized")) {
      rrs_tree_search_state(nodes, frontier_ids, lambda)
    } else {
      rrs_tree_epsilon_search_state(nodes, frontier_ids, epsilon)
    }
    common_trace <- list(
      action = actions,
      frontier_size = length(frontier_ids),
      evaluated_nodes = sum(vapply(
        nodes[frontier_ids],
        `[[`,
        logical(1),
        "evaluated"
      ))
    )
    trace[[length(trace) + 1L]] <- if (
        identical(selection_mode, "penalized")) {
      data.frame(
        common_trace,
        objective_lower = state$optimistic$objective,
        objective_upper = if (is.null(state$incumbent)) {
          Inf
        } else {
          state$incumbent$objective
        },
        objective_gap = state$objective_gap,
        omitted_mass_bound = state$omitted_mass_bound,
        score_evaluations = score_evaluations
      )
    } else {
      data.frame(
        common_trace,
        cost_lower = state$cost_lower,
        cost_upper = state$cost_upper,
        cost_gap = state$cost_gap,
        omitted_mass_bound = state$omitted_mass_bound,
        score_evaluations = score_evaluations
      )
    }
    bounds_met <- if (identical(selection_mode, "penalized")) {
      is.finite(state$objective_gap) &&
        state$objective_gap <= objective_tolerance &&
        state$omitted_mass_bound <= omitted_mass_tolerance
    } else {
      isTRUE(state$incumbent$feasible) &&
        state$omitted_mass_bound <= epsilon &&
        is.finite(state$cost_gap) &&
        state$cost_gap <= cost_tolerance
    }
    if (bounds_met) {
      status <- "certified"
      break
    }
    if (score_evaluations >= max_score_evaluations ||
        actions >= max_actions) {
      status <- "budget_exhausted"
      break
    }

    evaluated_frontier <- vapply(
      nodes[frontier_ids],
      `[[`,
      logical(1),
      "evaluated"
    )
    preflight_blocked <- vapply(nodes[frontier_ids], function(node) {
      isTRUE(node$preflight_blocked)
    }, logical(1))
    unresolved_ids <- frontier_ids[!evaluated_frontier & !preflight_blocked]
    if (length(unresolved_ids) == 0L) {
      status <- if (any(preflight_blocked) &&
          !isTRUE(state$incumbent$feasible)) {
        "uninformative_evidence"
      } else {
        "frontier_exhausted"
      }
      break
    }
    preferred <- intersect(state$optimistic_ids, unresolved_ids)
    candidates <- if (length(preferred) > 0L) preferred else unresolved_ids
    priorities <- vapply(
      nodes[candidates],
      function(node) {
        if (identical(selection_mode, "penalized")) {
          node$log_upper - lambda * node$cost_lower
        } else {
          node$log_upper - log(node$cost_lower)
        }
      },
      numeric(1)
    )
    node_id <- candidates[which.max(priorities)]
    node <- nodes[[node_id]]

    if (!node$terminal) {
      if (length(node$children) == 0L &&
          inherits(tree, "rrs_lazy_region_tree")) {
        generated <- tree$expand_node(node)
        if (length(generated) < 1L) {
          stop("A lazy nonterminal expansion returned no children.", call. = FALSE)
        }
        generated_ids <- vapply(generated, `[[`, character(1), "id")
        if (anyDuplicated(generated_ids) ||
            any(generated_ids %in% names(nodes))) {
          stop("A lazy expansion returned duplicate node ids.", call. = FALSE)
        }
        generated_cost_lower <- vapply(
          generated,
          `[[`,
          numeric(1),
          "cost_lower"
        )
        if (node$cost_lower > min(generated_cost_lower)) {
          stop(
            paste0(
              "A lazy parent descendant-cost lower bound cannot exceed the ",
              "smallest generated child descendant-cost lower bound."
            ),
            call. = FALSE
          )
        }
        node$children <- generated_ids
        nodes[[node_id]] <- node
        for (generated_node in generated) {
          nodes[[generated_node$id]] <- generated_node
        }
      }
      child_ids <- node$children
      child_budget_exhausted <- FALSE
      for (child_id in child_ids) {
        if (!is.finite(nodes[[child_id]]$log_upper)) {
          if (score_evaluations >= max_score_evaluations) {
            child_budget_exhausted <- TRUE
            break
          }
          nodes[[child_id]] <- rrs_bound_tree_node(
            nodes[[child_id]],
            log_score_fn
          )
          score_evaluations <- score_evaluations +
            nodes[[child_id]]$bound_score_evaluations
          bound_score_evaluations <- bound_score_evaluations +
            nodes[[child_id]]$bound_score_evaluations
        }
      }
      if (child_budget_exhausted) {
        status <- "budget_exhausted"
        break
      }
      frontier_ids <- c(
        frontier_ids[frontier_ids != node_id],
        child_ids
      )
      expansions <- expansions + 1L
    } else {
      evaluator_formals <- names(formals(evidence_evaluator))
      supports_budget <- "remaining_score_budget" %in% evaluator_formals ||
        "..." %in% evaluator_formals
      fit <- if (supports_budget) {
        evidence_evaluator(
          node,
          remaining_score_budget = max_score_evaluations - score_evaluations
        )
      } else {
        evidence_evaluator(node)
      }
      if (isTRUE(fit$preflight_failed)) {
        node$preflight_blocked <- TRUE
        node$evidence_preflight <- fit$preflight
        node$evaluated <- TRUE
        node$evidence_fit <- list(
          method = "deterministic_structural_fallback",
          log_lower = node$log_lower,
          log_upper = node$log_upper
        )
        nodes[[node_id]] <- node
        actions <- actions + 1L
        next
      }
      if (isTRUE(fit$budget_exhausted)) {
        score_evaluations <- score_evaluations + fit$score_evaluations
        evidence_score_evaluations <- evidence_score_evaluations +
          fit$score_evaluations
        status <- "budget_exhausted"
        break
      }
      required <- c(
        "log_evidence",
        "log_lower",
        "log_upper",
        "score_evaluations"
      )
      if (any(!required %in% names(fit))) {
        stop("The evidence evaluator returned an incomplete result.", call. = FALSE)
      }
      combined_lower <- max(node$log_lower, fit$log_lower)
      combined_upper <- min(node$log_upper, fit$log_upper)
      score_evaluations <- score_evaluations + fit$score_evaluations
      evidence_score_evaluations <- evidence_score_evaluations +
        fit$score_evaluations
      if (combined_lower > combined_upper + 1e-10) {
        node$evidence_fit <- fit$fit
        node$interval_conflict <- TRUE
        nodes[[node_id]] <- node
        actions <- actions + 1L
        status <- "evidence_interval_conflict"
        conflict_node_id <- node_id
        break
      }
      node$log_lower <- combined_lower
      node$log_upper <- combined_upper
      if (node$log_lower > node$log_upper) {
        midpoint <- (node$log_lower + node$log_upper) / 2
        node$log_lower <- midpoint
        node$log_upper <- midpoint
      }
      node$evaluated <- TRUE
      node$evidence_fit <- fit$fit
      node$evidence_fit$method <- if (
          is.null(node$evidence_fit$method)) {
        "stochastic"
      } else {
        node$evidence_fit$method
      }
      node$evidence_preflight <- fit$preflight
      node$interval_conflict <- FALSE
      nodes[[node_id]] <- node
    }
    actions <- actions + 1L
  }

  final_state <- if (identical(selection_mode, "penalized")) {
    rrs_tree_search_state(nodes, frontier_ids, lambda)
  } else {
    rrs_tree_epsilon_search_state(nodes, frontier_ids, epsilon)
  }
  budget_overshoot <- if (is.finite(max_score_evaluations)) {
    max(0, score_evaluations - max_score_evaluations)
  } else {
    0
  }
  epsilon_feasible <- identical(selection_mode, "epsilon") &&
    isTRUE(final_state$incumbent$feasible) &&
    final_state$omitted_mass_bound <= epsilon
  cost_optimality_met <- isTRUE(epsilon_feasible) &&
    is.finite(final_state$cost_gap) &&
    final_state$cost_gap <= cost_tolerance
  list(
    status = status,
    registered_status = rrs_registered_status(status),
    certified = identical(status, "certified"),
    bounds_met = identical(status, "certified"),
    epsilon_feasible = isTRUE(epsilon_feasible),
    cost_optimality_met = isTRUE(cost_optimality_met),
    selection_mode = selection_mode,
    epsilon = if (identical(selection_mode, "epsilon")) epsilon else NA_real_,
    selected_ids = final_state$selected_ids,
    objective_gap = if (identical(selection_mode, "penalized")) {
      final_state$objective_gap
    } else {
      NA_real_
    },
    cost_lower = if (identical(selection_mode, "epsilon")) {
      final_state$cost_lower
    } else {
      NA_real_
    },
    cost_upper = if (identical(selection_mode, "epsilon")) {
      final_state$cost_upper
    } else {
      NA_real_
    },
    cost_gap = if (identical(selection_mode, "epsilon")) {
      final_state$cost_gap
    } else {
      NA_real_
    },
    bracket_consistent = if (identical(selection_mode, "epsilon")) {
      final_state$bracket_consistent
    } else {
      NA
    },
    cost_tolerance = if (identical(selection_mode, "epsilon")) {
      cost_tolerance
    } else {
      NA_real_
    },
    omitted_mass_bound = final_state$omitted_mass_bound,
    frontier_ids = frontier_ids,
    nodes = nodes,
    trace = do.call(rbind, trace),
    total_score_evaluations = score_evaluations,
    bound_score_evaluations = bound_score_evaluations,
    evidence_score_evaluations = evidence_score_evaluations,
    score_evaluation_budget = max_score_evaluations,
    budget_overshoot = budget_overshoot,
    budget_checked_between_actions = TRUE,
    budget_checked_before_score_batches = TRUE,
    conflict_node_id = conflict_node_id,
    actions = actions,
    expansions = expansions,
    n_declared_nodes = tree$n_declared_nodes
  )
}

select_certified_region_tree <- certified_tree_region_search

select_minimum_cost_region_tree <- function(tree,
                                            log_score_fn,
                                            evidence_evaluator,
                                            epsilon,
                                            cost_tolerance = 0,
                                            ...) {
  certified_tree_region_search(
    tree = tree,
    log_score_fn = log_score_fn,
    evidence_evaluator = evidence_evaluator,
    selection_mode = "epsilon",
    epsilon = epsilon,
    cost_tolerance = cost_tolerance,
    ...
  )
}

fit_certified_region_selection <- function(tree,
                                           log_score_fn,
                                           lambda = NULL,
                                           delta = 0.05,
                                           n_pilot = 250L,
                                           n_paths = 1000L,
                                           ess_fraction = 0.8,
                                           resample_fraction = 0.5,
                                           pilot_mutation_steps = 1L,
                                           production_mutation_steps = 2L,
                                           maximum_relative_radius = NULL,
                                           maximum_paths = NULL,
                                           evidence_method = c("ais", "iid"),
                                           objective_tolerance = 0.05,
                                           omitted_mass_tolerance = 0.05,
                                           max_score_evaluations = Inf,
                                           max_actions = Inf,
                                           selection_mode = c(
                                             "penalized",
                                             "epsilon"
                                           ),
                                           epsilon = NULL,
                                           cost_tolerance = 0) {
  selection_mode <- match.arg(selection_mode)
  evidence_method <- match.arg(evidence_method)
  evaluator <- if (evidence_method == "iid") {
    make_iid_bounded_evidence_evaluator(
      log_score_fn = log_score_fn,
      n_declared_nodes = tree$n_declared_nodes,
      delta = delta,
      n_paths = n_paths,
      maximum_relative_radius = maximum_relative_radius,
      maximum_paths = maximum_paths
    )
  } else {
    make_pilot_fixed_evidence_evaluator(
      log_score_fn = log_score_fn,
      n_declared_nodes = tree$n_declared_nodes,
      delta = delta,
      n_pilot = n_pilot,
      n_paths = n_paths,
      ess_fraction = ess_fraction,
      resample_fraction = resample_fraction,
      pilot_mutation_steps = pilot_mutation_steps,
      production_mutation_steps = production_mutation_steps,
      maximum_relative_radius = maximum_relative_radius,
      maximum_paths = maximum_paths
    )
  }
  select_certified_region_tree(
    tree = tree,
    log_score_fn = log_score_fn,
    evidence_evaluator = evaluator,
    lambda = lambda,
    objective_tolerance = objective_tolerance,
    omitted_mass_tolerance = omitted_mass_tolerance,
    max_score_evaluations = max_score_evaluations,
    max_actions = max_actions,
    selection_mode = selection_mode,
    epsilon = epsilon,
    cost_tolerance = cost_tolerance
  )
}

approximate_selected_region_tree <- function(search,
                                             log_score_fn,
                                             n_local_particles = 1000L,
                                             ess_fraction = 0.8,
                                             resample_fraction = 0.5,
                                             local_mutation_steps = 2L) {
  if (!is.list(search) ||
      is.null(search$nodes) ||
      is.null(search$selected_ids) ||
      !is.function(log_score_fn)) {
    stop("Invalid certified selection or score function.", call. = FALSE)
  }
  selected_nodes <- if (isTRUE(search$certified)) {
    search$nodes[search$selected_ids]
  } else {
    list()
  }
  local_fits <- lapply(selected_nodes, function(node) {
    phi_schedule <- if (!is.null(node$evidence_fit$phi_schedule)) {
      node$evidence_fit$phi_schedule
    } else {
      c(0, 1)
    }
    region_annealed_smc(
      region = node$region,
      log_score_fn = log_score_fn,
      n_particles = n_local_particles,
      ess_fraction = ess_fraction,
      resample_fraction = resample_fraction,
      mutation_steps = local_mutation_steps,
      phi_schedule = phi_schedule
    )
  })
  local_score_evaluations <- if (length(local_fits) > 0L) {
    sum(vapply(local_fits, `[[`, numeric(1), "score_evaluations"))
  } else {
    0
  }
  regional_weights <- if (length(selected_nodes) > 0L) {
    selected_log_evidence <- vapply(
      seq_along(selected_nodes),
      function(index) {
        recorded <- selected_nodes[[index]]$evidence_fit$log_evidence
        if (length(recorded) == 1L && is.finite(recorded)) {
          recorded
        } else {
          local_fits[[index]]$log_evidence
        }
      },
      numeric(1)
    )
    exp(selected_log_evidence - rrs_log_sum_exp(selected_log_evidence))
  } else {
    numeric(0)
  }
  list(
    selected_ids = if (isTRUE(search$certified)) {
      search$selected_ids
    } else {
      character(0)
    },
    candidate_selected_ids = search$selected_ids,
    regional_weights = regional_weights,
    local_fits = local_fits,
    local_score_evaluations = local_score_evaluations
  )
}

fit_certified_region_tree <- function(tree,
                                      log_score_fn,
                                      lambda = NULL,
                                      delta = 0.05,
                                      n_pilot = 250L,
                                      n_paths = 1000L,
                                      n_local_particles = 1000L,
                                      ess_fraction = 0.8,
                                      resample_fraction = 0.5,
                                      pilot_mutation_steps = 1L,
                                      production_mutation_steps = 2L,
                                      maximum_relative_radius = NULL,
                                      maximum_paths = NULL,
                                      evidence_method = c("ais", "iid"),
                                      local_mutation_steps = 2L,
                                      objective_tolerance = 0.05,
                                      omitted_mass_tolerance = 0.05,
                                      max_score_evaluations = Inf,
                                      max_actions = Inf,
                                      selection_mode = c(
                                        "penalized",
                                        "epsilon"
                                      ),
                                      epsilon = NULL,
                                      cost_tolerance = 0) {
  selection_mode <- match.arg(selection_mode)
  evidence_method <- match.arg(evidence_method)
  search <- fit_certified_region_selection(
    tree = tree,
    log_score_fn = log_score_fn,
    lambda = lambda,
    delta = delta,
    n_pilot = n_pilot,
    n_paths = n_paths,
    ess_fraction = ess_fraction,
    resample_fraction = resample_fraction,
    pilot_mutation_steps = pilot_mutation_steps,
    production_mutation_steps = production_mutation_steps,
    maximum_relative_radius = maximum_relative_radius,
    maximum_paths = maximum_paths,
    evidence_method = evidence_method,
    objective_tolerance = objective_tolerance,
    omitted_mass_tolerance = omitted_mass_tolerance,
    max_score_evaluations = max_score_evaluations,
    max_actions = max_actions,
    selection_mode = selection_mode,
    epsilon = epsilon,
    cost_tolerance = cost_tolerance
  )
  downstream <- approximate_selected_region_tree(
    search = search,
    log_score_fn = log_score_fn,
    n_local_particles = n_local_particles,
    ess_fraction = ess_fraction,
    resample_fraction = resample_fraction,
    local_mutation_steps = local_mutation_steps
  )
  list(
    search = search,
    status = search$status,
    registered_status = search$registered_status,
    selected_ids = downstream$selected_ids,
    candidate_selected_ids = downstream$candidate_selected_ids,
    regional_weights = downstream$regional_weights,
    local_fits = downstream$local_fits,
    omitted_mass_bound = search$omitted_mass_bound,
    objective_gap = search$objective_gap,
    cost_lower = search$cost_lower,
    cost_upper = search$cost_upper,
    cost_gap = search$cost_gap,
    bracket_consistent = search$bracket_consistent,
    bounds_met = search$bounds_met,
    selection_mode = search$selection_mode,
    epsilon = search$epsilon,
    certified = search$certified,
    budget_overshoot = search$budget_overshoot,
    conflict_node_id = search$conflict_node_id,
    bound_score_evaluations = search$bound_score_evaluations,
    evidence_score_evaluations = search$evidence_score_evaluations,
    local_score_evaluations = downstream$local_score_evaluations,
    total_score_evaluations =
      search$total_score_evaluations + downstream$local_score_evaluations
  )
}

rrs_poisson_binomial_pmf <- function(probabilities) {
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) < 1L ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1)) {
    stop("Probabilities must lie strictly between zero and one.", call. = FALSE)
  }
  pmf <- c(1, rep(0, length(probabilities)))
  used <- 0L
  for (probability in probabilities) {
    used <- used + 1L
    previous <- pmf
    pmf[1L] <- previous[1L] * (1 - probability)
    for (k in seq_len(used)) {
      pmf[k + 1L] <-
        previous[k + 1L] * (1 - probability) +
        previous[k] * probability
    }
  }
  pmf
}

rrs_sample_bernoulli_given_count <- function(probabilities, count, n) {
  probabilities <- as.numeric(probabilities)
  count <- as.integer(count)
  n <- as.integer(n)
  m <- length(probabilities)
  if (count < 0L || count > m || n < 1L) {
    stop("Invalid conditional Bernoulli sample size or count.", call. = FALSE)
  }
  if (count == 0L) {
    return(matrix(0L, nrow = n, ncol = m))
  }
  odds <- probabilities / (1 - probabilities)
  suffix <- matrix(0, nrow = m + 1L, ncol = count + 1L)
  suffix[m + 1L, 1L] <- 1
  for (i in seq.int(m, 1L)) {
    suffix[i, 1L] <- 1
    for (k in seq_len(min(count, m - i + 1L))) {
      suffix[i, k + 1L] <-
        suffix[i + 1L, k + 1L] +
        odds[i] * suffix[i + 1L, k]
    }
  }
  out <- matrix(0L, nrow = n, ncol = m)
  for (row in seq_len(n)) {
    remaining <- count
    for (i in seq_len(m)) {
      if (remaining == 0L) {
        break
      }
      positions_left <- m - i + 1L
      include_probability <- if (remaining == positions_left) {
        1
      } else {
        odds[i] * suffix[i + 1L, remaining] /
          suffix[i, remaining + 1L]
      }
      if (runif(1) <= include_probability) {
        out[row, i] <- 1L
        remaining <- remaining - 1L
      }
    }
  }
  out
}

new_group_state_region <- function(inclusion_probabilities,
                                   groups,
                                   states,
                                   cost = 1L,
                                   label = NULL) {
  inclusion_probabilities <- as.numeric(inclusion_probabilities)
  p <- length(inclusion_probabilities)
  if (p < 1L ||
      any(!is.finite(inclusion_probabilities)) ||
      any(inclusion_probabilities <= 0 | inclusion_probabilities >= 1)) {
    stop("Inclusion probabilities must lie in (0, 1).", call. = FALSE)
  }
  groups <- lapply(groups, as.integer)
  if (!identical(sort(unname(unlist(groups))), seq_len(p))) {
    stop("groups must form a partition of the predictor indices.", call. = FALSE)
  }
  states <- as.character(states)
  allowed <- c("free", "zero", "one", "multi")
  if (length(states) != length(groups) || any(!states %in% allowed)) {
    stop("Each group needs a free, zero, one, or multi state.", call. = FALSE)
  }

  event_probabilities <- vapply(seq_along(groups), function(k) {
    probs <- inclusion_probabilities[groups[[k]]]
    pmf <- rrs_poisson_binomial_pmf(probs)
    switch(
      states[k],
      free = 1,
      zero = pmf[1L],
      one = if (length(pmf) >= 2L) pmf[2L] else 0,
      multi = if (length(pmf) >= 3L) sum(pmf[3:length(pmf)]) else 0
    )
  }, numeric(1))
  if (any(event_probabilities <= 0)) {
    stop("A declared group state has zero prior probability.", call. = FALSE)
  }

  sample_group_state <- function(n, probabilities, state) {
    n <- as.integer(n)
    m <- length(probabilities)
    if (state == "free") {
      return(matrix(
        rbinom(n * m, 1L, rep(probabilities, each = n)),
        nrow = n
      ))
    }
    if (state == "zero") {
      return(matrix(0L, nrow = n, ncol = m))
    }
    if (state == "one") {
      out <- matrix(0L, nrow = n, ncol = m)
      odds <- probabilities / (1 - probabilities)
      selected <- sample.int(m, size = n, replace = TRUE, prob = odds)
      out[cbind(seq_len(n), selected)] <- 1L
      return(out)
    }
    pmf <- rrs_poisson_binomial_pmf(probabilities)
    counts <- sample.int(
      m - 1L,
      size = n,
      replace = TRUE,
      prob = pmf[3:length(pmf)]
    ) + 1L
    out <- matrix(0L, nrow = n, ncol = m)
    for (count in unique(counts)) {
      rows <- which(counts == count)
      out[rows, ] <- rrs_sample_bernoulli_given_count(
        probabilities,
        count,
        length(rows)
      )
    }
    out
  }

  sample_region <- function(n) {
    n <- as.integer(n)
    out <- matrix(0L, nrow = n, ncol = p)
    for (k in seq_along(groups)) {
      indices <- groups[[k]]
      out[, indices] <- sample_group_state(
        n,
        inclusion_probabilities[indices],
        states[k]
      )
    }
    out
  }
  mutable_groups <- which(vapply(seq_along(groups), function(k) {
    state <- states[k]
    size <- length(groups[[k]])
    state == "free" ||
      (state == "one" && size > 1L) ||
      (state == "multi" && size > 2L)
  }, logical(1)))
  mutate_region <- function(particles) {
    particles <- as.matrix(particles)
    n <- nrow(particles)
    if (length(mutable_groups) == 0L) {
      return(list(
        particles = particles,
        proposed = rep(FALSE, n),
        type = "group-refresh"
      ))
    }
    selected_groups <- sample(mutable_groups, n, replace = TRUE)
    out <- particles
    for (k in unique(selected_groups)) {
      rows <- which(selected_groups == k)
      indices <- groups[[k]]
      out[rows, indices] <- sample_group_state(
        length(rows),
        inclusion_probabilities[indices],
        states[k]
      )
    }
    list(
      particles = out,
      proposed = rep(TRUE, n),
      type = "group-refresh"
    )
  }
  log_prob_inside <- function(gamma) {
    gamma <- as.matrix(gamma)
    rowSums(
      sweep(gamma, 2L, log(inclusion_probabilities), `*`) +
        sweep(1 - gamma, 2L, log1p(-inclusion_probabilities), `*`)
    ) - sum(log(event_probabilities))
  }
  log_prob <- function(gamma) {
    gamma <- as.matrix(gamma)
    if (ncol(gamma) != p) {
      stop("Support dimension does not match the region.", call. = FALSE)
    }
    valid <- rowSums(!(gamma == 0 | gamma == 1)) == 0L
    for (k in seq_along(groups)) {
      count <- rowSums(gamma[, groups[[k]], drop = FALSE])
      valid <- valid & switch(
        states[k],
        free = TRUE,
        zero = count == 0L,
        one = count == 1L,
        multi = count >= 2L
      )
    }
    out <- rep(-Inf, nrow(gamma))
    if (any(valid)) {
      out[valid] <- log_prob_inside(gamma[valid, , drop = FALSE])
    }
    out
  }
  anchor <- function() {
    out <- matrix(0L, nrow = 1L, ncol = p)
    for (k in seq_along(groups)) {
      indices <- groups[[k]]
      probs <- inclusion_probabilities[indices]
      state <- states[k]
      if (state == "free") {
        out[1L, indices] <- as.integer(probs >= 0.5)
      } else if (state == "one") {
        out[1L, indices[which.max(probs / (1 - probs))]] <- 1L
      } else if (state == "multi") {
        chosen <- which(probs >= 0.5)
        if (length(chosen) < 2L) {
          chosen <- order(probs, decreasing = TRUE)[1:2]
        }
        out[1L, indices[chosen]] <- 1L
      }
    }
    out
  }

  structure(
    list(
      sample = sample_region,
      mutate = mutate_region,
      log_prob = log_prob,
      log_prob_inside = log_prob_inside,
      anchor = anchor,
      cost = as.integer(cost),
      label = if (is.null(label)) {
        paste(states, collapse = "/")
      } else {
        as.character(label)
      },
      states = states,
      groups = groups,
      inclusion_probabilities = inclusion_probabilities,
      log_event_probability = sum(log(event_probabilities))
    ),
    class = "rrs_region"
  )
}

make_group_gibbs_target_mutation <- function(states,
                                             region,
                                             groups_per_sweep = NULL) {
  states <- as.character(states)
  groups <- region$groups
  p <- length(region$inclusion_probabilities)
  local_supports <- lapply(seq_along(groups), function(group) {
    size <- length(groups[[group]])
    supports <- as.matrix(expand.grid(rep(list(c(0L, 1L)), size)))
    count <- rowSums(supports)
    selected <- switch(
      states[group],
      free = seq_len(nrow(supports)),
      zero = which(count == 0L),
      one = which(count == 1L),
      multi = which(count >= 2L)
    )
    supports[selected, , drop = FALSE]
  })
  mutable <- which(vapply(local_supports, nrow, integer(1)) > 1L)
  if (is.null(groups_per_sweep)) {
    groups_per_sweep <- length(mutable)
  }
  groups_per_sweep <- as.integer(groups_per_sweep)
  if (length(mutable) > 0L && groups_per_sweep < 1L) {
    stop("groups_per_sweep must be positive.", call. = FALSE)
  }

  function(particles, phi, log_score_fn) {
    particles <- as.matrix(particles)
    n <- nrow(particles)
    if (ncol(particles) != p) {
      stop("Group Gibbs mutation received the wrong support dimension.", call. = FALSE)
    }
    if (length(mutable) == 0L) {
      log_nu <- region$log_prob(particles)
      log_score <- rrs_evaluate_log_score(log_score_fn, particles)
      return(list(
        particles = particles,
        log_g = log_score - log_nu,
        score_evaluations = n,
        accepted = n,
        proposed = n,
        type = "group-gibbs"
      ))
    }

    out <- particles
    final_log_score <- NULL
    final_log_nu <- NULL
    score_evaluations <- 0L
    sweep_groups <- if (groups_per_sweep >= length(mutable)) {
      mutable
    } else {
      sample(mutable, groups_per_sweep, replace = FALSE)
    }
    for (group in sweep_groups) {
      configurations <- local_supports[[group]]
      n_configurations <- nrow(configurations)
      candidate <- out[
        rep(seq_len(n), each = n_configurations),
        ,
        drop = FALSE
      ]
      candidate[, groups[[group]]] <- configurations[
        rep(seq_len(n_configurations), times = n),
        ,
        drop = FALSE
      ]
      candidate_log_nu <- region$log_prob_inside(candidate)
      candidate_log_score <- rrs_evaluate_log_score(
        log_score_fn,
        candidate
      )
      score_evaluations <- score_evaluations + nrow(candidate)
      log_target <- matrix(
        (1 - phi) * candidate_log_nu + phi * candidate_log_score,
        nrow = n,
        ncol = n_configurations,
        byrow = TRUE
      )
      selected <- integer(n)
      for (row in seq_len(n)) {
        probabilities <- exp(
          log_target[row, ] - rrs_log_sum_exp(log_target[row, ])
        )
        selected[row] <- sample.int(
          n_configurations,
          size = 1L,
          prob = probabilities
        )
      }
      selected_rows <- (seq_len(n) - 1L) * n_configurations + selected
      out <- candidate[selected_rows, , drop = FALSE]
      final_log_score <- candidate_log_score[selected_rows]
      final_log_nu <- candidate_log_nu[selected_rows]
    }
    list(
      particles = out,
      log_g = final_log_score - final_log_nu,
      score_evaluations = score_evaluations,
      accepted = score_evaluations,
      proposed = score_evaluations,
      type = "group-gibbs"
    )
  }
}

make_additive_group_state_rho_bound <- function(log_lower, log_upper) {
  required_states <- c("free", "zero", "one", "multi")
  log_lower <- as.matrix(log_lower)
  log_upper <- as.matrix(log_upper)
  if (!identical(dim(log_lower), dim(log_upper)) ||
      nrow(log_lower) < 1L ||
      !all(required_states %in% colnames(log_lower)) ||
      !all(required_states %in% colnames(log_upper))) {
    stop(
      "Lower and upper matrices need matching rows and named group-state columns.",
      call. = FALSE
    )
  }
  log_lower <- log_lower[, required_states, drop = FALSE]
  log_upper <- log_upper[, required_states, drop = FALSE]
  if (any(!is.finite(log_lower)) ||
      any(!is.finite(log_upper)) ||
      any(log_lower > log_upper)) {
    stop("Every declared group-state log interval must be finite and ordered.", call. = FALSE)
  }
  widths <- log_upper - log_lower
  function(states, depth = NULL) {
    states <- as.character(states)
    if (length(states) != nrow(widths) ||
        any(!states %in% required_states)) {
      stop("The state vector does not match the declared groups.", call. = FALSE)
    }
    sum(widths[cbind(seq_along(states), match(states, required_states))])
  }
}

new_group_state_tree <- function(inclusion_probabilities,
                                 groups,
                                 split_order = seq_along(groups),
                                 max_depth = length(split_order),
                                 rho = 0,
                                 log_g_bounds = NULL,
                                 log_evidence_bounds = NULL,
                                 target_mutation = NULL,
                                 region_factory = NULL,
                                 state_cost = c(
                                   zero = 1L,
                                   one = 2L,
                                   multi = 2L
                                 ),
                                 base_cost = 1L,
                                 materialize = FALSE) {
  split_order <- as.integer(split_order)
  max_depth <- as.integer(max_depth)
  if (max_depth < 0L || max_depth > length(split_order) ||
      anyDuplicated(split_order) ||
      any(!split_order %in% seq_along(groups))) {
    stop("Invalid group split order or depth.", call. = FALSE)
  }
  if (length(rho) != 1L && !is.function(rho)) {
    stop("rho must be one number or a node-specific function.", call. = FALSE)
  }
  if (!is.null(log_g_bounds) && !is.function(log_g_bounds)) {
    stop("log_g_bounds must be NULL or a node-specific function.", call. = FALSE)
  }
  if (!is.null(log_evidence_bounds) &&
      !is.function(log_evidence_bounds)) {
    stop(
      "log_evidence_bounds must be NULL or a node-specific function.",
      call. = FALSE
    )
  }
  if (!is.null(target_mutation) && !is.function(target_mutation)) {
    stop("target_mutation must be NULL or a region factory.", call. = FALSE)
  }
  if (!is.null(region_factory) && !is.function(region_factory)) {
    stop("region_factory must be NULL or a function.", call. = FALSE)
  }
  state_cost <- state_cost[c("zero", "one", "multi")]
  state_cost <- stats::setNames(
    as.integer(state_cost),
    c("zero", "one", "multi")
  )
  if (any(!is.finite(state_cost)) || any(state_cost < 1L)) {
    stop("Every group-state cost must be a positive integer.", call. = FALSE)
  }
  if (length(base_cost) != 1L ||
      !is.finite(base_cost) ||
      base_cost < 1L ||
      abs(base_cost - round(base_cost)) > sqrt(.Machine$double.eps)) {
    stop("base_cost must be a positive integer.", call. = FALSE)
  }
  base_cost <- as.integer(round(base_cost))

  node_from_states <- function(states, depth, id) {
    terminal <- depth >= max_depth
    remaining <- if (terminal) 0L else max_depth - depth
    realized_cost <- base_cost + sum(
      state_cost[states[states != "free"]]
    )
    cost_lower <- realized_cost + remaining * min(state_cost)
    node_rho <- if (is.function(rho)) rho(states, depth) else rho
    region_cost <- if (terminal) realized_cost else cost_lower
    region_label <- paste0("depth-", depth, ":", paste(states, collapse = "/"))
    region <- if (is.null(region_factory)) {
      new_group_state_region(
        inclusion_probabilities = inclusion_probabilities,
        groups = groups,
        states = states,
        cost = region_cost,
        label = region_label
      )
    } else {
      region_factory(states, region_cost, region_label)
    }
    rrs_validate_region(region)
    if (is.function(log_g_bounds)) {
      node_states <- states
      node_region <- region
      region$log_g_bounds_factory <- function(candidate_region) {
        log_g_bounds(node_states, candidate_region)
      }
      region$log_g_bounds <- function() region$log_g_bounds_factory(node_region)
    }
    if (is.function(log_evidence_bounds)) {
      node_states <- states
      node_region <- region
      region$log_evidence_bounds <- function() {
        log_evidence_bounds(node_states, node_region)
      }
    }
    if (is.function(target_mutation)) {
      region$target_mutate <- target_mutation(states, region)
    }
    new_region_tree_node(
      id = id,
      region = region,
      children = character(0),
      rho = node_rho,
      cost_lower = cost_lower,
      terminal = terminal
    )
  }

  expand_node <- function(node) {
    states <- node$region$states
    depth <- sum(states[split_order] != "free")
    if (depth >= max_depth) {
      stop("A terminal group-state node cannot be expanded.", call. = FALSE)
    }
    group_index <- split_order[depth + 1L]
    available_states <- if (length(groups[[group_index]]) >= 2L) {
      c("zero", "one", "multi")
    } else {
      c("zero", "one")
    }
    lapply(available_states, function(state) {
      child_states <- states
      child_states[group_index] <- state
      child_id <- paste0(node$id, "/g", group_index, "-", state)
      node_from_states(child_states, depth + 1L, child_id)
    })
  }

  branch_counts <- if (max_depth > 0L) {
    vapply(split_order[seq_len(max_depth)], function(group_index) {
      if (length(groups[[group_index]]) >= 2L) 3 else 2
    }, numeric(1))
  } else {
    numeric(0)
  }
  n_declared_nodes <- 1
  if (length(branch_counts) > 0L) {
    generation_size <- 1
    for (branch_count in branch_counts) {
      generation_size <- generation_size * branch_count
      n_declared_nodes <- n_declared_nodes + generation_size
    }
  }

  root <- node_from_states(
    rep("free", length(groups)),
    depth = 0L,
    id = "root"
  )
  lazy_tree <- new_lazy_region_tree(
    root_node = root,
    n_declared_nodes = n_declared_nodes,
    expand_node = expand_node
  )
  if (!isTRUE(materialize)) {
    return(lazy_tree)
  }

  nodes <- lazy_tree$nodes
  queue <- lazy_tree$root_id
  while (length(queue) > 0L) {
    node_id <- queue[1L]
    queue <- queue[-1L]
    node <- nodes[[node_id]]
    if (node$terminal) {
      next
    }
    generated <- expand_node(node)
    node$children <- vapply(generated, `[[`, character(1), "id")
    nodes[[node_id]] <- node
    for (child in generated) {
      nodes[[child$id]] <- child
    }
    queue <- c(queue, node$children)
  }
  out <- new_region_tree(nodes, lazy_tree$root_id)
  if (out$n_declared_nodes != n_declared_nodes) {
    stop("Materialized tree size does not match its declared size.", call. = FALSE)
  }
  out
}

gaussian_swap_log_score_bound <- function(column_distance,
                                          column_norm_bound,
                                          response_norm,
                                          sigma,
                                          slab_scale,
                                          log_prior_ratio_bound = 0) {
  inputs <- c(
    column_distance,
    column_norm_bound,
    response_norm,
    sigma,
    slab_scale,
    log_prior_ratio_bound
  )
  if (any(!is.finite(inputs)) ||
      any(inputs[1:3] < 0) ||
      sigma <= 0 ||
      slab_scale < 0 ||
      log_prior_ratio_bound < 0) {
    stop("Invalid Gaussian swap-bound inputs.", call. = FALSE)
  }
  tau2 <- slab_scale^2
  coefficient <-
    tau2 * column_norm_bound +
    tau2 * response_norm^2 * column_norm_bound / sigma^2 +
    tau2^2 * response_norm^2 * column_norm_bound^3 / sigma^2
  log_prior_ratio_bound + coefficient * column_distance
}
