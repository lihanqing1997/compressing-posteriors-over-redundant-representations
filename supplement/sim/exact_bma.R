log_sum_exp <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(-Inf)
  }
  m <- max(x[finite])
  m + log(sum(exp(x[finite] - m)))
}

log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  both_neg_inf <- !is.finite(m)
  out <- m + log(exp(a - m) + exp(b - m))
  out[both_neg_inf] <- -Inf
  out
}

safe_chol <- function(A, jitter = 1e-9, max_tries = 8L) {
  A <- as.matrix(A)
  A[!is.finite(A)] <- 0
  A <- (A + t(A)) / 2
  q <- nrow(A)
  for (k in seq_len(max_tries)) {
    ridge <- jitter * 10^(k - 1L)
    out <- tryCatch(chol(A + diag(ridge, q)), error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  eig <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  ridge <- max(jitter, -min(eig, na.rm = TRUE) + jitter)
  chol(A + diag(ridge, q))
}

enumerate_supports <- function(p) {
  ids <- 0:(2^p - 1)
  supports <- sapply(seq_len(p), function(j) {
    as.integer(bitwAnd(ids, bitwShiftL(1L, j - 1L)) != 0L)
  })
  colnames(supports) <- paste0("x", seq_len(p))
  supports
}

log_marginal_support <- function(gamma, X, y, tau2 = 4, a0 = 1, b0 = 1) {
  n <- nrow(X)
  active <- which(gamma == 1)
  q <- length(active)
  a_n <- a0 + n / 2
  yty <- sum(y^2)

  if (q == 0) {
    log_det <- 0
    quad <- yty
  } else {
    Xg <- X[, active, drop = FALSE]
    S <- diag(q) + tau2 * crossprod(Xg)
    R <- safe_chol(S)
    log_det <- 2 * sum(log(diag(R)))
    Xty <- crossprod(Xg, y)
    sol <- backsolve(R, forwardsolve(t(R), Xty))
    quad <- yty - tau2 * sum(Xty * sol)
    quad <- max(quad, 0)
  }

  b_gamma <- b0 + 0.5 * quad

  lgamma(a_n) - lgamma(a0) -
    (n / 2) * log(2 * pi) +
    a0 * log(b0) -
    0.5 * log_det -
    a_n * log(b_gamma)
}

fit_exact_bma <- function(X,
                          y,
                          theta = 0.08,
                          tau2 = 4,
                          a0 = 1,
                          b0 = 1,
                          supports = NULL) {
  p <- ncol(X)
  if (is.null(supports)) {
    supports <- enumerate_supports(p)
  }

  log_marginal <- numeric(nrow(supports))
  for (i in seq_len(nrow(supports))) {
    log_marginal[i] <- log_marginal_support(
      gamma = supports[i, ],
      X = X,
      y = y,
      tau2 = tau2,
      a0 = a0,
      b0 = b0
    )
  }
  size <- rowSums(supports)
  log_prior <- size * log(theta) + (p - size) * log1p(-theta)
  log_unnormalized <- log_marginal + log_prior
  log_evidence <- log_sum_exp(log_unnormalized)
  posterior <- exp(log_unnormalized - log_evidence)

  list(
    supports = supports,
    log_marginal = log_marginal,
    log_prior = log_prior,
    log_posterior = log(posterior),
    posterior = posterior,
    log_evidence = log_evidence,
    theta = theta,
    tau2 = tau2,
    a0 = a0,
    b0 = b0,
    n = nrow(X),
    p = p
  )
}

posterior_summary <- function(weights, mass = 0.95) {
  if (length(weights) == 0L || !is.finite(sum(weights)) || sum(weights) <= 0) {
    return(list(
      entropy = NA_real_,
      effective_count = NA_real_,
      n_mass = NA_integer_,
      max_probability = NA_real_
    ))
  }
  weights <- weights / sum(weights)
  positive <- weights[weights > 0]
  entropy <- -sum(positive * log(positive))
  sorted <- sort(weights, decreasing = TRUE)
  n_mass <- which(cumsum(sorted) >= mass)[1]

  list(
    entropy = entropy,
    effective_count = exp(entropy),
    n_mass = n_mass,
    max_probability = max(weights)
  )
}

predict_mixture <- function(fit,
                            X_train,
                            y_train,
                            X_test,
                            y_test = NULL,
                            weights = fit$posterior,
                            weight_tol = 1e-12) {
  supports <- fit$supports
  tau2 <- fit$tau2
  a0 <- fit$a0
  b0 <- fit$b0
  a_n <- a0 + nrow(X_train) / 2
  yty <- sum(y_train^2)
  n_test <- nrow(X_test)

  weights <- weights / sum(weights)
  keep <- which(weights > weight_tol)
  mean_mix <- rep(0, n_test)
  log_density_mix <- rep(-Inf, n_test)

  for (idx in keep) {
    gamma <- supports[idx, ]
    active <- which(gamma == 1)
    q <- length(active)

    if (q == 0) {
      mean_model <- rep(0, n_test)
      b_gamma <- b0 + 0.5 * yty
      scale <- rep(sqrt(b_gamma / a_n), n_test)
    } else {
      Xg <- X_train[, active, drop = FALSE]
      Xt <- X_test[, active, drop = FALSE]
      precision <- diag(1 / tau2, q) + crossprod(Xg)
      R <- safe_chol(precision)
      Sigma <- chol2inv(R)
      mu <- Sigma %*% crossprod(Xg, y_train)
      mean_model <- as.numeric(Xt %*% mu)
      quad <- yty - sum(crossprod(Xg, y_train) * mu)
      quad <- max(quad, 0)
      b_gamma <- b0 + 0.5 * quad
      leverage <- rowSums((Xt %*% Sigma) * Xt)
      scale <- sqrt((b_gamma / a_n) * (1 + leverage))
    }

    w <- weights[idx]
    mean_mix <- mean_mix + w * mean_model

    if (!is.null(y_test)) {
      log_density <- dt((y_test - mean_model) / scale,
                        df = 2 * a_n,
                        log = TRUE) - log(scale)
      log_density_mix <- log_add_exp(log_density_mix, log(w) + log_density)
    }
  }

  out <- list(mean = mean_mix)
  if (!is.null(y_test)) {
    out$rmse <- sqrt(mean((y_test - mean_mix)^2))
    out$mean_log_score <- mean(log_density_mix)
  }
  out
}
