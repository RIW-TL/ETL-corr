## ADT estimator for bandable covariance matrices.

options(stringsAsFactors = FALSE, warn = 1)

spectral_norm <- function(A) {
  if (!length(A) || all(A == 0)) return(0)
  svd(A, nu = 0L, nv = 0L)$d[1L]
}

sample_cov_n <- function(X) {
  Xc <- sweep(X, 2L, colMeans(X), FUN = "-")
  crossprod(Xc) / nrow(Xc)
}

block_ids <- function(p, bandwidth) {
  block_size <- max(1L, as.integer(floor(bandwidth)))
  1L + ((seq_len(p) - 1L) %/% block_size)
}

block_tridiagonal <- function(A, bandwidth) {
  ids <- block_ids(nrow(A), bandwidth)
  keep <- abs(outer(ids, ids, FUN = "-")) <= 1L
  A * keep
}

superband <- function(A, lower, upper) {
  p <- nrow(A)
  gap <- outer(
    seq_len(p), seq_len(p), FUN = function(j, l) l - j
  )
  A * (gap > lower & gap <= upper)
}

number_of_strips <- function(bandwidth, n_eff) {
  cap <- n_eff / log(n_eff)
  if (bandwidth > cap) return(0L)
  as.integer(floor(1 + log(cap / bandwidth, base = 2)))
}

strip_bias_proxy <- function(A, bandwidth, n_eff) {
  ell_max <- number_of_strips(bandwidth, n_eff)
  if (ell_max <= 0L) return(0)
  
  total <- 0
  for (ell in seq_len(ell_max)) {
    lower <- 2^(ell - 1L) * bandwidth
    upper <- 2^ell * bandwidth
    total <- total + spectral_norm(superband(A, lower, upper))
  }
  total
}

bandwidth_grid <- function(n_eff, p) {
  lower <- log(p)
  upper <- max(min(n_eff / log(n_eff), p), lower)
  
  integer_points <- seq.int(
    max(1L, as.integer(ceiling(lower))),
    min(
      p,
      max(as.integer(ceiling(lower)), as.integer(floor(upper)))
    )
  )
  
  unique(c(lower, integer_points[integer_points <= upper]))
}

select_bandwidth <- function(A, n_eff, log_term, multiplier, rho) {
  grid <- bandwidth_grid(n_eff, nrow(A))
  objective <- vapply(grid, function(bandwidth) {
    strip_bias_proxy(A, bandwidth, n_eff) +
      multiplier * rho^2 * sqrt((bandwidth + log_term) / n_eff)
  }, numeric(1L))
  
  selected <- which.min(objective)
  list(
    bandwidth = grid[selected],
    objective = objective[selected],
    grid = grid,
    objective_path = objective
  )
}

block_threshold_tridiagonal <- function(A, bandwidth, threshold) {
  p <- nrow(A)
  ids <- block_ids(p, bandwidth)
  groups <- split(seq_len(p), ids)
  result <- matrix(0, p, p)
  
  for (row_block in seq_along(groups)) {
    last_column_block <- min(length(groups), row_block + 1L)
    
    for (column_block in row_block:last_column_block) {
      rows <- groups[[row_block]]
      columns <- groups[[column_block]]
      current_block <- A[rows, columns, drop = FALSE]
      
      if (spectral_norm(current_block) >= threshold) {
        result[rows, columns] <- current_block
        
        if (column_block != row_block) {
          result[columns, rows] <- t(current_block)
        }
      }
    }
  }
  result
}

ADT <- function(target_data, source_data) {
  ## The paper only specifies that tau, xi1 and xi2 should be sufficiently
  ## large constants. These fixed values are the finite-sample calibration
  ## used in the accompanying reproduction.
  tau <- 15
  xi1 <- 15
  xi2 <- 2
  rho <- 1
  
  target_data <- as.matrix(target_data)
  source_data <- as.matrix(source_data)
  
  if (ncol(target_data) != ncol(source_data)) {
    stop("Target and source dimensions do not agree.")
  }
  
  n0 <- nrow(target_data)
  n1 <- nrow(source_data)
  n_total <- n0 + n1
  p <- ncol(target_data)
  
  S0 <- sample_cov_n(target_data)
  S1 <- sample_cov_n(source_data)
  S_pool <- (n0 / n_total) * S0 + (n1 / n_total) * S1
  D_hat <- S0 - S_pool
  
  fit_v <- select_bandwidth(
    A = S_pool,
    n_eff = n_total,
    log_term = log(p),
    multiplier = tau,
    rho = rho
  )
  
  fit_w <- select_bandwidth(
    A = D_hat,
    n_eff = n0,
    log_term = log(n_total * p),
    multiplier = xi1,
    rho = rho
  )
  
  pooled_component <- block_tridiagonal(S_pool, fit_v$bandwidth)
  correction_used <- log(p) <= n0 / log(n0)
  threshold <- xi2 * rho^2 * sqrt(
    (fit_w$bandwidth + log(n_total * p)) / n0
  )
  
  correction_component <- if (correction_used) {
    block_threshold_tridiagonal(D_hat, fit_w$bandwidth, threshold)
  } else {
    matrix(0, p, p)
  }
  
  Sigma_hat <- pooled_component + correction_component
  Sigma_hat <- (Sigma_hat + t(Sigma_hat)) / 2
  
  list(
    estimate = Sigma_hat,
    v = fit_v$bandwidth,
    w = fit_w$bandwidth,
    threshold = threshold,
    correction_used = correction_used,
    objective_v = fit_v$objective,
    objective_w = fit_w$objective
  )
}

fit_band <- function(target_data, lambda = 15, rho = 1) {
  S0 <- sample_cov_n(target_data)
  fit <- select_bandwidth(
    S0, nrow(target_data), log(ncol(target_data)), lambda, rho
  )
  list(
    estimate = block_tridiagonal(S0, fit$bandwidth),
    bandwidth = fit$bandwidth,
    objective = fit$objective
  )
}
