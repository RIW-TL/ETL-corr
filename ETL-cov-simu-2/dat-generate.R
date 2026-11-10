library(mvtnorm)

## Data generator for bandable covariance matrices.

options(stringsAsFactors = FALSE, warn = 1)

matrix_l1_norm <- function(A) {
  max(colSums(abs(A)))
}

min_eigenvalue <- function(A) {
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values)
}

rmvnorm_chol <- function(n, Sigma) {
  R <- chol(Sigma)
  matrix(rnorm(n * nrow(Sigma)), nrow = n, ncol = nrow(Sigma)) %*% R
}

symmetric_rademacher <- function(p) {
  A <- matrix(0, p, p)
  A[upper.tri(A)] <- sample(
    c(-1, 1), p * (p - 1L) / 2L, replace = TRUE
  )
  A + t(A)
}

symmetric_uniform <- function(p) {
  ## For j < l, U_jl ~ Uniform(0,1), and U_lj = U_jl.
  ## Therefore I(U_jl <= eta) ~ Bernoulli(eta), and use of the same U makes
  ## the contaminated-entry sets nested as eta increases.
  U <- matrix(0, p, p)
  U[upper.tri(U)] <- runif(p * (p - 1L) / 2L)
  U <- U + t(U)
  diag(U) <- 1
  U
}

eta_key <- function(eta) {
  formatC(eta, digits = 8L, format = "fg", flag = "#")
}

make_population_path <- function(
    p, eta_values, delta_max = 1.1, max_attempts = 5000L,
    eigen_tolerance = 1e-8 ,C_t_value) {
  distance <- abs(outer(seq_len(p), seq_len(p), FUN = "-"))
  off_diagonal <- distance > 0
  
  for (attempt in seq_len(max_attempts)) {
    Xi_t <- symmetric_rademacher(p)
    Xi_s <- symmetric_rademacher(p)
    U <- symmetric_uniform(p)
    
    Sigma0 <- diag(p)
    Sigma0[off_diagonal] <-
      C_t_value * distance[off_diagonal]^(-3) * Xi_t[off_diagonal]
    
    D <- matrix(0, p, p)
    D[off_diagonal] <-
      distance[off_diagonal]^(-2) * Xi_s[off_diagonal]
    
    c_value <- delta_max / matrix_l1_norm(D)
    
    Sigma1_by_eta <- vector("list", length(eta_values))
    M_by_eta <- vector("list", length(eta_values))
    names(Sigma1_by_eta) <- vapply(eta_values, eta_key, character(1L))
    names(M_by_eta) <- names(Sigma1_by_eta)
    
    all_positive_definite <- min_eigenvalue(Sigma0) > eigen_tolerance
    
    for (i in seq_along(eta_values)) {
      M <- 1L * (U <= eta_values[i])
      diag(M) <- 0L
      
      Sigma1 <- Sigma0 + c_value * (M * D)
      Sigma1 <- (Sigma1 + t(Sigma1)) / 2
      
      Sigma1_by_eta[[i]] <- Sigma1
      M_by_eta[[i]] <- M
      
      if (min_eigenvalue(Sigma1) <= eigen_tolerance) {
        all_positive_definite <- FALSE
        break
      }
    }
    
    if (all_positive_definite) {
      delta_values <- vapply(
        Sigma1_by_eta,
        function(Sigma1) matrix_l1_norm(Sigma1 - Sigma0),
        numeric(1L)
      )
      
      return(list(
        p = p,
        eta_values = eta_values,
        Sigma0 = Sigma0,
        Sigma1_by_eta = Sigma1_by_eta,
        M_by_eta = M_by_eta,
        Xi_t = Xi_t,
        Xi_s = Xi_s,
        U = U,
        D = D,
        c = c_value,
        delta = delta_values,
        attempts = attempt
      ))
    }
  }
  
  stop(
    "Unable to generate a positive-definite population path after ",
    max_attempts, " attempts for p = ", p, "."
  )
}

dat.g.band <- function(n0, n1, p, eta, C_t, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  population <- make_population_path(p = p, eta_values = eta, C_t_value = C_t)
  Sigma0 <- population$Sigma0
  Sigma1 <- population$Sigma1_by_eta[[1L]]
  
  target_data <- rmvnorm_chol(as.integer(n0), Sigma0)
  source_data <- rmvnorm_chol(as.integer(n1), Sigma1)
  
  list(
    target_data = target_data,
    source_data = source_data,
    Sigma0 = Sigma0,
    Sigma1 = Sigma1
  )
}
