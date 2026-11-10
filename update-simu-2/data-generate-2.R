library(mvtnorm)
#-----------------------------------------------------
generate_sim2_data <- function(
    p = 100,
    n0 = 100,
    n1 = 500,
    eta,
    tau = 0.2,
    C_t = 0.2,
    delta_max = 0.5,
    seed = NULL,
    max_attempts = 5000L,
    eigen_tolerance = 1e-8) {
  

  #==========================================================
  # Auxiliary quantities
  #==========================================================
  
  matrix_l1_norm <- function(A) {
    max(colSums(abs(A)))
  }
  
  min_eigenvalue <- function(A) {
    min(
      eigen(
        A,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    )
  }
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  idx <- 1:p
  gap <- abs(outer(idx, idx, "-"))
  off_diag <- gap > 0
  
  
  #==========================================================
  # Generate population covariance matrices
  #==========================================================
  
    #--------------------------------------------------------
    # Target Rademacher matrix
    #--------------------------------------------------------
    
    Xi_t <- matrix(0, p, p)
    
    Xi_t[upper.tri(Xi_t)] <- sample(
      c(-1, 1),
      p * (p - 1) / 2,
      replace = TRUE
    )
    
    Xi_t <- Xi_t + t(Xi_t)
    
    
    #--------------------------------------------------------
    # Source perturbation Rademacher matrix
    #--------------------------------------------------------
    
    Xi_s <- matrix(0, p, p)
    
    Xi_s[upper.tri(Xi_s)] <- sample(
      c(-1, 1),
      p * (p - 1) / 2,
      replace = TRUE
    )
    
    Xi_s <- Xi_s + t(Xi_s)
    
    
    #--------------------------------------------------------
    # Uniform matrix U
    #--------------------------------------------------------
    
    U <- matrix(0, p, p)
    
    U[upper.tri(U)] <- runif(
      p * (p - 1) / 2
    )
    
    U <- U + t(U)
    diag(U) <- 1
    
    
    #========================================================
    # Target covariance matrix
    #
    # Sigma0[j,j] = 1
    #
    # Sigma0[j,l]
    #   = C_t |j-l|^{-3} Xi_t[j,l]
    #========================================================
    
    Sigma0 <- diag(p)
    
    Sigma0[off_diag] <-
      C_t *
      gap[off_diag]^(-2) *
      Xi_t[off_diag]
    
    
    #========================================================
    # Perturbation matrix
    #
    # D[j,l]
    #   = |j-l|^{-2} Xi_s[j,l]
    #========================================================
    
    D <- matrix(0, p, p)
    
    D[off_diag] <-
      gap[off_diag]^(-2) *
      Xi_s[off_diag]
    
    
    #========================================================
    # Scaling constant
    #
    # ||c D||_1 = delta_max
    #========================================================
    
    c_value <-
      delta_max /
      matrix_l1_norm(D)
    
    
    #========================================================
    # Entry-wise perturbation indicators
    #
    # M[j,l](eta) = I(U[j,l] <= eta)
    #========================================================
    
    M <- 1L * (U <= eta)
    diag(M) <- 0L
    
    
    #========================================================
    # Source covariance matrix
    #
    # Sigma1 = Sigma0 + c (M o D)
    #========================================================
    
    perturbation <-
      c_value * (M * D)
    
    Sigma1 <-
      Sigma0 +
      perturbation + tau*diag(p)
    
    Sigma1 <-
      (Sigma1 + t(Sigma1)) / 2
    
    
    #========================================================
    # Positive-definiteness check
    #========================================================
    
  #   if (
  #     min_eigenvalue(Sigma0) > eigen_tolerance &&
  #     min_eigenvalue(Sigma1) > eigen_tolerance
  #   ) {
  #     break
  #   }
  #   
  #   
  #   if (attempt == max_attempts) {
  #     stop(
  #       "Unable to generate positive-definite covariance matrices ",
  #       "after ", max_attempts, " attempts."
  #     )
  #   }
  # }
  
  
  #==========================================================
  # Generate target and source data
  #==========================================================
  
  # X0 <- mvtnorm::rmvnorm(
  #   n = n0,
  #   mean = rep(0, p),
  #   sigma = Sigma0
  # )
  # 
  # X1 <- mvtnorm::rmvnorm(
  #   n = n1,
  #   mean = rep(0, p),
  #   sigma = Sigma1
  # )
  
  
  #==========================================================
  # Output
  #==========================================================
  
  list(
    # X0 = X0,
    # XS = list(X1),
    Sigma0 = Sigma0,
    Sigma1 = Sigma1
  )
}
