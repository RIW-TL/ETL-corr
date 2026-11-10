library(mvtnorm)
#============================================================
# Simulation 5: data generation
#============================================================

generate_sim5_data <- function(
    k0,
    d,
    p = 100,
    K = 5,
    s0 = 10,
    n0 = 100,
    n1 = 500,
    min_eig = 0.1,
    param_seed = 1,
    seed = 1) {
  

  
  set.seed(param_seed)
  #==========================================================
  # 1. Target precision matrix
  #==========================================================
  
  A <- diag(p) +
    matrix(
      runif(p * p, 0, 0.2),
      p, p
    )
  
  A <- (A + t(A)) / 2
  
  
  #----------------------------------------------------------
  # Keep s0 largest off-diagonal entries in each row
  #----------------------------------------------------------
  
  keep <- diag(TRUE, p)
  
  for (j in seq_len(p)) {
    
    ind <- setdiff(
      seq_len(p),
      j
    )
    
    id <- ind[
      order(
        abs(A[j, ind]),
        decreasing = TRUE
      )[seq_len(min(s0, p - 1))]
    ]
    
    keep[j, id] <- TRUE
  }
  
  keep <- keep | t(keep)
  
  A_sparse <- A * keep
  
  
  #----------------------------------------------------------
  # Positive-definite adjustment
  #----------------------------------------------------------
  
  lambda_min <- min(
    eigen(
      A_sparse,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  c0 <- max(
    0,
    min_eig - lambda_min
  )
  
  Omega0 <- A_sparse +
    c0 * diag(p)
  
  
  #==========================================================
  # 2. Reference precision matrix
  #==========================================================
  
  R0 <- cov2cor(Omega0)
  
  D0 <- diag(
    sqrt(diag(Omega0))
  )
  
  
  # Coordinate shift by p/2
  
  shift <- p / 2
  
  perm <- c(
    (shift + 1):p,
    1:shift
  )
  
  R_N <- R0[
    perm,
    perm
  ]
  
  Omega_N <-
    2 *
    D0 %*%
    R_N %*%
    D0
  
  Omega_N <-
    (Omega_N + t(Omega_N)) / 2
  
  
  #==========================================================
  # 3. Source precision matrices
  #==========================================================
  
  OmegaS <- vector(
    "list",
    K
  )
  
  for (k in seq_len(K)) {
    
    if (k <= k0) {
      
      OmegaS[[k]] <-
        (1 - d) * Omega0 +
        d * Omega_N
      
    } else {
      
      OmegaS[[k]] <-
        Omega_N
    }
    
    OmegaS[[k]] <-
      (OmegaS[[k]] +
         t(OmegaS[[k]])) / 2
  }
  
  
  #==========================================================
  # 4. Corresponding covariance matrices
  #==========================================================
  
  Sigma0 <- solve(Omega0)
  
  SigmaS <- lapply(
    OmegaS,
    solve
  )
  
  
  #==========================================================
  # 5. Generate target data
  #==========================================================
  set.seed(seed)
  X0 <- mvtnorm::rmvnorm(
    n = n0,
    mean = rep(0, p),
    sigma = Sigma0
  )
  
  
  #==========================================================
  # 6. Generate source data
  #==========================================================
  
  XS <- vector(
    "list",
    K
  )
  
  for (k in seq_len(K)) {
    
    XS[[k]] <- mvtnorm::rmvnorm(
      n = n1,
      mean = rep(0, p),
      sigma = SigmaS[[k]]
    )
  }
  
  
  #==========================================================
  # Output
  #==========================================================
  
  list(
    X0 = X0,
    XS = XS,
    Omega0 = Omega0,
    OmegaS = OmegaS
  )
}