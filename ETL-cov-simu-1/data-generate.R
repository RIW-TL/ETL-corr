library(mvtnorm)

#============================================================
# Simulation 1: parameter generation
#============================================================

generate_sim1_params <- function(
    p = 100,
    K = 5,
    k0 = 3,
    d = 0.4,
    s0 = 10,
    min_eig = 0.1,
    seed = 1) {
  
  #----------------------------------------------------------
  # Input checks
  #----------------------------------------------------------
  
  if (k0 < 0 || k0 > K) {
    stop("k0 must be between 0 and K.")
  }
  
  if (d < 0 || d > 1) {
    stop("d must be between 0 and 1.")
  }
  
  if (p %% 2 != 0) {
    stop("p must be even for the p/2 coordinate shift.")
  }
  
  set.seed(seed)
  
  
  #==========================================================
  # 1. Target covariance matrix
  #==========================================================
  
  A <- diag(p) + matrix(
    runif(p * p, 0, 0.2),
    p, p
  )
  
  A <- (A + t(A)) / 2
  
  
  #----------------------------------------------------------
  # Keep s0 largest off-diagonal entries in each row
  #----------------------------------------------------------
  
  keep <- diag(TRUE, p)
  
  for (j in seq_len(p)) {
    
    ind <- setdiff(seq_len(p), j)
    
    id <- ind[
      order(
        abs(A[j, ind]),
        decreasing = TRUE
      )[seq_len(min(s0, p - 1))]
    ]
    
    keep[j, id] <- TRUE
  }
  
  # Keep symmetric counterparts
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
  
  c0 <- max(0, min_eig - lambda_min)
  
  Sigma0 <- A_sparse + c0 * diag(p)
  
  
  #==========================================================
  # 2. Target correlation matrix and marginal variances
  #==========================================================
  
  R0 <- cov2cor(Sigma0)
  var0 <- diag(Sigma0)
  
  D0 <- diag(sqrt(var0))
  
  
  #==========================================================
  # 3. Reference covariance matrix Sigma_N
  #
  # P shifts coordinates by p/2:
  # (1,...,p/2,p/2+1,...,p)
  # ->
  # (p/2+1,...,p,1,...,p/2)
  #==========================================================
  
  shift <- p / 2
  
  perm <- c(
    (shift + 1):p,
    1:shift
  )
  
  # P R0 P^T
  R_N <- R0[perm, perm]
  
  # Reference covariance:
  # different correlation structure
  # + marginal variances twice those of the target
  Sigma_N <- 2 * D0 %*% R_N %*% D0
  
  # Numerical symmetrization
  Sigma_N <- (Sigma_N + t(Sigma_N)) / 2
  
  
  #==========================================================
  # 4. Source covariance matrices
  #==========================================================
  
  SigmaS <- vector("list", K)
  
  for (k in seq_len(K)) {
    
    if (k <= k0) {
      
      # Potentially informative sources
      SigmaS[[k]] <-
        (1 - d) * Sigma0 +
        d * Sigma_N
      
    } else {
      
      # Non-informative sources
      SigmaS[[k]] <- Sigma_N
    }
    
    # Numerical symmetrization
    SigmaS[[k]] <-
      (SigmaS[[k]] + t(SigmaS[[k]])) / 2
  }
  
  
  #==========================================================
  # 5. Population correlation matrices and variances
  #==========================================================
  
  RS <- lapply(
    SigmaS,
    cov2cor
  )
  
  varS <- lapply(
    SigmaS,
    diag
  )
  
  
  #==========================================================
  # Return population quantities
  #==========================================================
  
  list(
    Sigma0 = Sigma0,
    SigmaS = SigmaS,
    
    R0 = R0,
    RS = RS,
    
    var0 = var0,
    varS = varS,
    
    Sigma_N = Sigma_N,
    R_N = R_N,
    
    p = p,
    K = K,
    k0 = k0,
    d = d
  )
}


#============================================================
# Simulation 1: data generation
#============================================================

generate_sim1_data <- function(
    param,
    n0 = 100,
    n_s_vec = rep(500, param$K),
    seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  p <- param$p
  K <- param$K
  
  
  #----------------------------------------------------------
  # Check source sample sizes
  #----------------------------------------------------------
  
  if (length(n_s_vec) != K) {
    stop("length(n_s_vec) must equal K.")
  }
  
  
  #==========================================================
  # Target data
  #==========================================================
  
  X0 <- mvtnorm::rmvnorm(
    n = n0,
    mean = rep(0, p),
    sigma = param$Sigma0
  )
  
  
  #==========================================================
  # Source data
  #==========================================================
  
  XS <- vector("list", K)
  
  for (k in seq_len(K)) {
    
    XS[[k]] <- mvtnorm::rmvnorm(
      n = n_s_vec[k],
      mean = rep(0, p),
      sigma = param$SigmaS[[k]]
    )
  }
  
  
  #==========================================================
  # Output
  #==========================================================
  
  list(
    X0 = X0,
    XS = XS
  )
}
