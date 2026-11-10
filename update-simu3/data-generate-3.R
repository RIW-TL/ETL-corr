
library(mvtnorm)

generate_sim3_data <- function(
    p = 100,
    n0 = 100,
    n1 = 500,
    K = 5,
    eta = 0.4,
    d = 0.6,
    tau = 0,
    L = 10,
    seed = NULL) {
  
  #==========================================================
  # 0. Input checks
  #==========================================================
  
  stopifnot(
    p >= 2,
    p %% L == 0,
    n0 > 0,
    n1 > 0,
    K >= 1,
    eta >= 0,
    eta <= 1,
    d >= 0,
    d < 1,
    tau >= 0,
    tau < 1
  )
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  #==========================================================
  # 1. Target correlation matrix
  #==========================================================
  
  block_size <- p / L
  block_id <- rep(1:L, each = block_size)
  
  # u_j ~ Unif(0.4, 0.6)
  u <- runif(
    p,
    min = 0.6,
    max = 0.8
  )
  
  # Rademacher variables
  xi0 <- sample(
    c(-1, 1),
    p,
    replace = TRUE
  )
  
  # B^(0)
  B0 <- matrix(
    0,
    nrow = p,
    ncol = L
  )
  
  B0[
    cbind(1:p, block_id)
  ] <- xi0 * u
  
  # R^(0)
  R0 <-
    B0 %*% t(B0) +
    diag(1 - u^2)
  
  
  #==========================================================
  # 2. Target marginal variances and covariance matrix
  #==========================================================
  
  v <- runif(
    p,
    min = -0.3,
    max = 0.3
  )
  
  var0 <- exp(v)
  
  Lambda0 <- diag(var0)
  
  D0 <- diag(sqrt(var0))
  
  Sigma0 <-
    D0 %*%
    R0 %*%
    D0
  
  
  #==========================================================
  # 3. Source covariance matrices
  #==========================================================
  
  RS <- vector("list", K)
  varS <- vector("list", K)
  LambdaS <- vector("list", K)
  SigmaS <- vector("list", K)
  
  S_R <- vector("list", K)
  S_Lambda <- vector("list", K)
  
  # approximately eta proportion of correlation pairs
  # and eta proportion of marginal variances are matched
  s_R <- floor(p * sqrt(eta))
  s_Lambda <- floor(p * eta)
  
  
  for (k in 1:K) {
    
    #--------------------------------------------------------
    # Locally informative coordinates
    #--------------------------------------------------------
    
    S_R[[k]] <- if (s_R > 0) {
      sample(
        1:p,
        size = s_R,
        replace = FALSE
      )
    } else {
      integer(0)
    }
    
    S_Lambda[[k]] <- if (s_Lambda > 0) {
      sample(
        1:p,
        size = s_Lambda,
        replace = FALSE
      )
    } else {
      integer(0)
    }
    
    
    #--------------------------------------------------------
    # Source correlation matrix
    #--------------------------------------------------------
    
    Ck <- matrix(
      1 - d,
      nrow = p,
      ncol = p
    )
    
    diag(Ck) <- 1
    
    if (s_R > 0) {
      Ck[
        S_R[[k]],
        S_R[[k]]
      ] <- 1
    }
    
    Rk <- R0 * Ck
    
    diag(Rk) <- 1
    
    RS[[k]] <- Rk
    
    
    #--------------------------------------------------------
    # Source marginal variances
    #--------------------------------------------------------
    
    xi_k <- sample(
      c(-1, 1),
      p,
      replace = TRUE
    )
    
    unmatched <- !(
      1:p %in% S_Lambda[[k]]
    )
    
    var_k <-
      var0 *
      exp(
        d *
          xi_k *
          as.numeric(unmatched)
      )
    
    varS[[k]] <- var_k
    LambdaS[[k]] <- diag(var_k)
    
    
    #--------------------------------------------------------
    # Source covariance matrix
    #--------------------------------------------------------
    
    Dk <- diag(sqrt(var_k))
    
    SigmaS[[k]] <-
      Dk %*%
      Rk %*%
      Dk
  }
  
  
  #==========================================================
  # 4. Gaussian / Gaussian-mixture generator
  #
  # tau = 0:
  #     X ~ N_p(0, Sigma)
  #
  # tau > 0:
  #     X ~ 1/2 N_p(0, (1-tau) Sigma)
  #       + 1/2 N_p(0, (1+tau) Sigma)
  #==========================================================
  
  generate_X <- function(
    n,
    Sigma,
    tau) {
    
    # Standard Gaussian observations
    Z <- matrix(
      rnorm(n * p),
      nrow = n,
      ncol = p
    )
    
    # chol(Sigma): t(A) %*% A = Sigma
    A <- chol(Sigma)
    
    X <- Z %*% A
    
    if (tau > 0) {
      
      # Each observation independently chooses
      # one of the two Gaussian components
      component <- sample(
        c(-1, 1),
        size = n,
        replace = TRUE
      )
      
      # component = -1: covariance (1-tau) Sigma
      # component =  1: covariance (1+tau) Sigma
      scale_i <- sqrt(
        1 + tau * component
      )
      
      X <- sweep(
        X,
        MARGIN = 1,
        STATS = scale_i,
        FUN = "*"
      )
    }
    
    X
  }
  
  
  #==========================================================
  # 5. Generate target and source data
  #==========================================================
  
  X0 <- generate_X(
    n = n0,
    Sigma = Sigma0,
    tau = tau
  )
  
  XS <- lapply(
    SigmaS,
    function(Sigma_k) {
      
      generate_X(
        n = n1,
        Sigma = Sigma_k,
        tau = tau
      )
    }
  )
  
  
  #==========================================================
  # 6. Output
  #==========================================================
  
  list(
    X0 = X0,
    XS = XS,
    
    Sigma0 = Sigma0,
    SigmaS = SigmaS,
    
    R0 = R0,
    RS = RS,
    
    var0 = var0,
    varS = varS
  )
}
