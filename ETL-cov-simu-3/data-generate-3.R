############################################################
## 1. Generate target/source covariance parameters
############################################################

generate_etl_params <- function(
    p = 50,
    K = 5,
    G = 5,
    eta = 0.4,
    d = 0.5,
    u_range = c(0.4, 0.6),
    v_range = c(-0.3, 0.3),
    seed = 2026
) {
  
  set.seed(seed)
  
  #----------------------------------------------------------
  # Basic checks
  #----------------------------------------------------------
  stopifnot(
    p >= 2,
    K >= 1,
    G >= 1,
    G <= p,
    eta >= 0,
    eta <= 1,
    d >= 0,
    d <= 1
  )
  
  ############################################################
  ## Target correlation matrix R0
  ############################################################
  
  # Approximately equal-sized blocks
  block_id <- rep(1:G, length.out = p)
  
  # Random loadings
  u <- runif(
    p,
    min = u_range[1],
    max = u_range[2]
  )
  
  sign_vec <- sample(
    c(-1, 1),
    size = p,
    replace = TRUE
  )
  
  a <- sign_vec * u
  
  # B0: each coordinate loads on only one block factor
  B0 <- matrix(
    0,
    nrow = p,
    ncol = G
  )
  
  for (j in 1:p) {
    B0[j, block_id[j]] <- a[j]
  }
  
  # Correlation matrix:
  # R0 = B0 B0^T + diag(1-u_j^2)
  R0 <- tcrossprod(B0) +
    diag(1 - u^2, p)
  
  # Numerical symmetry
  R0 <- (R0 + t(R0)) / 2
  
  
  ############################################################
  ## Target variance matrix Lambda0
  ############################################################
  
  v <- runif(
    p,
    min = v_range[1],
    max = v_range[2]
  )
  
  var0 <- exp(v)
  
  Lambda0 <- diag(var0, p)
  
  # Target covariance
  sqrtLambda0 <- diag(sqrt(var0), p)
  
  Sigma0 <- sqrtLambda0 %*%
    R0 %*%
    sqrtLambda0
  
  Sigma0 <- (Sigma0 + t(Sigma0)) / 2
  
  
  ############################################################
  ## Source covariance matrices
  ############################################################
  
  # Number of coordinates involved in informative
  # correlation and variance components
  s_R <- floor(p * sqrt(eta))
  s_Lambda <- floor(p * eta)
  
  R_source <- vector(
    "list",
    K
  )
  
  Lambda_source <- vector(
    "list",
    K
  )
  
  Sigma_source <- vector(
    "list",
    K
  )
  
  var_source <- vector(
    "list",
    K
  )
  
  # Store informative coordinate sets
  S_R <- vector(
    "list",
    K
  )
  
  S_Lambda <- vector(
    "list",
    K
  )
  
  # Entry-wise informative masks
  informative_R_mask <- vector(
    "list",
    K
  )
  
  # Coordinate-wise informative masks
  informative_Lambda_mask <- vector(
    "list",
    K
  )
  
  
  for (k in 1:K) {
    
    ##########################################################
    ## Correlation component
    ##########################################################
    
    # Informative coordinate subset for source k
    if (s_R > 0) {
      S_R[[k]] <- sort(
        sample(
          1:p,
          size = s_R,
          replace = FALSE
        )
      )
    } else {
      S_R[[k]] <- integer(0)
    }
    
    # C^(k)
    # c_jl = 1 if j=l or both j,l are informative
    #        = 1-d otherwise
    Ck <- matrix(
      1 - d,
      nrow = p,
      ncol = p
    )
    
    diag(Ck) <- 1
    
    if (length(S_R[[k]]) > 0) {
      Ck[
        S_R[[k]],
        S_R[[k]]
      ] <- 1
    }
    
    # Source correlation matrix
    Rk <- R0 * Ck
    
    Rk <- (Rk + t(Rk)) / 2
    
    R_source[[k]] <- Rk
    
    # Informative entry mask
    mask_R <- matrix(
      FALSE,
      nrow = p,
      ncol = p
    )
    
    if (length(S_R[[k]]) > 0) {
      mask_R[
        S_R[[k]],
        S_R[[k]]
      ] <- TRUE
    }
    
    diag(mask_R) <- TRUE
    
    informative_R_mask[[k]] <- mask_R
    
    
    ##########################################################
    ## Variance component
    ##########################################################
    
    # Independently generate informative variance coordinates
    if (s_Lambda > 0) {
      S_Lambda[[k]] <- sort(
        sample(
          1:p,
          size = s_Lambda,
          replace = FALSE
        )
      )
    } else {
      S_Lambda[[k]] <- integer(0)
    }
    
    informative_var <- rep(
      FALSE,
      p
    )
    
    informative_var[
      S_Lambda[[k]]
    ] <- TRUE
    
    informative_Lambda_mask[[k]] <-
      informative_var
    
    # +/- 1 perturbation direction
    xi <- sample(
      c(-1, 1),
      size = p,
      replace = TRUE
    )
    
    # Informative coordinates:
    # sigma_jj^(k) = sigma_jj^(0)
    #
    # Noninformative coordinates:
    # sigma_jj^(k) = sigma_jj^(0) exp(d xi_jk)
    var_k <- var0 *
      exp(
        d *
          xi *
          (!informative_var)
      )
    
    Lambda_k <- diag(
      var_k,
      p
    )
    
    var_source[[k]] <- var_k
    Lambda_source[[k]] <- Lambda_k
    
    
    ##########################################################
    ## Source covariance matrix
    ##########################################################
    
    sqrtLambda_k <- diag(
      sqrt(var_k),
      p
    )
    
    Sigma_k <- sqrtLambda_k %*%
      Rk %*%
      sqrtLambda_k
    
    Sigma_k <- (
      Sigma_k +
        t(Sigma_k)
    ) / 2
    
    Sigma_source[[k]] <- Sigma_k
  }
  
  
  ############################################################
  ## Some useful diagnostic quantities
  ############################################################
  
  # Number of informative sources for each correlation entry
  N_R <- matrix(
    0,
    nrow = p,
    ncol = p
  )
  
  for (k in 1:K) {
    N_R <- N_R +
      informative_R_mask[[k]]
  }
  
  # Number of informative sources for each variance coordinate
  N_Lambda <- rep(
    0,
    p
  )
  
  for (k in 1:K) {
    N_Lambda <- N_Lambda +
      informative_Lambda_mask[[k]]
  }
  
  
  ############################################################
  ## Positive-definiteness diagnostics
  ############################################################
  
  min_eig_R0 <- min(
    eigen(
      R0,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  min_eig_Sigma0 <- min(
    eigen(
      Sigma0,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  min_eig_R_source <- sapply(
    R_source,
    function(A) {
      min(
        eigen(
          A,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      )
    }
  )
  
  min_eig_Sigma_source <- sapply(
    Sigma_source,
    function(A) {
      min(
        eigen(
          A,
          symmetric = TRUE,
          only.values = TRUE
        )$values
      )
    }
  )
  
  
  ############################################################
  ## Output
  ############################################################
  
  return(
    list(
      
      # Basic settings
      p = p,
      K = K,
      G = G,
      eta = eta,
      d = d,
      
      # Target parameters
      R0 = R0,
      Lambda0 = Lambda0,
      Sigma0 = Sigma0,
      var0 = var0,
      
      # Source parameters
      R_source = R_source,
      Lambda_source = Lambda_source,
      Sigma_source = Sigma_source,
      var_source = var_source,
      
      # Informative sets
      S_R = S_R,
      S_Lambda = S_Lambda,
      
      # Informative masks
      informative_R_mask = informative_R_mask,
      informative_Lambda_mask =
        informative_Lambda_mask,
      
      # Number of informative sources
      N_R = N_R,
      N_Lambda = N_Lambda,
      
      # Auxiliary quantities
      block_id = block_id,
      B0 = B0,
      u = u,
      sign_vec = sign_vec,
      
      # Actual informative-set sizes
      s_R = s_R,
      s_Lambda = s_Lambda,
      
      # Positive-definiteness checks
      min_eig_R0 = min_eig_R0,
      min_eig_Sigma0 = min_eig_Sigma0,
      min_eig_R_source =
        min_eig_R_source,
      min_eig_Sigma_source =
        min_eig_Sigma_source
    )
  )
}
############################################################
## 2. Generate target/source observations
############################################################

generate_etl_data <- function(
    params,
    n0 = 100,
    nk = 500,
    distribution = c(
      "scale_mixture",
      "gaussian"
    ),
    tau = 0.4,
    seed = 1
) {
  
  set.seed(seed)
  
  distribution <- match.arg(
    distribution
  )
  
  p <- params$p
  K <- params$K
  
  Sigma0 <- params$Sigma0
  Sigma_source <- params$Sigma_source
  
  
  ############################################################
  ## Source sample sizes
  ############################################################
  
  if (length(nk) == 1) {
    nk <- rep(
      nk,
      K
    )
  }
  
  if (length(nk) != K) {
    stop(
      "nk must be either a scalar or a vector of length K."
    )
  }
  
  
  ############################################################
  ## Helper: generate observations with covariance Sigma
  ############################################################
  
  generate_one_domain <- function(
    n,
    Sigma
  ) {
    
    # Cholesky factor:
    # t(L) %*% L = Sigma
    L <- chol(Sigma)
    
    Z <- matrix(
      rnorm(n * p),
      nrow = n,
      ncol = p
    )
    
    if (distribution == "gaussian") {
      
      # X ~ N(0, Sigma)
      X <- Z %*% L
      
    } else {
      
      # Bounded Gaussian scale mixture
      #
      # V takes 1-tau and 1+tau with equal probability.
      # E(V) = 1, hence Cov(X) = Sigma.
      
      if (tau <= 0 || tau >= 1) {
        stop(
          "For the scale-mixture design, tau must satisfy 0 < tau < 1."
        )
      }
      
      V <- sample(
        c(
          1 - tau,
          1 + tau
        ),
        size = n,
        replace = TRUE,
        prob = c(0.5, 0.5)
      )
      
      X <- Z %*% L
      
      X <- X *
        sqrt(V)
    }
    
    return(X)
  }
  
  
  ############################################################
  ## Target data
  ############################################################
  
  X0 <- generate_one_domain(
    n = n0,
    Sigma = Sigma0
  )
  
  
  ############################################################
  ## Source data
  ############################################################
  
  XS <- vector(
    "list",
    K
  )
  
  for (k in 1:K) {
    
    XS[[k]] <- generate_one_domain(
      n = nk[k],
      Sigma = Sigma_source[[k]]
    )
  }
  
  
  ############################################################
  ## Output
  ############################################################
  
  return(
    list(
      X0 = X0,
      XS = XS,
      n0 = n0,
      nk = nk,
      distribution = distribution,
      tau = if (
        distribution == "scale_mixture"
      ) tau else NA,
      
      # True parameters
      Sigma0 = Sigma0,
      Sigma_source = Sigma_source,
      R0 = params$R0,
      R_source = params$R_source,
      Lambda0 = params$Lambda0,
      Lambda_source =
        params$Lambda_source
    )
  )
}


# params = generate_etl_params(eta = 0.4, d = 0.5)
# fit = generate_etl_data(params, distribution = "gaussian", tau = 0, seed = ii)
# X0 = fit$X0
# XS = fit$XS
