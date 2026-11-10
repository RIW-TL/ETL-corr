library(ks)
library(mvtnorm)
library(parallel)

#============================================================
# Positive-definite adjustment
#============================================================

make_pd <- function(S, min_eig = 0.1, eps = 1e-6) {
  
  # Symmetrize
  S <- (S + t(S)) / 2
  
  # Minimum eigenvalue
  lambda_min <- min(
    eigen(S, symmetric = TRUE, only.values = TRUE)$values
  )
  
  # Add c * I if necessary
  c_add <- max(0, min_eig - lambda_min + eps)
  
  S <- S + c_add * diag(nrow(S))
  
  return(S)
}


#============================================================
# Loss function
#============================================================

cov_loss_fun <- function(Sigma_hat,
                         Sigma_true,
                         corr_hat = NULL,
                         var_hat = NULL) {
  
  p <- nrow(Sigma_true)
  
  #----------------------------------------------------------
  # True correlation matrix and variance vector
  #----------------------------------------------------------
  R_true <- cov2cor(Sigma_true)
  var_true <- diag(Sigma_true)
  
  #----------------------------------------------------------
  # For ordinary covariance estimators:
  # obtain correlation/variance by decomposition
  #----------------------------------------------------------
  if (is.null(corr_hat)) {
    corr_hat <- cov2cor(Sigma_hat)
  }
  
  if (is.null(var_hat)) {
    var_hat <- diag(Sigma_hat)
  }
  
  #----------------------------------------------------------
  # Losses
  #----------------------------------------------------------
  
  cov_loss <- sum((Sigma_hat - Sigma_true)^2) / p
  
  cor_loss <- sum((corr_hat - R_true)^2) / p
  
  var_loss <- sum((var_hat - var_true)^2)
  
  return(
    c(
      cov = cov_loss,
      cor = cor_loss,
      var = var_loss
    )
  )
}

compare.fun <- function(eta,
                        c_A1 = 0.9,
                        c_M1 = 0.5,
                        c_M2 = 0.1,
                        iter = 50,
                        ncores = 8) {
  
  #==========================================================
  # Persistent parallel cluster
  #
  # Used only by ETL-corr.
  # Create once and reuse for all Monte Carlo replications.
  #==========================================================
  
  cl <- NULL
  
  if (ncores > 1L) {
    
    cl <- parallel::makeCluster(ncores)
    
    user_lib <- path.expand("~/R/library")
    
    parallel::clusterCall(
      cl,
      function(lib) {
        .libPaths(c(lib, .libPaths()))
        invisible(NULL)
      },
      user_lib
    )
  }
  
  on.exit(
    {
      if (!is.null(cl)) {
        parallel::stopCluster(cl)
      }
    },
    add = TRUE
  )
  
  
  #==========================================================
  # Method names
  #==========================================================
  
  cov_methods <- c(
    "ETL-cov",
    "ETL-Scov",
    "ETL-Scov-adaptive",
    "Target-band",
    "ADT"
  )
  
  cor_methods <- c(
    "ETL-cov",
    "ETL-Scov",
    "ETL-Scov-adaptive",
    "Target-band",
    "ADT"
  )
  
  var_methods <- c(
    "ETL-cov",
    "Target-band",
    "ADT"
  )
  
  
  #==========================================================
  # Containers
  #==========================================================
  
  result.cov <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(cov_methods),
    dimnames = list(NULL, cov_methods)
  )
  
  result.cor <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(cor_methods),
    dimnames = list(NULL, cor_methods)
  )
  
  result.var <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(var_methods),
    dimnames = list(NULL, var_methods)
  )
  
  
  #==========================================================
  # Monte Carlo replications
  #==========================================================
  
  for (ii in seq_len(iter)) {
    
    #========================================================
    # Data generation
    #========================================================
    
    dat <- generate_sim2_data(
      eta = eta,
      seed = ii
    )
    
    X0 <- dat$X0
    XS <- dat$XS
    Sigma_true <- dat$Sigma0
    
    
    #========================================================
    # 1. ETL-corr
    #========================================================
    
    # ETL-corr
    fit <- fit.ETL.corr(X0, XS, c_A1 = c_A1, c_M1 = c_M1, c_M2 = c_M2, ncores = ncores)
    ETL_corr_hat = fit$R_hat

    #========================================================
    # 2. ETL-Scor
    #========================================================
    H_min = fit$H_min
    ETL_scor_hat = fit.ETL.Scor(ETL_corr_hat,H_min,c_T = 1)
    
    H_pair = fit$H_pair
    pair_id = fit$pair_id
    ETL_scor_hat_adaptive = fit.ETL.Scor.adaptive(R_hat = ETL_corr_hat, H_pair, pair_id)
    
    
    #========================================================
    # 3. ETL-var
    #========================================================
    
    fit <- ETL_var(
      X0 = X0,
      XS = XS,
      M = M,
      ncores = ncores
    )
    ETL_var_hat = fit$ETL_var_hat
    
    
    #========================================================
    # 4. ETL-cov
    #========================================================
    
    D_ETL <- diag(
      sqrt(ETL_var_hat)
    )
    
    ETL_cov_hat <-
      D_ETL %*%
      ETL_corr_hat %*%
      D_ETL
    
    ETL_cov_hat <- make_pd(
      ETL_cov_hat,
      min_eig = 0.1
    )
    
    
    #========================================================
    # 5. ETL-Scov
    #========================================================
    
    ETL_scov_hat <-
      D_ETL %*%
      ETL_scor_hat %*%
      D_ETL
    
    ETL_scov_hat <- make_pd(
      ETL_scov_hat,
      min_eig = 0.1
    )
    
    #========================================================
    # 5. ETL-Scov-adaptive
    #========================================================
    
    ETL_scov_hat_adaptive <-
      D_ETL %*%
      ETL_scor_hat_adaptive %*%
      D_ETL
    
    ETL_scov_hat_adaptive <- make_pd(
      ETL_scov_hat_adaptive,
      min_eig = 0.1
    )
    
    
    #========================================================
    # 6. Target-band
    #========================================================
    
    fit_target_band <- fit_band(
      target_data = X0
    )
    
    Target_band_hat <- fit_target_band$estimate
    Target_band_cor_hat <- cov2cor(Target_band_hat)
    Target_band_var_hat <- diag(Target_band_hat)
    
    
    #========================================================
    # 7. ADT
    #========================================================
    
    fit_adt <- ADT(
      target_data = X0,
      source_data = XS[[1L]]
    )
    
    ADT_hat <- fit_adt$estimate
    ADT_cor_hat <- cov2cor(ADT_hat)
    ADT_var_hat <- diag(ADT_hat)
    
    
    #========================================================
    # Compute losses
    #========================================================
    
    p <- ncol(Sigma_true)
    
    R_true <- cov2cor(Sigma_true)
    var_true <- diag(Sigma_true)
    
    
    #========================================================
    # ETL-cov: direct calculation
    #========================================================
    
    loss_ETL_cov <- c(
      cov = norm(
        ETL_cov_hat - Sigma_true,
        type = "F"
      )^2 / p,
      
      cor = norm(
        ETL_corr_hat - R_true,
        type = "F"
      )^2 / p,
      
      var = sum(
        (ETL_var_hat - var_true)^2
      )
    )
    
    
    #========================================================
    # ETL-Scov: direct calculation
    #========================================================
    
    loss_ETL_scov <- c(
      cov = norm(
        ETL_scov_hat - Sigma_true,
        type = "F"
      )^2 / p,
      
      cor = norm(
        ETL_scor_hat - R_true,
        type = "F"
      )^2 / p,
      
      var = sum(
        (ETL_var_hat - var_true)^2
      )
    )
    
    #========================================================
    # ETL-Scov-adaptive: direct calculation
    #========================================================
    
    loss_ETL_scov_adaptive <- c(
      cov = norm(
        ETL_scov_hat_adaptive - Sigma_true,
        type = "F"
      )^2 / p,
      
      cor = norm(
        ETL_scor_hat_adaptive - R_true,
        type = "F"
      )^2 / p,
      
      var = sum(
        (ETL_var_hat - var_true)^2
      )
    )
    
    
    #========================================================
    # Target-band: keep original loss function
    #========================================================
    
    loss_Target_band <- cov_loss_fun(
      Sigma_hat = Target_band_hat,
      Sigma_true = Sigma_true,
      corr_hat = Target_band_cor_hat,
      var_hat = Target_band_var_hat
    )
    
    
    #========================================================
    # ADT: keep original loss function
    #========================================================
    
    loss_ADT <- cov_loss_fun(
      Sigma_hat = ADT_hat,
      Sigma_true = Sigma_true,
      corr_hat = ADT_cor_hat,
      var_hat = ADT_var_hat
    )
    
    
    #========================================================
    # Store covariance losses
    #========================================================
    
    result.cov[ii, ] <- c(
      loss_ETL_cov["cov"],
      loss_ETL_scov["cov"],
      loss_ETL_scov_adaptive["cov"],
      loss_Target_band["cov"],
      loss_ADT["cov"]
    )
    
    
    #========================================================
    # Store correlation losses
    #========================================================
    
    result.cor[ii, ] <- c(
      loss_ETL_cov["cor"],
      loss_ETL_scov["cor"],
      loss_ETL_scov_adaptive["cor"],
      loss_Target_band["cor"],
      loss_ADT["cor"]
    )
    
    
    #========================================================
    # Store variance losses
    #========================================================
    
    result.var[ii, ] <- c(
      loss_ETL_cov["var"],
      loss_Target_band["var"],
      loss_ADT["var"]
    )
  }
  
  
  #==========================================================
  # Mean losses
  #==========================================================
  
  mean.cov <- colMeans(
    result.cov,
    na.rm = TRUE
  )
  
  mean.cor <- colMeans(
    result.cor,
    na.rm = TRUE
  )
  
  mean.var <- colMeans(
    result.var,
    na.rm = TRUE
  )
  
  
  #==========================================================
  # Output
  #==========================================================
  
  list(
    cov_loss = result.cov,
    cor_loss = result.cor,
    var_loss = result.var,
    
    mean_cov_loss = mean.cov,
    mean_cor_loss = mean.cor,
    mean_var_loss = mean.var
  )
}

s1  <- compare.fun(eta = 0.1)
# s2  <- compare.fun(eta = 0.2)
s3  <- compare.fun(eta = 0.3)
# s4  <- compare.fun(eta = 0.4)
s5  <- compare.fun(eta = 0.5)
# s6  <- compare.fun(eta = 0.6)
s7  <- compare.fun(eta = 0.7)
# s8  <- compare.fun(eta = 0.8)
s9  <- compare.fun(eta = 0.9)
# s10 <- compare.fun(eta = 1.0)

rbind(
  s1$mean_cor_loss,
  # s2$mean_cor_loss,
  s3$mean_cor_loss,
  # s4$mean_cor_loss,
  s5$mean_cor_loss,
  # s6$mean_cor_loss,
  s7$mean_cor_loss,
  # s8$mean_cor_loss,
  s9$mean_cor_loss
  # s10$mean_cor_loss
)

rbind(
  s1$mean_cov_loss,
  # s2$mean_cov_loss,
  s3$mean_cov_loss,
  # s4$mean_cov_loss,
  s5$mean_cov_loss,
  # s6$mean_cov_loss,
  s7$mean_cov_loss,
  # s8$mean_cov_loss,
  s9$mean_cov_loss
)

rbind(
  s1$mean_var_loss,
  # s2$mean_var_loss,
  s3$mean_var_loss,
  # s4$mean_var_loss,
  s5$mean_var_loss,
  # s6$mean_var_loss,
  s7$mean_var_loss,
  # s8$mean_var_loss,
  s9$mean_var_loss
)
