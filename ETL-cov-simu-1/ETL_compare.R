library(ks)
library(mvtnorm)
library(sparseCov)
library(parallel)

#============================================================
# Source project files
#============================================================
# setwd("C:/Users/AMA/Nutstore/1/Nutstore/numerical studies/Trans-ECM/ETL-cov-simu-1")

# source("data-generate.R")
# source("ETL-functions.R")
# source("ETL-var.R")
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


#============================================================
# Main comparison function
#============================================================

compare.fun <- function(k0, d,
                        c_A1 = 0.5,
                        c_M1 = 0.5,
                        c_M2 = 0.1,
                        iter = 50,
                        ncores = 18) {
  
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
    "Target-cov",
    "Target-Scov",
    "Pool-cov"
  )
  
  var_methods <- c(
    "ETL-var",
    "Target-var",
    "Pool-var"
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
    ncol = length(cov_methods),
    dimnames = list(NULL, cov_methods)
  )
  
  result.var <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(var_methods),
    dimnames = list(NULL, var_methods)
  )
  
  
  #==========================================================
  # Parameter setting
  #==========================================================
  
  param <- generate_sim1_params(
    p = 100,
    K = 5,
    k0 = k0,
    d = d,
    s0 = 10,
    seed = 1
  )
  
  p <- param$p
  K <- param$K
  
  Sigma_true <- param$Sigma0
  
  
  #==========================================================
  # Monte Carlo replications
  #==========================================================
  
  for (ii in seq_len(iter)) {
    
    #========================================================
    # Data generation
    #========================================================
    
    dat <- generate_sim1_data(
      param = param,
      n0 = 100,
      n_s_vec = rep(500, K),
      seed = ii
    )
    
    X0 <- dat$X0
    XS <- dat$XS
    
    
    #========================================================
    # 1. ETL-corr
    #========================================================
    
    fit_corr <- fit.ETL.corr(
      X0 = X0,
      XS = XS,
      c_A1 = c_A1,
      c_M1 = c_M1,
      c_M2 = c_M2,
      cl = cl
    )
    
    ETL_corr_hat <- fit_corr$R_hat
    H_min <- fit_corr$H_min
    
    
    #========================================================
    # 2. ETL-Scor
    #========================================================
    
    ETL_scor_hat <- fit.ETL.Scor(
      R_hat = ETL_corr_hat,
      H_min = H_min,
      c_T = 0.5
    )
    
    
    #========================================================
    # 3. ETL-var
    #========================================================
    
    fit <- ETL_var(
      X0 = X0,
      XS = XS
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
    # 6. Target-cov
    #========================================================
    
    Target_cov_hat <- cov(X0)
    
    
    #========================================================
    # 7. Pool-cov
    #========================================================
    
    X_all <- do.call(
      rbind,
      c(list(X0), XS)
    )
    
    Pool_cov_hat <- cov(X_all)
    
    
    #========================================================
    # 8. Target-Scov
    #========================================================
    
    Target_scov_hat <- fit.target.scor(X0, corr = FALSE)
    
    Target_scov_hat <- as.matrix(
      Target_scov_hat
    )
    
    
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
    
    loss_Target_cov <- cov_loss_fun(
      Sigma_hat = Target_cov_hat,
      Sigma_true = Sigma_true
    )
    
    loss_Target_scov <- cov_loss_fun(
      Sigma_hat = Target_scov_hat,
      Sigma_true = Sigma_true
    )
    
    loss_Pool_cov <- cov_loss_fun(
      Sigma_hat = Pool_cov_hat,
      Sigma_true = Sigma_true
    )
    
    
    #========================================================
    # Store covariance losses
    #========================================================
    
    result.cov[ii, ] <- c(
      loss_ETL_cov["cov"],
      loss_ETL_scov["cov"],
      loss_Target_cov["cov"],
      loss_Target_scov["cov"],
      loss_Pool_cov["cov"]
    )
    
    
    #========================================================
    # Store correlation losses
    #========================================================
    
    result.cor[ii, ] <- c(
      loss_ETL_cov["cor"],
      loss_ETL_scov["cor"],
      loss_Target_cov["cor"],
      loss_Target_scov["cor"],
      loss_Pool_cov["cor"]
    )
    
    
    #========================================================
    # Store variance losses
    #========================================================
    
    result.var[ii, ] <- c(
      loss_ETL_cov["var"],
      loss_Target_cov["var"],
      loss_Pool_cov["var"]
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


#============================================================
# Simulation 1
#============================================================


a0 <- compare.fun(k0 = 0, d = 0.2)
a1 <- compare.fun(k0 = 1, d = 0.2)
a2 <- compare.fun(k0 = 2, d = 0.2)
a3 <- compare.fun(k0 = 3, d = 0.2)
a4 <- compare.fun(k0 = 4, d = 0.2)
a5 <- compare.fun(k0 = 5, d = 0.2)

# b0 <- compare.fun(k0 = 0, d = 0.4)
# b1 <- compare.fun(k0 = 1, d = 0.4)
# b2 <- compare.fun(k0 = 2, d = 0.4)
# b3 <- compare.fun(k0 = 3, d = 0.4)
# b4 <- compare.fun(k0 = 4, d = 0.4)
# b5 <- compare.fun(k0 = 5, d = 0.4)

c0 <- compare.fun(k0 = 0, d = 0.6)
c1 <- compare.fun(k0 = 1, d = 0.6)
c2 <- compare.fun(k0 = 2, d = 0.6)
c3 <- compare.fun(k0 = 3, d = 0.6)
c4 <- compare.fun(k0 = 4, d = 0.6)
c5 <- compare.fun(k0 = 5, d = 0.6)

# save.image(
#   file = file.path(script_dir, "result-simu-1.RData")
# )
# 
# message(
#   "Simulation finished. Results saved to: ",
#   file.path(script_dir, "result-simu-1.RData")
# )
