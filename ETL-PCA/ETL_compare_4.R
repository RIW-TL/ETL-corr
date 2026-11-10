library(ks)
library(mvtnorm)
library(sparseCov)
library(parallel)

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

cov_loss_fun <- function(Sigma_hat, U_t,r) {
  
  U_hat = eigen(Sigma_hat)$vectors[,1:r]
  
  index = norm(U_hat %*% t(U_hat)  - U_t %*% t(U_t),"F")
  return(index)
}


#============================================================
# Main comparison function
#============================================================

compare.fun <- function(k0,
                        d = 6,
                        c_A1 = 0.5,
                        c_M1 = 0.05,
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
    "Target-cov",
    "Target-Scov",
    "Pool-cov",
    "GB-1",
    "GB-2",
    "GB-3"
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
  
  #==========================================================
  # Monte Carlo replications
  #==========================================================
  
  for (ii in seq_len(iter)) {
    
    #========================================================
    # Data generation
    #========================================================
    
    dat <- generate_sim4_data(k0 = k0, d = d, seed = ii)
    
    X0 <- dat$X0
    XS <- dat$XS
    p = ncol(dat$X0)
    K = length(XS)
    r0 = dat$r0
    U_t = dat$U_t
    rk_list = dat$rk_list
    
    Sigma_true = dat$Sig_t
    
    
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
      c_T = 1
    )
    
    
    #========================================================
    # 3. ETL-var
    #========================================================
    
    # fit <- ETL_var(
    #   X0 = X0,
    #   XS = XS
    # )
    # 
    # ETL_var_hat = fit$ETL_var_hat
    ETL_var_hat = diag(Sigma_true)
    
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
    # GB-PCA
    #========================================================
    
    GB_hat_1 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.25)
    GB_hat_2 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.5)
    GB_hat_3 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.75)
    
    #========================================================
    # Compute losses
    #========================================================
    
    loss_ETL_cov <- cov_loss_fun(
      Sigma_hat = ETL_cov_hat,
      U_t = U_t,
      r = r0
    )
    
    loss_ETL_scov <- cov_loss_fun(
      Sigma_hat = ETL_scov_hat,
      U_t = U_t,
      r = r0
    )
    
    loss_Target_cov <- cov_loss_fun(
      Sigma_hat = Target_cov_hat,
      U_t = U_t,
      r = r0
    )
    
    loss_Target_scov <- cov_loss_fun(
      Sigma_hat = Target_scov_hat,
      U_t = U_t,
      r = r0
    )
    
    loss_Pool_cov <- cov_loss_fun(
      Sigma_hat = Pool_cov_hat,
      U_t = U_t,
      r = r0
    )
    
    loss_GB1_cov <- norm(GB_hat_1  - U_t %*% t(U_t),"F")
    loss_GB2_cov <- norm(GB_hat_2  - U_t %*% t(U_t),"F")
    loss_GB3_cov <- norm(GB_hat_3  - U_t %*% t(U_t),"F")
    #========================================================
    # Store covariance losses
    #========================================================
    
    result.cov[ii, ] <- c(
      loss_ETL_cov,
      loss_ETL_scov,
      loss_Target_cov,
      loss_Target_scov,
      loss_Pool_cov,
      loss_GB1_cov,
      loss_GB2_cov,
      loss_GB3_cov
    )
  }
  
  
  #==========================================================
  # Mean losses
  #==========================================================
  
  mean.cov <- colMeans(
    result.cov,
    na.rm = TRUE
  )
  
  
  #==========================================================
  # Output
  #==========================================================
  
  list(
    cov_loss = result.cov,
    mean_cov_loss = mean.cov
  )
}


#============================================================
# Simulation 3
#============================================================


# a0 <- compare.fun(k0 = 0,d = 4)
# a1 <- compare.fun(k0 = 1,d = 4)
a2 <- compare.fun(k0 = 2,d = 4)
# a3 <- compare.fun(k0 = 3,d = 4)
a4 <- compare.fun(k0 = 4,d = 4)
# a5 <- compare.fun(k0 = 5,d = 4)

# b0 <- compare.fun(k0 = 0,d = 6)
# b1 <- compare.fun(k0 = 1,d = 6)
b2 <- compare.fun(k0 = 2,d = 6)
# b3 <- compare.fun(k0 = 3,d = 6)
b4 <- compare.fun(k0 = 4,d = 6)
# b5 <- compare.fun(k0 = 5,d = 6)

# c0 <- compare.fun(k0 = 0,d = 8)
# c1 <- compare.fun(k0 = 1,d = 8)
c2 <- compare.fun(k0 = 2,d = 8)
# c3 <- compare.fun(k0 = 3,d = 8)
c4 <- compare.fun(k0 = 4,d = 8)
# c5 <- compare.fun(k0 = 5,d = 8)

# d0 <- compare.fun(k0 = 0,d = 20)
# d1 <- compare.fun(k0 = 1,d = 20)
d2 <- compare.fun(k0 = 2,d = 20)
# d3 <- compare.fun(k0 = 3,d = 20)
d4 <- compare.fun(k0 = 4,d = 20)
# d5 <- compare.fun(k0 = 5,d = 20)

rbind(a0$mean_cov_loss,a1$mean_cov_loss,a3$mean_cov_loss,a5$mean_cov_loss)
rbind(b0$mean_cov_loss,b1$mean_cov_loss,b3$mean_cov_loss,b5$mean_cov_loss)
rbind(c0$mean_cov_loss,c1$mean_cov_loss,c3$mean_cov_loss,c5$mean_cov_loss)
rbind(d0$mean_cov_loss,d1$mean_cov_loss,d3$mean_cov_loss,d5$mean_cov_loss)
