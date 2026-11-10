#============================================================
# 0. Packages
#============================================================

library(mvtnorm)
library(ks)
library(parallel)
library(glasso)
#============================================================
# 1. Source project files
#============================================================

source("data-generate.R")
source("ETL-functions.R")
source("ETL-var.R")
source("Trans-CLIME-functions.R")
  #========================================================
  # Data generation
  #========================================================
  
  dat <- generate_sim5_data(
    K = 10,
    k0 = 10,
    d = 0.6
  )
  
  
  X0 <- dat$X0
  XS <- dat$XS
  
  Omega0 <- dat$Omega0
  OmegaS <- dat$Omega0
  
  #========================================================
  # 8. Target-Glasso
  #========================================================
  
  fit = fit.glasso.cv(X0,rho_grid = seq(0.1, 0.3, length.out = 10))
  Target_Glasso_hat = fit$Omega_hat
  
  #========================================================
  # 11. Pool-Glasso
  #========================================================
  
  fit = fit.pool.glasso.cv(
    X0,
    XS,
    rho_grid = seq(0.05, 0.3, length.out = 10)
    )
  
  Pool_Glasso_hat <- fit$Omega_hat

  
  #========================================================
  # 13. Trans-CLIME
  #
  # Standard Trans-CLIME:
  # all source samples are used; k0 is NOT supplied.
  #========================================================
  
  Trans_CLIME_hat <-
    fit.Trans.CLIME(
      X0 = X0,
      XS = XS,
      const = 0.5
    )
  
  #========================================================
  # 3. ETL-var
  #========================================================
  
  fit_var <- ETL_var(
    X0 = X0,
    XS = XS,
    A_grid = c(3,3.5,4,4.5,5),
    M_grid = c(0.2,0.4,0.6,0.8),
    ncores = 8
  )


  ETL_var_hat <- fit_var$ETL_var_hat
  # ETL_var_hat = diag(solve(dat$Omega0))
  
  #========================================================
  # 1. ETL-corr
  #========================================================
  
  fit_corr <- fit.ETL.corr(
    X0 = X0,
    XS = XS,
    c_A1 = 0.8,
    M = 0.3,
    ncores = 8
  )

  ETL_corr_hat <- fit_corr$R_hat
  
  #========================================================
  # 4. ETL-cov
  #========================================================
  
  D_ETL <-
    diag(
      sqrt(
        ETL_var_hat
      )
    )
  
  
  ETL_cov_hat <-
    D_ETL %*%
    ETL_corr_hat %*%
    D_ETL
  
  
  ETL_cov_hat <-
    make_pd(
      ETL_cov_hat,
      min_eig = 0.1
    )

  #========================================================
  # 9. ETL-Glasso
  #========================================================
  rho.range = seq(0.01,0.2,length.out = 10)
  loss.ETL = rep(0,length(rho.range))
  
  for (r in seq_along(rho.range)) {
    
    ETL_Glasso_hat <- glasso(s = ETL_cov_hat,rho = rho.range[r],penalize.diagonal = FALSE)$wi
    loss.ETL[r] <- loss.fun(Omega_hat = ETL_Glasso_hat, Omega_true = Omega0)[1]
  }
  rho.opt = rho.range[which.min(loss.ETL)];rho.opt
  ETL_Glasso_hat <- glasso(s = ETL_cov_hat,rho = 0.005,penalize.diagonal = FALSE)$wi
  
  
  
  #========================================================
  # Losses
  #========================================================
  loss_ETL_Glasso <- loss.fun(
    Omega_hat = ETL_Glasso_hat,
    Omega_true = Omega0
  )
  
 loss_Target_Glasso <- loss.fun(
    Omega_hat = Target_Glasso_hat,
    Omega_true = Omega0
  )
  
  loss_Pool_Glasso <- loss.fun(
    Omega_hat = Pool_Glasso_hat,
    Omega_true = Omega0
  )
  
  loss_Trans_CLIME <- loss.fun(
    Omega_hat = Trans_CLIME_hat,
    Omega_true = Omega0
  )
  
  result <- rbind(
    "ETL-Glasso"     = loss_ETL_Glasso,
    "Target-Glasso"  = loss_Target_Glasso,
    "Pool-Glasso"    = loss_Pool_Glasso,
    "Trans-CLIME"    = loss_Trans_CLIME
  )
  
  result
  