#============================================================
# Simulation 5: ETL-Glasso comparison

#============================================================
# 0. Packages
#============================================================

library(mvtnorm)
library(ks)
library(parallel)
library(glasso)
library(reticulate)


#============================================================
# 1. Source project files
#============================================================

source("data-generate.R")
source("ETL-functions.R")
source("ETL-var.R")
source("Trans-CLIME-functions.R")


#============================================================
# 2. Positive-definite adjustment
#============================================================

make_pd <- function(S, min_eig = 0.1) {
  
  S <- as.matrix(S)
  S <- (S + t(S)) / 2
  
  lambda_min <- min(
    eigen(
      S,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  if (lambda_min < min_eig) {
    S <- S +
      (min_eig - lambda_min) *
      diag(nrow(S))
  }
  
  S
}

#============================================================
# 4. Target Glasso
#============================================================

fit.glasso.cv <- function(
    X0,
    rho_grid = seq(0.05, 0.5, length.out = 20),
    nfold = 5,
    seed = 1
) {
  
  X0 <- as.matrix(X0)
  n <- nrow(X0)
  
  #========================================================
  # 1. Cross-validation
  #========================================================
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  
  nrho <- length(rho_grid)
  cv_loss <- matrix(Inf, nfold, nrho)
  
  for (k in seq_len(nfold)) {
    
    X_train <- X0[fold_id != k, , drop = FALSE]
    X_val   <- X0[fold_id == k, , drop = FALSE]
    
    mu_train <- colMeans(X_train)
    
    X_train_c <- sweep(X_train, 2, mu_train, "-")
    X_val_c   <- sweep(X_val, 2, mu_train, "-")
    
    S_train <- crossprod(X_train_c) / nrow(X_train)
    S_val   <- crossprod(X_val_c) / nrow(X_val)
    
    fit <- glasso::glassopath(
      s = S_train,
      rholist = rho_grid,
      penalize.diagonal = FALSE,
      trace = 0
    )
    
    for (j in seq_len(nrho)) {
      
      if (fit$errflag[j] != 0) next
      
      Omega <- fit$wi[, , j]
      
      R <- tryCatch(
        chol(Omega),
        error = function(e) NULL
      )
      
      if (!is.null(R)) {
        cv_loss[k, j] <-
          sum(S_val * Omega) -
          2 * sum(log(diag(R)))
      }
    }
  }
  
  #========================================================
  # 2. Optimal penalty
  #========================================================
  
  mean_loss <- colMeans(cv_loss)
  
  if (all(!is.finite(mean_loss))) {
    stop("All cross-validation fits failed.")
  }
  
  rho_opt <- rho_grid[which.min(mean_loss)]
  
  #========================================================
  # 3. Final estimation
  #========================================================
  
  Xc <- sweep(X0, 2, colMeans(X0), "-")
  S0 <- crossprod(Xc) / n
  
  fit_final <- glasso::glasso(
    s = S0,
    rho = rho_opt,
    penalize.diagonal = FALSE
  )
  
  return(list(
    Omega_hat = fit_final$wi,
    Sigma_hat = fit_final$w,
    rho_opt = rho_opt,
    rho_grid = rho_grid,
    cv_loss = mean_loss
  ))
}

fit.pool.glasso.cv <- function(
    X0,
    XS,
    rho_grid = seq(0.05, 0.5, length.out = 20),
    nfold = 5,
    seed = 1
) {
  
  X0 <- as.matrix(X0)
  
  if (!is.list(XS)) XS <- list(XS)
  XS <- lapply(XS, as.matrix)
  
  n <- nrow(X0)
  p <- ncol(X0)
  
  stopifnot(
    nfold >= 2,
    nfold < n,
    length(XS) >= 1,
    all(is.finite(X0)),
    all(vapply(XS, ncol, integer(1)) == p),
    all(vapply(XS, function(x) all(is.finite(x)), logical(1))),
    all(is.finite(rho_grid)),
    all(rho_grid > 0),
    all(diff(rho_grid) > 0)
  )
  
  X_source <- do.call(rbind, XS)
  
  #========================================================
  # 1. Cross-validation
  #========================================================
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  
  nrho <- length(rho_grid)
  cv_loss <- matrix(Inf, nfold, nrho)
  
  for (k in seq_len(nfold)) {
    
    X_train <- X0[fold_id != k, , drop = FALSE]
    X_val   <- X0[fold_id == k, , drop = FALSE]
    
    # Pool target training data and all sources
    X_pool <- rbind(X_train, X_source)
    
    S_train <- cov(X_pool)
    
    mu_train <- colMeans(X_pool)
    X_val_c <- sweep(X_val, 2, mu_train, "-")
    
    S_val <- crossprod(X_val_c) / nrow(X_val)
    
    fit <- glasso::glassopath(
      s = S_train,
      rholist = rho_grid,
      penalize.diagonal = FALSE,
      trace = 0
    )
    
    for (j in seq_len(nrho)) {
      
      if (fit$errflag[j] != 0) next
      
      Omega <- fit$wi[, , j]
      
      R <- tryCatch(
        chol(Omega),
        error = function(e) NULL
      )
      
      if (!is.null(R)) {
        cv_loss[k, j] <-
          sum(S_val * Omega) -
          2 * sum(log(diag(R)))
      }
    }
  }
  
  #========================================================
  # 2. Optimal penalty
  #========================================================
  
  mean_loss <- colMeans(cv_loss)
  
  if (all(!is.finite(mean_loss))) {
    stop("All cross-validation fits failed.")
  }
  
  rho_opt <- rho_grid[which.min(mean_loss)]
  
  #========================================================
  # 3. Final Pool-Glasso estimation
  #========================================================
  
  X_all <- rbind(X0, X_source)
  S_all <- cov(X_all)
  
  fit_final <- glasso::glasso(
    s = S_all,
    rho = rho_opt,
    penalize.diagonal = FALSE
  )
  
  return(list(
    Omega_hat = fit_final$wi,
    Sigma_hat = fit_final$w,
    rho_opt = rho_opt,
    rho_grid = rho_grid,
    cv_loss = mean_loss
  ))
}


#============================================================
# 5. Loss function
#============================================================

loss.fun <- function(
    Omega_hat,
    Omega_true) {
  
  Omega_hat <- as.matrix(Omega_hat)
  Omega_true <- as.matrix(Omega_true)
  
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2
  
  p <- ncol(Omega_true)
  
  Sigma_true <- solve(Omega_true)
  
  
  #----------------------------------------------------------
  # Precision estimation loss
  #----------------------------------------------------------
  
  precision_loss <-
    norm(
      Omega_hat - Omega_true,
      type = "F"
    )^2 / p
  
  
  #----------------------------------------------------------
  # Negative log-likelihood loss
  #----------------------------------------------------------
  
  det_Omega <- determinant(
    Omega_hat,
    logarithm = TRUE
  )
  
  if (det_Omega$sign <= 0) {
    
    nll_loss <- Inf
    
  } else {
    
    logdet <- as.numeric(
      det_Omega$modulus
    )
    
    nll_loss <-
      (
        sum(
          ETL_cov_hat * Omega_hat
        ) -
          logdet
      ) / p
  }
  
  
  c(
    precision = precision_loss,
    nll = nll_loss
  )
}


#============================================================
# 6. Main comparison function
#============================================================

compare.fun <- function(k0,d,iter = 10,ncores = 8) {
  
  
  #==========================================================
  # Method names
  #==========================================================
  
  methods <- c(
    "ETL-Glasso",
    "Target-Glasso",
    "Pool-Glasso",
    "Trans-CLIME"
  )
  
  
  #==========================================================
  # Containers
  #==========================================================
  
  result.precision <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(methods),
    dimnames = list(
      paste0(
        "iter",
        seq_len(iter)
      ),
      methods
    )
  )
  
  
  result.nll <- matrix(
    NA_real_,
    nrow = iter,
    ncol = length(methods),
    dimnames = list(
      paste0(
        "iter",
        seq_len(iter)
      ),
      methods
    )
  )
  
  
  #==========================================================
  # Persistent parallel cluster for ETL
  #==========================================================
  
  cl <- NULL
  
  if (ncores > 1L) {
    
    cl <- parallel::makeCluster(
      ncores
    )
    
    
    master_lib <- .libPaths()
    
    
    parallel::clusterCall(
      cl,
      function(lib) {
        
        .libPaths(
          unique(
            c(
              lib,
              .libPaths()
            )
          )
        )
        
        invisible(NULL)
      },
      master_lib
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
  # Monte Carlo replications
  #==========================================================
  
  for (ii in seq_len(iter)) {
    
    cat(
      "k0 =", k0,
      " d =", d,
      " iteration =", ii,
      "/", iter,
      "\n"
    )
    
    
    #========================================================
    # Data generation
    #========================================================
    
    dat <- generate_sim5_data(
      K = 10,
      k0 = k0,
      d = d,
      seed = ii
    )
    
    
    X0 <- dat$X0
    XS <- dat$XS
    
    Omega0 <- dat$Omega0
    OmegaS <- dat$OmegaS
    
    #========================================================
    # 1. ETL-corr
    #========================================================
    
    fit_corr <- fit.ETL.corr(
      X0 = X0,
      XS = XS,
      ncores = ncores,
      cl = cl
    )
    
    
    ETL_corr_hat <-
      fit_corr$R_hat
    
    #========================================================
    # 3. ETL-var
    #========================================================
    
    fit_var <- ETL_var(
      X0 = X0,
      XS = XS,
      A_grid = c(3,3.5,4,4.5,5),
      M_grid = c(0.2,0.4,0.6,0.8),
      ncores = ncores,
      cl = cl
    )
  
    ETL_var_hat <- fit_var$ETL_var_hat
    
    
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
        ETL_cov_hat
      )
    
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
    # 9. ETL-Glasso
    #========================================================
    # rho.range = seq(0.01,0.2,length.out = 10)
    # loss.ETL = rep(0,length(rho.range))
    # 
    # for (r in seq_along(rho.range)) {
    #   
    #   ETL_Glasso_hat <- glasso(s = ETL_cov_hat,rho = rho.range[r])$wi
    #   loss.ETL[r] <- loss.fun(Omega_hat = ETL_Glasso_hat, Omega_true = Omega0)[2]
    # }
    # rho.opt = rho.range[which.min(loss.ETL)]
    ETL_Glasso_hat <- glasso(s = ETL_cov_hat,rho = 0.03,penalize.diagonal = FALSE)$wi
    
    #========================================================
    # 13. Trans-CLIME
    #========================================================
    
    Trans_CLIME_hat <-
      fit.Trans.CLIME(
        X0 = X0,
        XS = XS,
        const = 0.5
      )
    
    
    #========================================================
    # Collect all estimators
    #========================================================
    
    Omega_hat_list <- list(
      "ETL-Glasso" =
        ETL_Glasso_hat,
      
      "Target-Glasso" =
        Target_Glasso_hat,
      
      "Pool-Glasso" =
        Pool_Glasso_hat,
      
      "Trans-CLIME" =
        Trans_CLIME_hat
    )
    
    
    #========================================================
    # Losses
    #========================================================
    
    for (m in seq_along(methods)) {
      
      loss <-
        loss.fun(
          Omega_hat = Omega_hat_list[[methods[m]]],
          Omega_true = Omega0
        )
      
      
      result.precision[
        ii,
        m
      ] <-
        loss["precision"]
      
      
      result.nll[
        ii,
        m
      ] <-
        loss["nll"]
    }
  }
  
  
  #==========================================================
  # Mean losses over all iter replications
  #==========================================================
  
  mean.precision <-
    colMeans(
      result.precision
    )
  
  
  mean.nll <-
    colMeans(
      result.nll
    )
  
  
  #==========================================================
  # Output:
  #==========================================================
  
  list(
    
    precision_loss =
      result.precision,
    
    nll_loss =
      result.nll,
    
    mean_precision_loss =
      mean.precision,
    
    mean_nll_loss =
      mean.nll
  )
}


#============================================================
# 8. Example
#============================================================

a1 <- compare.fun(k0 = 2, d = 0.2)
a2 <- compare.fun(k0 = 4, d = 0.2)
a3 <- compare.fun(k0 = 6, d = 0.2)
a4 <- compare.fun(k0 = 8, d = 0.2)
a5 <- compare.fun(k0 = 10, d = 0.2)

b1 <- compare.fun(k0 = 2, d = 0.3)
b2 <- compare.fun(k0 = 4, d = 0.3)
b3 <- compare.fun(k0 = 6, d = 0.3)
b4 <- compare.fun(k0 = 8, d = 0.3)
b5 <- compare.fun(k0 = 10, d = 0.3)

c1 <- compare.fun(k0 = 2, d = 0.4)
c2 <- compare.fun(k0 = 4, d = 0.4)
c3 <- compare.fun(k0 = 6, d = 0.4)
c4 <- compare.fun(k0 = 8, d = 0.4)
c5 <- compare.fun(k0 = 10, d = 0.4)

d1 <- compare.fun(k0 = 2, d = 0.5)
d2 <- compare.fun(k0 = 4, d = 0.5)
d3 <- compare.fun(k0 = 6, d = 0.5)
d4 <- compare.fun(k0 = 8, d = 0.5)
d5 <- compare.fun(k0 = 10, d = 0.5)

lapply(c("a", "b", "c", "d"), function(x) {
  rbind(
    get(paste0(x, "1"))$mean_precision_loss,
    get(paste0(x, "2"))$mean_precision_loss,
    get(paste0(x, "3"))$mean_precision_loss,
    get(paste0(x, "4"))$mean_precision_loss,
    get(paste0(x, "5"))$mean_precision_loss
  )
})