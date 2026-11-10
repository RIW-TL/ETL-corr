library(ks)
library(parallel)

ETL_var <- function(
    X0,
    XS,
    A_grid = c(2, 2.5, 3),
    M = 0.5,
    nfold = 5,
    ncores = 8,
    seed = 1,
    cl = NULL
) {
  
  #============================================================
  # Basic setup
  #============================================================
  
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Please install the 'ks' package.")
  }
  
  X0 <- as.matrix(X0)
  XS <- lapply(XS, as.matrix)
  
  n0 <- nrow(X0)
  p  <- ncol(X0)
  
  A_grid <- sort(unique(as.numeric(A_grid)))
  
  if (length(XS) == 0L) {
    stop("XS must be a nonempty list.")
  }
  
  if (any(vapply(XS, ncol, integer(1)) != p)) {
    stop("X0 and all source matrices must have the same number of columns.")
  }
  
  if (nfold < 2L || nfold >= n0) {
    stop("nfold must satisfy 2 <= nfold < nrow(X0).")
  }
  
  
  #============================================================
  # Five-fold split
  #============================================================
  
  set.seed(seed)
  
  fold_id <- sample(
    rep(seq_len(nfold), length.out = n0)
  )
  
  folds <- lapply(
    seq_len(nfold),
    function(f) {
      list(
        train = which(fold_id != f),
        val   = which(fold_id == f)
      )
    }
  )
  
  
  #============================================================
  # One coordinate
  #============================================================
  
  fit_coordinate <- function(j) {
    
    x0 <- X0[, j]
    
    
    #----------------------------------------------------------
    # Source quantities
    #----------------------------------------------------------
    
    source_info <- lapply(
      XS,
      function(Xk) {
        
        xk <- Xk[, j]
        uk <- xk - mean(xk)
        
        list(
          uk   = uk,
          vk   = var(xk),
          fhat = ks::kde(x = uk)
        )
      }
    )
    
    
    #==========================================================
    # ETL-var estimator
    #==========================================================
    
    estimate_var <- function(x_target, A) {
      
      nt <- length(x_target)
      v0 <- var(x_target)
      
      if (!is.finite(v0) || v0 <= 0) {
        return(NA_real_)
      }
      
      target_sum <- nt * v0
      
      source_sum        <- 0
      source_weight_sum <- 0
      
      
      for (k in seq_along(source_info)) {
        
        uk   <- source_info[[k]]$uk
        vk   <- source_info[[k]]$vk
        fhat <- source_info[[k]]$fhat
        
        
        #------------------------------------------------------
        # Source-level screening
        #------------------------------------------------------
        
        variance_ratio <- vk / v0
        
        if (!is.finite(variance_ratio) ||
            abs(variance_ratio - 1) > M) {
          next
        }
        
        
        #------------------------------------------------------
        # Residual truncation
        #------------------------------------------------------
        
        id <- which(abs(uk) <= A)
        
        if (length(id) == 0L) {
          next
        }
        
        u <- uk[id]
        
        
        #------------------------------------------------------
        # Density-ratio weight
        #------------------------------------------------------
        
        rho_hat <- sqrt(variance_ratio)
        
        denominator <- as.numeric(
          predict(fhat, x = u)
        )
        
        numerator <- rho_hat * as.numeric(
          predict(fhat, x = rho_hat * u)
        )
        
        weight <- numerator / denominator
        
        valid <- is.finite(weight)
        
        if (!any(valid)) {
          next
        }
        
        u      <- u[valid]
        weight <- weight[valid]
        
        
        #------------------------------------------------------
        # Source contribution
        #------------------------------------------------------
        
        source_sum <- source_sum +
          sum(weight * u^2)
        
        source_weight_sum <- source_weight_sum +
          sum(weight)
      }
      
      
      #--------------------------------------------------------
      # ETL variance
      #--------------------------------------------------------
      
      out <- (target_sum + source_sum) /
        (nt + source_weight_sum)
      
      if (!is.finite(out) || out <= 0) {
        return(NA_real_)
      }
      
      out
    }
    
    
    #==========================================================
    # Cross-validation loss
    #
    # Gaussian quasi-likelihood:
    #
    # log(v_hat) + (X - mu_train)^2 / v_hat
    #==========================================================
    
    cv_loss <- rep(0, length(A_grid))
    
    for (a in seq_along(A_grid)) {
      
      A <- A_grid[a]
      
      loss_sum <- 0
      n_val_sum <- 0
      
      for (f in seq_len(nfold)) {
        
        train_id <- folds[[f]]$train
        val_id   <- folds[[f]]$val
        
        x_train <- x0[train_id]
        x_val   <- x0[val_id]
        
        v_hat <- estimate_var(
          x_target = x_train,
          A = A
        )
        
        if (!is.finite(v_hat) || v_hat <= 0) {
          loss_sum <- Inf
          break
        }
        
        mu_train <- mean(x_train)
        
        loss_f <- log(v_hat) +
          (x_val - mu_train)^2 / v_hat
        
        loss_sum <- loss_sum + sum(loss_f)
        n_val_sum <- n_val_sum + length(x_val)
      }
      
      if (is.finite(loss_sum)) {
        cv_loss[a] <- loss_sum / n_val_sum
      } else {
        cv_loss[a] <- Inf
      }
    }
    
    
    #==========================================================
    # Full-data estimators for all A
    #==========================================================
    
    full_hat <- vapply(
      A_grid,
      function(A) {
        estimate_var(
          x_target = x0,
          A = A
        )
      },
      numeric(1)
    )
    
    
    list(
      cv_loss = cv_loss,
      full_hat = full_hat
    )
  }
  
  
  #============================================================
  # Parallel cluster
  #============================================================
  
  own_cluster <- FALSE
  
  if (is.null(cl) && ncores > 1L) {
    
    cl <- parallel::makeCluster(
      min(ncores, p)
    )
    
    master_lib <- .libPaths()
    
    parallel::clusterCall(
      cl,
      function(lib) {
        .libPaths(c(lib, .libPaths()))
        invisible(NULL)
      },
      master_lib
    )
    
    own_cluster <- TRUE
  }
  
  on.exit(
    {
      if (own_cluster && !is.null(cl)) {
        parallel::stopCluster(cl)
      }
    },
    add = TRUE
  )
  
  
  #============================================================
  # Parallel over coordinates
  #============================================================
  
  if (is.null(cl)) {
    
    fit_list <- lapply(
      seq_len(p),
      fit_coordinate
    )
    
  } else {
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "X0",
        "XS",
        "A_grid",
        "M",
        "nfold",
        "folds",
        "fit_coordinate"
      ),
      envir = environment()
    )
    
    fit_list <- parallel::parLapply(
      cl,
      seq_len(p),
      fit_coordinate
    )
  }
  
  
  #============================================================
  # Global CV criterion
  #============================================================
  
  cv_loss_coordinate <- do.call(
    rbind,
    lapply(
      fit_list,
      function(z) z$cv_loss
    )
  )
  
  cv_loss <- colMeans(
    cv_loss_coordinate
  )
  
  names(cv_loss) <- paste0(
    "A=", A_grid
  )
  
  best_id <- which.min(cv_loss)
  A_hat   <- A_grid[best_id]
  
  
  #============================================================
  # Final ETL-var
  #============================================================
  
  ETL_var_hat <- vapply(
    fit_list,
    function(z) z$full_hat[best_id],
    numeric(1)
  )
  
  
  list(
    ETL_var_hat = ETL_var_hat,
    A_hat = A_hat,
    cv_loss = cv_loss,
    A_grid = A_grid,
    cv_loss_coordinate = cv_loss_coordinate
  )
}

#============================================================
# Implementation
#============================================================

compare_var = function(eta = 0.9, M = 0.5, iter = 2){
  
  result = matrix(0,iter,4)
  colnames(result) = c("ETL","Target","ADT","pool")
  
  for (ii in 1:iter) {
    
    #==============================================
    # data generation
    #==============================================
    dat <- generate_sim2_data(eta = eta, tau = 0.1, seed = ii)
    
    X0 <- dat$X0
    XS <- dat$XS
    
    var_true = diag(dat$Sigma0)
    
    ##================Implementation=========================
    ##=======================================================
    fit <- ETL_var(
      X0 = X0,
      XS = XS,
      M = M
    )
    ETL_var_hat = fit$ETL_var_hat
    
    
    ADT.sigma.hat = ADT(X0,XS[[1]])$estimate
    pool_var_hat = diag(cov(rbind(X0,XS[[1]])))
    
    
    result[ii,] = c(
      sum((ETL_var_hat - var_true)^2),
      sum((diag(cov(X0)) - var_true)^2),
      sum((diag(ADT.sigma.hat) - var_true)^2),
      sum((pool_var_hat - var_true)^2))
  }
  
  return(colMeans(result))
}

v1 = compare_var(eta = 0.1)
v2 = compare_var(eta = 0.3)
v3 = compare_var(eta = 0.5)
v4 = compare_var(eta = 0.7)
v5 = compare_var(eta = 0.9)
rbind(v1,v2,v3,v4,v5)

#=====================Comparison=======================
#======================================================
# result = c(sum(abs(ETL_var_hat - 1))/p,
#            sum(abs(target_var_hat - 1))/p,
#            sum(abs(pool_var_hat - 1))/p)
# names(result) = c("ETL-var","Target-var","Pool-var");result

