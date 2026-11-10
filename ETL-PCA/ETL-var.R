library(ks)
library(parallel)

ETL_var <- function(
    X0,
    XS,
    A_grid = c(2.5, 3, 3.5),
    M_grid = c(0.2, 0.4,0.6,0.8,1,1.2),
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
  M_grid <- sort(unique(as.numeric(M_grid)))
  
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
  # Candidate (A, M) combinations
  #============================================================
  
  param_grid <- expand.grid(
    A = A_grid,
    M = M_grid
  )
  
  n_param <- nrow(param_grid)
  
  
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
    
    estimate_var <- function(
    x_target,
    A,
    M) {
      
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
        
        id <- which(
          abs(uk) <= A
        )
        
        if (length(id) == 0L) {
          next
        }
        
        u <- uk[id]
        
        
        #------------------------------------------------------
        # Density-ratio weight
        #------------------------------------------------------
        
        rho_hat <- sqrt(
          variance_ratio
        )
        
        denominator <- as.numeric(
          predict(
            fhat,
            x = u
          )
        )
        
        numerator <- rho_hat *
          as.numeric(
            predict(
              fhat,
              x = rho_hat * u
            )
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
      
      out <-
        (target_sum + source_sum) /
        (nt + source_weight_sum)
      
      if (!is.finite(out) || out <= 0) {
        return(NA_real_)
      }
      
      out
    }
    
    
    #==========================================================
    # Cross-validation over (A, M)
    #==========================================================
    
    cv_loss <- rep(
      Inf,
      n_param
    )
    
    for (g in seq_len(n_param)) {
      
      A <- param_grid$A[g]
      M <- param_grid$M[g]
      
      loss_sum <- 0
      n_val_sum <- 0
      
      
      for (f in seq_len(nfold)) {
        
        train_id <- folds[[f]]$train
        val_id   <- folds[[f]]$val
        
        x_train <- x0[train_id]
        x_val   <- x0[val_id]
        
        
        #------------------------------------------------------
        # Estimate variance using training target data
        #------------------------------------------------------
        
        v_hat <- estimate_var(
          x_target = x_train,
          A = A,
          M = M
        )
        
        if (!is.finite(v_hat) ||
            v_hat <= 0) {
          
          loss_sum <- Inf
          break
        }
        
        
        #------------------------------------------------------
        # Validation loss
        #------------------------------------------------------
        
        mu_train <- mean(
          x_train
        )
        
        loss_f <-
          log(v_hat) +
          (x_val - mu_train)^2 / v_hat
        
        loss_sum <- loss_sum +
          sum(loss_f)
        
        n_val_sum <- n_val_sum +
          length(x_val)
      }
      
      
      if (is.finite(loss_sum)) {
        
        cv_loss[g] <-
          loss_sum /
          n_val_sum
      }
    }
    
    
    #==========================================================
    # Full-data estimators for all (A, M)
    #==========================================================
    
    full_hat <- vapply(
      seq_len(n_param),
      function(g) {
        
        estimate_var(
          x_target = x0,
          A = param_grid$A[g],
          M = param_grid$M[g]
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
        "M_grid",
        "param_grid",
        "n_param",
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
      function(z) {
        z$cv_loss
      }
    )
  )
  
  
  cv_loss <- colMeans(
    cv_loss_coordinate
  )
  
  
  #============================================================
  # Select optimal (A, M)
  #============================================================
  
  best_id <- which.min(
    cv_loss
  )
  
  A_hat <- param_grid$A[
    best_id
  ]
  
  M_hat <- param_grid$M[
    best_id
  ]
  
  
  #============================================================
  # Final ETL-var
  #============================================================
  
  ETL_var_hat <- vapply(
    fit_list,
    function(z) {
      z$full_hat[best_id]
    },
    numeric(1)
  )
  
  
  #============================================================
  # CV loss matrix
  #
  # rows    = A
  # columns = M
  #============================================================
  
  cv_loss_matrix <- matrix(
    NA_real_,
    nrow = length(A_grid),
    ncol = length(M_grid),
    dimnames = list(
      paste0("A=", A_grid),
      paste0("M=", M_grid)
    )
  )
  
  for (g in seq_len(n_param)) {
    
    a_id <- match(
      param_grid$A[g],
      A_grid
    )
    
    m_id <- match(
      param_grid$M[g],
      M_grid
    )
    
    cv_loss_matrix[
      a_id,
      m_id
    ] <- cv_loss[g]
  }
  
  
  #============================================================
  # Output
  #============================================================
  
  list(
    ETL_var_hat = ETL_var_hat,
    A_hat = A_hat,
    M_hat = M_hat,
    cv_loss = cv_loss_matrix,
    A_grid = A_grid,
    M_grid = M_grid,
    cv_loss_coordinate = cv_loss_coordinate
  )
}
#============================================================
# Implementation
#============================================================

#==============================================
# data generation
#==============================================
ii = ii + 1
dat <- generate_sim3_data(d = 0.6, eta = 0, seed = ii)

X0 <- dat$X0
XS <- dat$XS
p = ncol(dat$X0)
K = length(XS)
var_true = diag(dat$Sigma0)

##================Implementation=========================
##=======================================================
fit <- ETL_var(
  X0 = X0,
  XS = XS
)

ETL_var_hat = fit$ETL_var_hat

X_all = X0
for (k in 1:K) {
  X_all = rbind(X_all,XS[[k]])
}


c(sum((ETL_var_hat - var_true)^2),
  sum((diag(cov(X0)) - var_true)^2),
  sum((diag(cov(X_all)) - var_true)^2))
fit$A_hat;fit$M_hat
