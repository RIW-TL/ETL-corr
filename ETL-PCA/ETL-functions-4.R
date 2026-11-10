library(ks)
library(parallel)
library(sparseCov)

#============================================================
# 1. ETL-corr for a grid of M. Parallelization only over (j,l)
#============================================================

fit.ETL.corr <- function(
    X0,
    XS,
    c_A1 = 0.9,
    c_M1 = 0.05,
    c_M2 = 0.1,
    c_A2 = 0.95,
    c_H = 1,
    density_floor = 1e-10,
    ncores = 8,
    cl = NULL) {
  
  p <- ncol(X0)
  K <- length(XS)
  
  #==========================================================
  # Initial estimators
  #==========================================================
  
  # R_t_initial <- fit.target.scor(X0,corr = TRUE)
  # R_s_initial <- list()
  # for (k in 1:K) {
  #   R_s_initial[[k]] = fit.target.scor(Xs[[k]],corr = TRUE)
  # }
  # R_t_initial <- cov2cor(Sigma_true)
  R_t_initial <- cor(X0)
  R_s_initial <- list()
  for (k in 1:K) {
    R_s_initial[[k]] = cor(XS[[k]])
  }
  
  #==========================================================
  # Standardization
  #==========================================================
  
  Y0 <- scale(X0)
  YS <- lapply(XS, scale)
  
  pair_id <- which(
    lower.tri(matrix(0, p, p)),
    arr.ind = TRUE
  )
  
  target_cross <- crossprod(Y0)
  
  
  #==========================================================
  # Precompute A2
  #
  # A2 depends only on source k and coordinate l,
  # not on j, so it should NOT be recomputed for every pair.
  #==========================================================
  
  A2_list <- lapply(
    YS,
    function(Yk) {
      apply(
        abs(Yk),
        2,
        quantile,
        probs = c_A2,
        names = FALSE
      )
    }
  )
  
  
  #==========================================================
  # One pair
  #==========================================================
  
  fit_one_pair <- function(ii) {
    
    j <- pair_id[ii, 1]
    l <- pair_id[ii, 2]
    
    #--------------------------------------------------------
    # Target contribution
    #--------------------------------------------------------
    
    num <- target_cross[l, j]
    den <- target_cross[l, l]
    
    
    #========================================================
    # Sources
    #========================================================
    
    for (k in seq_len(K)) {
      
      yj <- YS[[k]][, j]
      yl <- YS[[k]][, l]
      
      nk <- length(yl)
      
      rho_t <- R_t_initial[j, l]
      rho_s <- R_s_initial[[k]][j, l]
      
      
      #------------------------------------------------------
      # Residuals
      #------------------------------------------------------
      
      error_1 <- yj - rho_s * yl
      error_0 <- yj - rho_t * yl
      
      delta_abs <- abs(
        (rho_s - rho_t) * yl
      )
      
      
      #------------------------------------------------------
      # Sample selection
      #------------------------------------------------------
      
      A1 <- as.numeric(
        quantile(
          abs(error_1),
          probs = c_A1,
          names = FALSE
        )
      )
      
      # Precomputed
      A2 <- A2_list[[k]][l]
      
      # Adaptive M
      M <- c_M1 /
        (c_M2 + abs(rho_s - rho_t))
      
      id <-
        abs(error_1) <= A1 &
        abs(yl) <= A2 &
        delta_abs <= M
      
      
      #------------------------------------------------------
      # Skip KDE if nothing is selected
      #------------------------------------------------------
      
      if (!any(id))
        next
      
      id_sel <- which(id)
      m <- length(id_sel)
      
      
      #------------------------------------------------------
      # KDE training sample
      #------------------------------------------------------
      
      Z <- cbind(
        error_1,
        yl
      )
      
      
      #------------------------------------------------------
      # Bandwidth matrix
      #
      # Because Y_j and Y_l are standardized and rho_s is
      # their sample correlation:
      #
      # Var(error_1) = 1 - rho_s^2
      # Cov(error_1, yl) = 0
      # Var(yl) = 1
      #
      # Therefore cov(Z) need not be recomputed.
      #------------------------------------------------------
      
      H <-
        c_H *
        nk^(-1 / 3) *
        diag(
          c(
            max(1 - rho_s^2, 1e-8),
            1
          )
        )
      
      H <- H + 1e-8 * diag(2)
      
      
      #------------------------------------------------------
      # Evaluate KDE ONLY at selected observations
      #------------------------------------------------------
      
      eval_points <- rbind(
        cbind(
          error_1[id_sel],
          yl[id_sel]
        ),
        cbind(
          error_0[id_sel],
          yl[id_sel]
        ),
        cbind(
          -error_0[id_sel],
          yl[id_sel]
        )
      )
      
      kde_fit <- ks::kde(
        x = Z,
        H = H,
        eval.points = eval_points,
        binned = TRUE,
        bgridsize = c(41, 41),
        compute.cont = FALSE
      )
      
      f_all <- as.numeric(
        kde_fit$estimate
      )
      
      ii_sel <- seq_len(m)
      
      #------------------------------------------------------
      # Density ratio
      #------------------------------------------------------
      
      f_den <- f_all[ii_sel]
      
      f_num1 <- f_all[
        m + ii_sel
      ]
      
      f_num2 <- f_all[
        2 * m + ii_sel
      ]
      
      f_den <- pmax(
        f_den,
        density_floor
      )
      
      weight <-
        (f_num1 + f_num2) /
        (2 * f_den)
      
      
      #------------------------------------------------------
      # Weighted source contribution
      #------------------------------------------------------
      
      num <- num +
        sum(
          weight *
            yj[id_sel] *
            yl[id_sel]
        )
      
      den <- den +
        sum(
          weight *
            yl[id_sel]^2
        )
    }
    
    
    list(
      rho = num / den,
      H = den
    )
  }
  
  
  #==========================================================
  # Parallel cluster
  #==========================================================
  
  own_cluster <- FALSE
  
  if (is.null(cl) && ncores > 1L) {
    
    cl <- parallel::makeCluster(
      ncores
    )
    
    user_lib <- path.expand(
      "~/R/library"
    )
    
    parallel::clusterCall(
      cl,
      function(lib) {
        .libPaths(
          c(lib, .libPaths())
        )
        invisible(NULL)
      },
      user_lib
    )
    
    own_cluster <- TRUE
  }
  
  on.exit(
    {
      if (own_cluster &&
          !is.null(cl)) {
        
        parallel::stopCluster(cl)
      }
    },
    add = TRUE
  )
  
  
  #==========================================================
  # Parallelization over pairs
  #==========================================================
  
  if (is.null(cl)) {
    
    pair_fit <- lapply(
      seq_len(nrow(pair_id)),
      fit_one_pair
    )
    
  } else {
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "YS",
        "R_t_initial",
        "R_s_initial",
        "pair_id",
        "target_cross",
        "A2_list",
        "K",
        "c_A1",
        "c_M1",
        "c_M2",
        "c_H",
        "density_floor",
        "fit_one_pair"
      ),
      envir = environment()
    )
    
    pair_fit <- parallel::parLapply(
      cl,
      seq_len(nrow(pair_id)),
      fit_one_pair
    )
  }
  
  
  #==========================================================
  # Collect
  #==========================================================
  
  rho_pair <- vapply(
    pair_fit,
    `[[`,
    numeric(1),
    "rho"
  )
  
  H_pair <- vapply(
    pair_fit,
    `[[`,
    numeric(1),
    "H"
  )
  
  
  #==========================================================
  # Correlation matrix
  #==========================================================
  
  R_hat <- diag(p)
  
  for (ii in seq_len(nrow(pair_id))) {
    
    j <- pair_id[ii, 1]
    l <- pair_id[ii, 2]
    
    R_hat[j, l] <- rho_pair[ii]
    R_hat[l, j] <- rho_pair[ii]
  }
  
  
  #==========================================================
  # Minimum denominator
  #==========================================================
  
  H_min <- min(H_pair)
  
  
  list(
    R_hat = R_hat,
    H_min = H_min,
    rho_pair = rho_pair,
    H_pair = H_pair,
    pair_id = pair_id
  )
}

# system.time({f1 = fit.ETL.corr(X0,XS,ncores = 18)})

#============================================================
# Sparse ETL-corr
#============================================================

fit.ETL.Scor <- function(
    R_hat,
    H_min,
    c_T = 1) {
  
  p <- ncol(R_hat)
  
  threshold <-
    c_T *
    sqrt(
      log(p) / H_min
    )
  
  R_scor <- R_hat
  
  id <- upper.tri(
    R_scor
  )
  
  R_scor[
    id &
      abs(R_scor) < threshold
  ] <- 0
  
  
  # Enforce symmetry
  R_scor[
    lower.tri(R_scor)
  ] <-
    t(R_scor)[
      lower.tri(R_scor)
    ]
  
  diag(R_scor) <- 1
  return(R_scor)
}

fit.target.scor = function(X0, corr = TRUE){
  
  sigma.hat = cov(X0)
  
  c_T = 1
  threshold <-
    c_T *
    sqrt(
      log(p) / nrow(X0)
    )
  sigma.thre.hat = sigma.hat
  sigma.thre.hat[which(abs(sigma.hat) < threshold)] = 0
  
  if(corr == TRUE){
    R.thre.hat = cov2cor(sigma.thre.hat)
  }else{
    
    R.thre.hat = sigma.thre.hat
  }
  
  return(R.thre.hat)
}


