library(ks)
library(parallel)

#============================================================
# 1. ETL-corr for a grid of M. Parallelization only over (j,l)
#============================================================

fit.ETL.corr <- function(
    X0,
    XS,
    c_A1 = 0.8,
    c_A2 = 0.95,
    M = 0.3,
    c_H = 1,
    density_floor = 1e-10,
    ncores = 8,
    cl = NULL) {
  
  p <- ncol(X0)
  K <- length(XS)
  
  #==========================================================
  # Initial estimators
  #==========================================================

  R_t_initial <- cov2cor(fit_band(X0)$estimate)

  R_s_initial <- lapply(
    XS,
    function(X) cov2cor(fit_band(X)$estimate)
  )
  
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
        "M",
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

# system.time({f1 = fit.ETL.corr(X0,XS,ncores = 8)})

#============================================================
# Sparse ETL-corr
#============================================================

fit.ETL.Scor <- function(
    R_hat,
    H_min,
    c_T = 0.5) {
  
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

fit.ETL.Scor.adaptive <- function(
    R_hat,
    H_pair,
    pair_id,
    c_T = 1
) {
  
  p <- ncol(R_hat)
  
  R_scor <- R_hat
  
  for (ii in seq_len(nrow(pair_id))) {
    
    j <- pair_id[ii, 1]
    l <- pair_id[ii, 2]
    
    threshold <- c_T *
      sqrt(
        log(p) / H_pair[ii]
      )
    
    if (abs(R_hat[j, l]) < threshold) {
      
      R_scor[j, l] <- 0
      R_scor[l, j] <- 0
    }
  }
  
  diag(R_scor) <- 1
  
  R_scor
}

#============================================================
# Implementation
#============================================================

#==============================================
# data generation
#==============================================
ii = ii + 1
dat <- generate_sim2_data(eta = 0.5, seed = ii)

X0 <- dat$X0
XS <- dat$XS
K = length(XS)
R0 = cov2cor(dat$Sigma0)


# ETL-corr
fit <- fit.ETL.corr(X0, XS, c_A1 = 0.6, M = 0.8)

ETL_corr_hat = fit$R_hat
H_min = fit$H_min



# ETL-Scor
ETL_Scor_hat = fit.ETL.Scor(R_hat = ETL_corr_hat, H_min, c_T = 1)


H_pair = fit$H_pair
pair_id = fit$pair_id
ETL_Scor_hat_2 = fit.ETL.Scor.adaptive(R_hat = ETL_corr_hat, H_pair, pair_id)

# Target-corr
Target_band_hat = cov2cor(fit_band(X0)$estimate)

# Pool-corr
X_all <- do.call(rbind,
                 c(list(X0), XS))

Pool_corr_hat <- cor(X_all)

ADT.sigma.hat = ADT(X0,XS[[1]])$estimate
ADR.cor.hat = cov2cor(ADT.sigma.hat)

p = dim(X0)[2]
result = c(
  norm(ETL_corr_hat - R0, "F") ^ 2 / p,
  norm(ETL_Scor_hat - R0, "F") ^ 2 / p,
  norm(ETL_Scor_hat_2 - R0, "F") ^ 2 / p,
  norm(ADR.cor.hat - R0, "F") ^ 2 / p,
  norm(Target_band_hat - R0, "F") ^ 2 / p,
  norm(Pool_corr_hat - R0, "F") ^ 2 / p)
names(result) = c("ETL-corr", "ETL-Scor","ETL-Scor-2","ADT","Target_band", "Pool_corr_hat")
result


