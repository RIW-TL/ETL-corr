library(ks)
library(parallel)

#============================================================
# 1. ETL-corr for a grid of M. Parallelization only over (j,l)
#============================================================

fit.ETL.corr <- function(
    X0,
    XS,
    c_A1 = 0.4,
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
  
  R_t_initial <- fit.domain.scor(X0, corr = TRUE)
  R_s_initial <- list()
  for (k in 1:K) {
    R_s_initial[[k]] = fit.domain.scor(XS[[k]], corr = TRUE)
  }
  # R_s_initial <- R_s_true
  # R_t_initial <- cor(X0)
  # R_s_initial <- list()
  # for (k in 1:K) {
  #   R_s_initial[[k]] = cor(XS[[k]])
  # }
  
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


fit.ETL.corr.cv <- function(
    X0,
    XS,
    c_A1 = 0.8,
    c_A2 = 0.95,
    M_grid = c(0.3, 0.6, 0.9),
    c_H = 1,
    density_floor = 1e-10,
    ncores = 8,
    cl = NULL,
    cv_folds = 5) {
  
  M_grid <- sort(unique(M_grid))
  if (anyNA(M_grid) || !length(M_grid) || any(M_grid < 0)) {
    stop("M_grid must contain nonnegative finite values.")
  }
  if (any(!is.finite(M_grid))) {
    stop("M_grid must contain finite values.")
  }
  if (length(cv_folds) != 1L || is.na(cv_folds) ||
      cv_folds < 2L || cv_folds != as.integer(cv_folds) ||
      cv_folds >= nrow(X0)) {
    stop("cv_folds must be an integer between 2 and nrow(X0)-1.")
  }
  cv_folds <- as.integer(cv_folds)
  
  p <- ncol(X0)
  K <- length(XS)
  nM <- length(M_grid)
  max_M <- max(M_grid)
  
  #==========================================================
  # Initial estimators (full data)
  #==========================================================
  R_t_initial <- fit.domain.scor(X0, corr = TRUE)
  R_s_initial <- lapply(XS, function(Xk) {
    fit.domain.scor(Xk, corr = TRUE)
  })
  
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
  # A2: precomputed once per source and coordinate
  #==========================================================
  A2_list <- lapply(YS, function(Yk) {
    apply(abs(Yk), 2, quantile, probs = c_A2, names = FALSE)
  })
  
  #==========================================================
  # Target-only cross-validation splits
  #==========================================================
  fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(X0)))
  cv_target_cross <- vector("list", cv_folds)
  cv_valid_cross <- vector("list", cv_folds)
  cv_R_t_initial <- vector("list", cv_folds)
  
  for (f in seq_len(cv_folds)) {
    Xtr <- X0[fold_id != f, , drop = FALSE]
    Xva <- X0[fold_id == f, , drop = FALSE]
    
    Ytr <- scale(Xtr)
    Yva <- scale(
      Xva,
      center = attr(Ytr, "scaled:center"),
      scale = attr(Ytr, "scaled:scale")
    )
    
    cv_target_cross[[f]] <- crossprod(Ytr)
    cv_valid_cross[[f]] <- crossprod(Yva)
    cv_R_t_initial[[f]] <- fit.domain.scor(Xtr, corr = TRUE)
  }
  
  #==========================================================
  # One pair: all M values and all CV folds together
  #==========================================================
  fit_one_pair <- function(ii) {
    j <- pair_id[ii, 1]
    l <- pair_id[ii, 2]
    
    # Configuration 1 is the full-data fit; 2,... are CV folds.
    rho_target <- c(
      R_t_initial[j, l],
      vapply(cv_R_t_initial, function(R) R[j, l], numeric(1))
    )
    target_num <- c(
      target_cross[l, j],
      vapply(cv_target_cross, function(G) G[l, j], numeric(1))
    )
    target_den <- c(
      target_cross[l, l],
      vapply(cv_target_cross, function(G) G[l, l], numeric(1))
    )
    nc <- length(rho_target)
    
    num <- matrix(rep(target_num, each = nM), nrow = nM)
    den <- matrix(rep(target_den, each = nM), nrow = nM)
    
    for (k in seq_len(K)) {
      yj <- YS[[k]][, j]
      yl <- YS[[k]][, l]
      nk <- length(yl)
      rho_s <- R_s_initial[[k]][j, l]
      
      error_1 <- yj - rho_s * yl
      A1 <- as.numeric(quantile(
        abs(error_1), probs = c_A1, names = FALSE
      ))
      A2 <- A2_list[[k]][l]
      
      # First two selection constraints do not involve M.
      base_id <- which(abs(error_1) <= A1 & abs(yl) <= A2)
      if (!length(base_id)) next
      
      # Retain the largest-M set for each target fit.
      id_by_config <- lapply(rho_target, function(rho_t) {
        base_id[
          abs((rho_s - rho_t) * yl[base_id]) <= max_M
        ]
      })
      union_id <- sort(unique(unlist(id_by_config, use.names = FALSE)))
      if (!length(union_id)) next
      
      # Same source KDE and bandwidth as the original method.
      Z <- cbind(error_1, yl)
      H <- c_H * nk^(-1 / 3) * diag(c(
        max(1 - rho_s^2, 1e-8), 1
      ))
      H <- H + 1e-8 * diag(2)
      
      # Evaluate denominator once per union observation.
      # Evaluate numerator densities for all target fits.
      eval_list <- list(cbind(error_1[union_id], yl[union_id]))
      offset <- length(union_id)
      start_num <- integer(nc)
      
      for (a in seq_len(nc)) {
        ids <- id_by_config[[a]]
        if (!length(ids)) next
        
        m <- length(ids)
        error_0 <- yj[ids] - rho_target[a] * yl[ids]
        start_num[a] <- offset
        
        eval_list[[length(eval_list) + 1L]] <- cbind(
          error_0, yl[ids]
        )
        eval_list[[length(eval_list) + 1L]] <- cbind(
          -error_0, yl[ids]
        )
        offset <- offset + 2L * m
      }
      
      f_all <- as.numeric(ks::kde(
        x = Z,
        H = H,
        eval.points = do.call(rbind, eval_list),
        binned = TRUE,
        bgridsize = c(41, 41),
        compute.cont = FALSE
      )$estimate)
      f_den <- pmax(f_all[seq_along(union_id)], density_floor)
      
      for (a in seq_len(nc)) {
        ids <- id_by_config[[a]]
        m <- length(ids)
        if (!m) next
        
        s <- start_num[a]
        f_num1 <- f_all[s + seq_len(m)]
        f_num2 <- f_all[s + m + seq_len(m)]
        weight <- (f_num1 + f_num2) /
          (2 * f_den[match(ids, union_id)])
        
        wxy <- weight * yj[ids] * yl[ids]
        wyy <- weight * yl[ids]^2
        delta_abs <- abs((rho_s - rho_target[a]) * yl[ids])
        
        for (b in seq_len(nM)) {
          keep <- delta_abs <= M_grid[b]
          num[b, a] <- num[b, a] + sum(wxy[keep])
          den[b, a] <- den[b, a] + sum(wyy[keep])
        }
      }
    }
    
    # Validation prediction loss: Y_j predicted by rho * Y_l.
    cv_loss <- numeric(nM)
    for (f in seq_len(cv_folds)) {
      rho_f <- num[, f + 1L] / den[, f + 1L]
      G <- cv_valid_cross[[f]]
      cv_loss <- cv_loss +
        G[j, j] - 2 * rho_f * G[j, l] + rho_f^2 * G[l, l]
    }
    
    list(
      rho = num[, 1L] / den[, 1L],
      H = den[, 1L],
      cv_loss = cv_loss
    )
  }
  
  #==========================================================
  # Parallel cluster
  #==========================================================
  own_cluster <- FALSE
  
  if (is.null(cl) && ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    
    user_lib <- path.expand("~/R/library")
    parallel::clusterCall(cl, function(lib) {
      .libPaths(c(lib, .libPaths()))
      invisible(NULL)
    }, user_lib)
    
    own_cluster <- TRUE
  }
  
  on.exit({
    if (own_cluster && !is.null(cl)) {
      parallel::stopCluster(cl)
    }
  }, add = TRUE)
  
  #==========================================================
  # Parallelization over pairs
  #==========================================================
  if (is.null(cl)) {
    pair_fit <- lapply(seq_len(nrow(pair_id)), fit_one_pair)
  } else {
    parallel::clusterExport(
      cl,
      varlist = c(
        "YS", "R_t_initial", "R_s_initial", "pair_id",
        "target_cross", "A2_list", "K", "c_A1",
        "M_grid", "max_M", "nM", "c_H", "density_floor",
        "cv_target_cross", "cv_valid_cross", "cv_R_t_initial",
        "cv_folds", "fit_one_pair"
      ),
      envir = environment()
    )
    
    pair_fit <- parallel::parLapplyLB(
      cl, seq_len(nrow(pair_id)), fit_one_pair
    )
  }
  
  #==========================================================
  # Choose M by target CV prediction loss
  #==========================================================
  cv_loss <- Reduce(`+`, lapply(pair_fit, `[[`, "cv_loss")) /
    (nrow(X0) * nrow(pair_id))
  names(cv_loss) <- as.character(M_grid)
  best <- which.min(cv_loss)
  M_selected <- M_grid[best]
  
  #==========================================================
  # Collect the full-data fit for the selected M
  #==========================================================
  rho_pair <- vapply(pair_fit, function(z) z$rho[best], numeric(1))
  H_pair <- vapply(pair_fit, function(z) z$H[best], numeric(1))
  
  R_hat <- diag(p)
  for (ii in seq_len(nrow(pair_id))) {
    j <- pair_id[ii, 1]
    l <- pair_id[ii, 2]
    R_hat[j, l] <- rho_pair[ii]
    R_hat[l, j] <- rho_pair[ii]
  }
  
  H_min <- min(H_pair)
  
  list(
    R_hat = R_hat,
    H_min = H_min,
    rho_pair = rho_pair,
    H_pair = H_pair,
    pair_id = pair_id,
    M_selected = M_selected,
    M_grid = M_grid,
    cv_loss = cv_loss
  )
}

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

fit.domain.scor = function(X, corr = TRUE){
  
  sigma.hat = cov(X)
  p = ncol(X)
  
  c_T = 1
  threshold <-
    c_T *
    sqrt(
      log(p) / nrow(X)
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


# #============================================================
# # Implementation
# #============================================================
# 
# #==============================================
# # data generation
# #==============================================
# ii = ii + 1
# dat <- generate_sim3_data(eta = 0.7, tau = 0, seed = ii)
# 
# X0 <- dat$X0
# XS <- dat$XS
# p = ncol(dat$X0)
# R0 = dat$R0
# 
# fit <- fit.ETL.corr(X0, XS, c_A1 = 0.6, M = 0.4)
# ETL_corr_hat = fit$R_hat
# H_min = fit$H_min
# 
# fit.cv <- fit.ETL.corr(X0, XS, c_A1 = 0.6)
# ETL_corr_hat_cv = fit.cv$R_hat
# 
# # ETL-Scor
# ETL_Scor_hat = fit.ETL.Scor(R_hat = ETL_corr_hat, H_min, c_T = 1)
# 
# # Target-corr
# Target_corr_hat = cor(X0)
# 
# # Target-Scor
# Target_scor_hat = fit.domain.scor(X0, corr = TRUE)
# 
# 
# # Pool-corr
# X_all <- do.call(rbind,
#                  c(list(X0), XS))
# 
# Pool_corr_hat <- cor(X_all)
# 
# result = c(
#   norm(ETL_corr_hat - R0, "F") ^ 2 / p,
#   norm(ETL_corr_hat_cv - R0,"F")^2/p,
#   norm(ETL_Scor_hat - R0, "F") ^ 2 / p,
#   norm(Target_corr_hat - R0, "F") ^ 2 / p,
#   norm(Pool_corr_hat - R0, "F") ^ 2 / p,
#   norm(Target_scor_hat - R0, "F") ^ 2 / p
# )
# names(result) = c("ETL-corr","ETL-corr-2",  "ETL-Scor","Target-corr", "Pool_corr_hat","Target_scor_hat")
# result


