
library(mvtnorm)

generate_sim4_data <- function(n0 = 100, n1 = 500, p = 50, 
                               k0 = 2, K = 5, d = 10, seed = NULL) {
  
  #==========================================================
  # Target covariance matrix
  #==========================================================
  
  U1 <- matrix(rnorm(p * p), p, p)
  qr_result_U <- qr(U1)
  U2 <- qr.Q(qr_result_U)
  
  r0 <- 5
  
  sigma0 <- diag(
    c(
      rep(5, r0),
      rep(0.5, p - r0)
    )
  )
  
  Sig_t <- as.matrix(
    U2 %*% sigma0 %*% t(U2)
  )
  
  
  #==========================================================
  # Source covariance matrices
  #==========================================================
  
  Sig_s <- vector("list", K)
  rk_list = list()
  rk_list[[1]] = r0
  
  for (k in 1:K) {
    
    if (k <= k0) {
      
      # Informative source
      U1k <- U2 + matrix(
        runif(p * p, 0, d / p),
        p, p
      )
      
      qr_result_Uk <- qr(U1k)
      U2k <- qr.Q(qr_result_Uk)
      
      sigmak <- sigma0
      rk_list[[k + 1]] = r0
      
    } else {
      
      # Non-informative source
      U1k <- U2 + matrix(
        runif(p * p, 1, 2),
        p, p
      )
      
      qr_result_Uk <- qr(U1k)
      U2k <- qr.Q(qr_result_Uk)
      
      rk <- 6
      rk_list[[k + 1]] = rk
      
      sigmak <- diag(
        c(
          rep(5, rk),
          rep(0.5, p - rk)
        )
      )
    }
    
    Sig_s[[k]] <- as.matrix(
      U2k %*% sigmak %*% t(U2k)
    )
  }
  
  
  #==========================================================
  # Data generation
  #==========================================================
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mu_0 <- rep(0, p)
  mu_1 <- rep(0, p)
  
  X0 <- mvtnorm::rmvnorm(
    n0,
    mu_0,
    Sig_t
  )
  
  XS <- vector("list", K)
  
  for (k in 1:K) {
    
    XS[[k]] <- mvtnorm::rmvnorm(
      n1,
      mu_1,
      Sig_s[[k]]
    )
  }
  
  U_t = U2[,1:r0]
  #==========================================================
  # Output
  #==========================================================
  
  list(
    X0 = X0,
    XS = XS,
    Sig_t = Sig_t,
    Sig_s = Sig_s,
    r0 = r0,
    rk_list = rk_list,
    U_t = U_t
  )
}
