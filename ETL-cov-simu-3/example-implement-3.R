
cov.loss <- function(Sigma.hat, Sigma.true) {
  
  # dimension
  p <- nrow(Sigma.true)
  
  # basic checks
  if (!all(dim(Sigma.hat) == dim(Sigma.true))) {
    stop("Sigma.hat and Sigma.true must have the same dimensions.")
  }
  
  # -------------------------------------------------------
  # 1. Covariance matrix: scaled squared Frobenius loss
  # -------------------------------------------------------
  loss.cov <- norm(Sigma.hat - Sigma.true, type = "F")^2 / p
  
  # -------------------------------------------------------
  # 2. Correlation matrix: scaled squared Frobenius loss
  # -------------------------------------------------------
  R.hat  <- cov2cor(Sigma.hat)
  R.true <- cov2cor(Sigma.true)
  
  loss.cor <- norm(R.hat - R.true, type = "F")^2 / p
  
  # -------------------------------------------------------
  # 3. Diagonal variance matrix: scaled squared Frobenius loss
  # -------------------------------------------------------
  Lambda.hat  <- diag(diag(Sigma.hat))
  Lambda.true <- diag(diag(Sigma.true))
  
  loss.diag <- norm(Lambda.hat - Lambda.true, type = "F")^2
  
  # output
  return(c(
    Covariance = loss.cov,
    Correlation = loss.cor,
    Diagonal = loss.diag
  ))
}

new.fit = function(n0 = 100,n1 = 500, p, K = 5, eta = 0.4, d = 0.5, tau = 0,
                   c_A1 = 0.4, m2_A = 1, iter = 50){
  
  out.loss.ave = out.loss.corr.ave = out.loss.var.ave = matrix(0,iter,5)
  out.time.ave = matrix(0,iter,3)
  
  colnames(out.loss.ave) = colnames(out.loss.corr.ave) = colnames(out.loss.var.ave) =
    c("ETL-cov","ETL-scov", "Target-cov","Pool-cov","Target-scov")
  
  colnames(out.time.ave) = c("ETL-corr","ETL-scor", "ETL-var")
  
  for (ii in 1:iter) {
    
    ##==============================
    ## data generation
    ##==============================
    params = generate_etl_params(K = K, eta = eta, d = d)
    fit = generate_etl_data(params, nk = n1, distribution = "gaussian", tau = tau, seed = ii)
    X0 = fit$X0
    XS = fit$XS
    
    XS <- simplify2array(XS)
    
    Sigma0 = params$Sigma0
    diag.elements = diag(Sigma0)
    R_t = diag(diag.elements^(-1/2)) %*% Sigma0 %*% diag(diag.elements^(-1/2))
    
    
    ##=============================
    ## Tar-cov
    ##=============================
    sigma0.target.hat = cov(X0)
    
    ##=============================
    ## Target-Scov
    ##=============================
    target_Scov_hat = est_sparseCov(X0, method = 'cv', operator = 'hard', corr = FALSE)
    diag(target_Scov_hat) <- apply(X0, 2, var)
    
    ##=============================
    ## ETL-cov
    ##=============================
    
    # Step 1: ETL-corr estimator
    #-----------------------------------------
    
    ETL.start.time = Sys.time()
    
    R_t_initial = cor(X0)
    R_s_initial = array(0,c(p,p,K))
    for (k in 1:K) {R_s_initial[,,k] = cor(XS[,,k])}
    
    precomputed = precompute.ETL.corr.parallel(
      X0,
      XS,
      R_t_initial,
      R_s_initial,
      c_h2 = 1,
      density_floor = 1e-12,
      ncores = parallel::detectCores() - 2)
    
    fit.ETL.corr = fit.ETL.corr.fast(precomputed,
                                     c_A1, # c_A1 \in \{0.4,0.6,0.8\}
                                     c_A2 = 1,   # c_A2 \in \{0.6,0.8,1\}
                                     c_M1 = 0.5,
                                     c_M2 = 0.1)
    trans_corr_hat = fit.ETL.corr[[1]]
    
    ETL.corr.end.time = Sys.time()

    #-----------------------------------------------
    
    S_num_min = fit.ETL.corr[[2]]
    
    
    ##=============================
    ## ETL-Scov
    ##=============================
    
    fit.scor = fit.ETL.Scor.cv(trans_corr_hat, X0, R_true = R_t, 
                               lambda.range = seq(0.02,0.3,0.04),
                               rule = "theory", c_0 = 0.5, S_num_min)
    trans_scor_hat = fit.scor$R_scor_hat
    
    # norm(trans_scor_hat - R_t,"F")^2/p
    
    ETL.scor.end.time = Sys.time()
    
    
    # Step 2: ETL-var estimator
    #-----------------------------------------
    ETL.var.start.time = Sys.time()
    
    trans.m1.hat = fit.mean.trans(X0, XS,
                                  c_m1_A = 0.6, # c_m1_A \in \{0.4,0.6,0.8\}
                                  M1 = 0.3)
    
    M2.range = c(0,0.2,0.4,0.6)
    err.var = rep(0,length(M2.range))
    result.var = matrix(0,p,length(M2.range))
    for (l in seq(M2.range)) {
      
      M2 = M2.range[l]
      
      trans.m2.hat = fit.m2.trans(X0, XS,
                                  m2_A, # c_m1_A \in \{0.4,0.6,0.8\}
                                  M2)
      result.var[,l] = trans.m2.hat - trans.m1.hat^2
      err.var[l] = sum((result.var[,l] - diag.elements)^2)
    }
    ETL.var.hat = result.var[, which.min(err.var)]
    
    # sum((ETL.var.hat - diag(Sigma0))^2)
    # sum((diag(cov(X0)) - diag(Sigma0))^2)
    
    ETL.var.end.time = Sys.time()
    #-----------------------------------------
    
    
    ## Step 3: correlation-variance combination
    ##-----------------------------------------
    ETL.cov.hat = diag(ETL.var.hat^{1/2}) %*% trans_corr_hat %*% diag(ETL.var.hat^{1/2})
    ETL.scov.hat = diag(ETL.var.hat^{1/2}) %*% trans_scor_hat %*% diag(ETL.var.hat^{1/2})
    
    
    
    ##=============================
    ## Pool-cov
    ##=============================
    X.A = X0
    for (k in 1:K) {
      X.A = rbind(X.A,XS[,,k])
    }
    sigma0.pool.hat = cov(X.A)
    
    ##=============================
    ## loss
    ##=============================
    
    # loss of all methods on covariance matrix estimation
    out.loss = c(
      norm(ETL.cov.hat - Sigma0,"F")^2/p,
      norm(ETL.scov.hat - Sigma0,"F")^2/p,
      # norm(sigma0.ADT.hat - Sigma0,"F")^2/p,
      norm(sigma0.target.hat - Sigma0,"F")^2/p,
      norm(sigma0.pool.hat - Sigma0,"F")^2/p,
      norm(target_Scov_hat - Sigma0,"F")^2/p)
      # norm(target_band_hat - Sigma0,"F")^2/p)
    
    out.loss.corr = c(cov.loss(ETL.cov.hat, Sigma0)[2],
                      cov.loss(ETL.scov.hat, Sigma0)[2],
                      # cov.loss(sigma0.ADT.hat, Sigma0)[2],
                      cov.loss(sigma0.target.hat, Sigma0)[2],
                      cov.loss(sigma0.pool.hat, Sigma0)[2],
                      cov.loss(target_Scov_hat, Sigma0)[2])
                      # cov.loss(target_band_hat, Sigma0)[2])
    
    out.loss.var = c(cov.loss(ETL.cov.hat, Sigma0)[3],
                     cov.loss(ETL.scov.hat, Sigma0)[3],
                     # cov.loss(sigma0.ADT.hat, Sigma0)[3],
                     cov.loss(sigma0.target.hat, Sigma0)[3],
                     cov.loss(sigma0.pool.hat, Sigma0)[3],
                     cov.loss(target_Scov_hat, Sigma0)[3])
                     # cov.loss(target_band_hat, Sigma0)[3])
    
    out.time = c(
      as.numeric(ETL.corr.end.time - ETL.start.time),
      as.numeric(ETL.scor.end.time - ETL.start.time),
      as.numeric(ETL.var.end.time - ETL.var.start.time)
    )
    names(out.time) <- c("ETL-corr","ETL-scor", "ETL-var")
    
    out.loss.ave[ii,] = out.loss
    out.loss.corr.ave[ii,] = out.loss.corr
    out.loss.var.ave[ii,] = out.loss.var
    
    out.time.ave[ii,] = out.time
  }
  
  result.list = list()
  result.list[[1]] = out.loss.ave
  result.list[[2]] = out.loss.corr.ave
  result.list[[3]] = out.loss.var.ave
  
  return(list(result.all = result.list, 
              out.cov.F.ave = colMeans(out.loss.ave),
              out.corr.F.ave = colMeans(out.loss.corr.ave),
              out.loss.var.ave = colMeans(out.loss.var.ave),
              out.time1.ave = colMeans(out.time.ave)))
}


aa1 = new.fit(eta = 0,d = 0.3,m2_A = 1)
aa2 = new.fit(eta = 0.2,d = 0.3,m2_A = 1)
aa3 = new.fit(eta = 0.4,d = 0.3,m2_A = 1)
aa4 = new.fit(eta = 0.6,d = 0.3,m2_A = 1)
aa5 = new.fit(eta = 0.8,d = 0.3,m2_A = 1)

bb1 = new.fit(eta = 0,d = 0.5,m2_A = 1)
bb2 = new.fit(eta = 0.2,d = 0.5,m2_A = 1)
bb3 = new.fit(eta = 0.4,d = 0.5,m2_A = 1)
bb4 = new.fit(eta = 0.6,d = 0.5,m2_A = 1)
bb5 = new.fit(eta = 0.8,d = 0.5,m2_A = 1)

cc1 = new.fit(eta = 0,d = 0.7,m2_A = 1)
cc2 = new.fit(eta = 0.2,d = 0.7,m2_A = 1)
cc3 = new.fit(eta = 0.4,d = 0.7,m2_A = 1)
cc4 = new.fit(eta = 0.6,d = 0.7,m2_A = 1)
cc5 = new.fit(eta = 0.8,d = 0.7,m2_A = 1)

aaa1 = new.fit(eta = 0,d = 0.3,m2_A = 2)
aaa2 = new.fit(eta = 0.2,d = 0.3,m2_A = 2)
aaa3 = new.fit(eta = 0.4,d = 0.3,m2_A = 2)
aaa4 = new.fit(eta = 0.6,d = 0.3,m2_A = 2)
aaa5 = new.fit(eta = 0.8,d = 0.3,m2_A = 2)

bbb1 = new.fit(eta = 0,d = 0.5,m2_A = 2)
bbb2 = new.fit(eta = 0.2,d = 0.5,m2_A = 2)
bbb3 = new.fit(eta = 0.4,d = 0.5,m2_A = 2)
bbb4 = new.fit(eta = 0.6,d = 0.5,m2_A = 2)
bbb5 = new.fit(eta = 0.8,d = 0.5,m2_A = 2)

ccc1 = new.fit(eta = 0,d = 0.7,m2_A = 2)
ccc2 = new.fit(eta = 0.2,d = 0.7,m2_A = 2)
ccc3 = new.fit(eta = 0.4,d = 0.7,m2_A = 2)
ccc4 = new.fit(eta = 0.6,d = 0.7,m2_A = 2)
ccc5 = new.fit(eta = 0.8,d = 0.7,m2_A = 2)
