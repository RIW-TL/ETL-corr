time.fit = function(n0 = 100,n1 = 500, p, K, eta = 0.4, d = 0.5, tau = 0,
                   c_A1 = 0.4, m2_A = 1, iter = 5){
  
  out.time.ave = matrix(0,iter,3)
  
  colnames(out.time.ave) = c("ETL-corr","ETL-scor", "ETL-var")
  
  for (ii in 1:iter) {
    
    ##==============================
    ## data generation
    ##==============================
    params = generate_etl_params(K = K, p = p, eta = eta, d = d)
    fit = generate_etl_data(params, nk = n1, distribution = "gaussian", tau = tau, seed = ii)
    X0 = fit$X0
    XS = fit$XS
    
    XS <- simplify2array(XS)
    
    Sigma0 = params$Sigma0
    diag.elements = diag(Sigma0)
    R_t = diag(diag.elements^(-1/2)) %*% Sigma0 %*% diag(diag.elements^(-1/2))
    
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
                                     c_A1, 
                                     c_A2 = 1,
                                     c_M1 = 0.5,
                                     c_M2 = 0.1)
    trans_corr_hat = fit.ETL.corr[[1]]
    
    ETL.corr.end.time = Sys.time()
    # ETL.corr.end.time - ETL.start.time
    #-----------------------------------------------
    
    S_num_min = fit.ETL.corr[[2]]
    
    
    ##=============================
    ## ETL-Scov
    ##=============================
    
    fit.scor = fit.ETL.Scor.cv(trans_corr_hat, X0, R_true = R_t, 
                               lambda.range = seq(0.02,0.3,0.04),
                               rule = "theory", c_0 = 0.5, S_num_min)
    trans_scor_hat = fit.scor$R_scor_hat
    
    ETL.scor.end.time = Sys.time()
    
    
    # Step 2: ETL-var estimator
    #-----------------------------------------
    ETL.var.start.time = Sys.time()
    
    trans.m1.hat = fit.mean.trans(X0, XS,
                                  c_m1_A = 0.6, # c_m1_A \in \{0.4,0.6,0.8\}
                                  M1 = 0.3)

      trans.m2.hat = fit.m2.trans(X0, XS,
                                  m2_A, # c_m1_A \in \{0.4,0.6,0.8\}
                                  M2 = 0.4)
      ETL.var.hat = trans.m2.hat - trans.m1.hat^2
    
    ETL.var.end.time = Sys.time()
    #-----------------------------------------
    
    
    ## Step 3: correlation-variance combination
    ##-----------------------------------------
    ETL.cov.hat = diag(ETL.var.hat^{1/2}) %*% trans_corr_hat %*% diag(ETL.var.hat^{1/2})
    ETL.scov.hat = diag(ETL.var.hat^{1/2}) %*% trans_scor_hat %*% diag(ETL.var.hat^{1/2})
    
    ##=============================
    ## loss
    ##=============================
    out.time = c(
      as.numeric(ETL.corr.end.time - ETL.start.time),
      as.numeric(ETL.scor.end.time - ETL.start.time),
      as.numeric(ETL.var.end.time - ETL.var.start.time)
    )
    names(out.time) <- c("ETL-corr","ETL-scor", "ETL-var")
    
    out.time.ave[ii,] = out.time
  }
  
  return(colMeans(out.time.ave))
}


a1 = time.fit(p = 50,K = 1)
a2 = time.fit(p = 100,K = 1)
a4 = time.fit(p = 200,K = 1)

b1 = time.fit(p = 50,K = 3)
b2 = time.fit(p = 100,K = 3)
b4 = time.fit(p = 200,K = 3)

c1 = time.fit(p = 50,K = 5)
c2 = time.fit(p = 100,K = 5)
c4 = time.fit(p = 200,K = 5)