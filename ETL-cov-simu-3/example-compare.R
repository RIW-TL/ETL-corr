
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

band.fit = function(n0 = 100,n1 = 500,p = 100, eta, C_t = 0.2,
                    c_m1_A = 0.7, c_m2_A = 0.7, M1 = 0.3, M2 = 0.6, iter = 100){
  
  out.loss.ave = out.loss.corr.ave = out.loss.var.ave = matrix(0,iter,7)
  out.time.ave = matrix(0,iter,3)
  
  colnames(out.loss.ave) = colnames(out.loss.corr.ave) = colnames(out.loss.var.ave) =
    c("ETL-cov","ETL-scov", "ADT", "Target-cov","Pool-cov","Target-scov","Target-band")
  
  colnames(out.time.ave) = c("ETL-corr","ETL-scor", "ETL-var")
  
  for (ii in 1:iter) {
    
    ##==============================
    ## data generation
    ##==============================
    fit = dat.g.band(n0, n1, p, eta, C_t, seed = 100 + ii) 
    X0 = fit$target_data
    XS = fit$source_data
    
    Sigma0 = fit$Sigma0
    diag.elements = diag(Sigma0)
    R_t = diag(diag.elements^(-1/2)) %*% Sigma0 %*% diag(diag.elements^(-1/2))
    
    
    ##=============================
    ## ADT
    ##=============================
    fit.ADT = ADT(target_data = X0, source_data = XS)
    sigma0.ADT.hat = fit.ADT$estimate
    
    # norm(sigma0.ADT.hat - Sigma0,"F")^2/p
    
    ##============================
    ## Target-band
    ##============================
    target_band_hat = fit_band(X0)$estimate
    
    
    ##=============================
    ## ETL-cov
    ##=============================
    
    # Step 1: ETL-corr estimator
    #-----------------------------------------
    
    ETL.start.time = Sys.time()
    
    R_t_initial = fit_band(X0)$estimate
    R_s_initial = fit_band(XS)$estimate
    
    precomputed = precompute.ETL.corr.parallel(
      X0,
      XS,
      R_t_initial,
      R_s_initial,
      c_h2 = 1,
      density_floor = 1e-12,
      ncores = parallel::detectCores() - 2)
    
    fit.ETL.corr = fit.ETL.corr.fast(precomputed,
                                     c_A1 = 0.6, # c_A1 \in \{0.4,0.6,0.8\}
                                     c_A2 = 1,   # c_A2 \in \{0.6,0.8,1\}
                                     c_M1 = 0.5,
                                     c_M2 = 1)
    trans_corr_hat = fit.ETL.corr[[1]]
    
    ETL.corr.end.time = Sys.time()
    #-----------------------------------------------
    
    S_num_min = fit.ETL.corr[[2]]
    
    
    ##=============================
    ## ETL-Scov
    ##=============================
    
    fit.scor = fit.ETL.Scor.cv(trans_corr_hat, X0, R_true = R_t, 
                               lambda.range = seq(0.02,0.3,0.04),
                               rule = "theory", c_0 = 1, S_num_min)
    trans_scor_hat = fit.scor$R_scor_hat
    
    ETL.scor.end.time = Sys.time()
    
    
    # Step 2: ETL-var estimator
    #-----------------------------------------
    ETL.var.start.time = Sys.time()
    
    # trans.m1.hat = fit.mean.trans(X0, XS, 
    #                               c_m1_A = 0.8, # c_m1_A \in \{0.4,0.6,0.8\}
    #                               c_m1_h = 1, 
    #                               M = 0.3)  
    # 
    # trans.m2.hat = fit.m2.trans(X0, XS, 
    #                             c_m2_A = 0.6, # c_m1_A \in \{0.4,0.6,0.8\}
    #                             c_m2_h = 0.5, 
    #                             M = 0.6)
    trans.m1.hat = fit.mean.trans(X0, XS,
                                  c_m1_A, # c_m1_A \in \{0.4,0.6,0.8\}
                                  M1)
    
    trans.m2.hat = fit.m2.trans(X0, XS,
                                c_m2_A, # c_m1_A \in \{0.4,0.6,0.8\}
                                M2)
    ETL.var.hat = trans.m2.hat - trans.m1.hat^2
    
    ETL.var.end.time = Sys.time()
    #-----------------------------------------
    
    
    ## Step 3: correlation-variance combination
    ##-----------------------------------------
    ETL.cov.hat = diag(ETL.var.hat^{1/2}) %*% trans_corr_hat %*% diag(ETL.var.hat^{1/2})
    ETL.scov.hat = diag(ETL.var.hat^{1/2}) %*% trans_scor_hat %*% diag(ETL.var.hat^{1/2})
    
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
    ## Pool-cov
    ##=============================
    sigma0.pool.hat = cov(rbind(X0,XS))
    
    ##=============================
    ## loss
    ##=============================
    
    # loss of all methods on covariance matrix estimation
    out.loss = c(
      norm(ETL.cov.hat - Sigma0,"F")^2/p,
      norm(ETL.scov.hat - Sigma0,"F")^2/p,
      norm(sigma0.ADT.hat - Sigma0,"F")^2/p,
      norm(sigma0.target.hat - Sigma0,"F")^2/p,
      norm(sigma0.pool.hat - Sigma0,"F")^2/p,
      norm(target_Scov_hat - Sigma0,"F")^2/p,
      norm(target_band_hat - Sigma0,"F")^2/p)
    
    out.loss.corr = c(cov.loss(ETL.cov.hat, Sigma0)[2],
                      cov.loss(ETL.scov.hat, Sigma0)[2],
                      cov.loss(sigma0.ADT.hat, Sigma0)[2],
                      cov.loss(sigma0.target.hat, Sigma0)[2],
                      cov.loss(sigma0.pool.hat, Sigma0)[2],
                      cov.loss(target_Scov_hat, Sigma0)[2],
                      cov.loss(target_band_hat, Sigma0)[2])
    
    out.loss.var = c(cov.loss(ETL.cov.hat, Sigma0)[3],
                     cov.loss(ETL.scov.hat, Sigma0)[3],
                     cov.loss(sigma0.ADT.hat, Sigma0)[3],
                     cov.loss(sigma0.target.hat, Sigma0)[3],
                     cov.loss(sigma0.pool.hat, Sigma0)[3],
                     cov.loss(target_Scov_hat, Sigma0)[3],
                     cov.loss(target_band_hat, Sigma0)[3])
    
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


b1 = band.fit(eta = 0.1,c_m2_A = 0.6)
b2 = band.fit(eta = 0.3,c_m2_A = 0.6)
b3 = band.fit(eta = 0.5,c_m2_A = 0.6)
b4 = band.fit(eta = 0.7,c_m2_A = 0.6)
b5 = band.fit(eta = 0.9,c_m2_A = 0.6)

b11 = band.fit(eta = 0.1,c_m2_A = 0.7)
b12 = band.fit(eta = 0.3,c_m2_A = 0.7)
b13 = band.fit(eta = 0.5,c_m2_A = 0.7)
b14 = band.fit(eta = 0.7,c_m2_A = 0.7)
b15 = band.fit(eta = 0.9,c_m2_A = 0.7)

b21 = band.fit(eta = 0.1,c_m2_A = 0.8)
b22 = band.fit(eta = 0.3,c_m2_A = 0.8)
b23 = band.fit(eta = 0.5,c_m2_A = 0.8)
b24 = band.fit(eta = 0.7,c_m2_A = 0.8)
b25 = band.fit(eta = 0.9,c_m2_A = 0.8)

c1 = band.fit(eta = 0.1,c_m2_A = 0.3)
c2 = band.fit(eta = 0.3,c_m2_A = 0.3)
c3 = band.fit(eta = 0.5,c_m2_A = 0.3)
c4 = band.fit(eta = 0.7,c_m2_A = 0.3)
c5 = band.fit(eta = 0.9,c_m2_A = 0.3)

c11 = band.fit(eta = 0.1,c_m2_A = 0.4)
c12 = band.fit(eta = 0.3,c_m2_A = 0.4)
c13 = band.fit(eta = 0.5,c_m2_A = 0.4)
c14 = band.fit(eta = 0.7,c_m2_A = 0.4)
c15 = band.fit(eta = 0.9,c_m2_A = 0.4)

c21 = band.fit(eta = 0.1,c_m2_A = 0.5)
c22 = band.fit(eta = 0.3,c_m2_A = 0.5)
c23 = band.fit(eta = 0.5,c_m2_A = 0.5)
c24 = band.fit(eta = 0.7,c_m2_A = 0.5)
c25 = band.fit(eta = 0.9,c_m2_A = 0.5)




