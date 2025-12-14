


fit.compare.PCA = function(num.info,d){
  
  n0 = 100;p = 100;n1 = 500;K = 5
  
  iter = 20
  out.cor = matrix(0,iter,4)
  out.cov = matrix(0,iter,5)
  
  colnames(out.cor) = c("single","pool","Trans-ECM","Trans-SECM")
  colnames(out.cov) = c("single","pool","Trans-ECM","Trans-SECM","GB")
  
  for (ii in 1:iter) {
    
    # data generation
    #-----------------------------------------
    data = dat.g.PCA(n0,n1,p,num.info,K,d)
    
    Xt = data$Xt;XS = data$XS
    rk_list = data$rk_list
    r0 = rk_list[[1]]
    
    Sig_t = data$Sig_t
    Sig_s = data$Sig_s
    
    R_t = diag(diag(Sig_t)^{-1/2}) %*% Sig_t %*% diag(diag(Sig_t)^{-1/2})
    R_s = array(0,c(p,p,K))
    for (k in 1:K) {R_s[,,k] = diag(diag(Sig_s[,,k])^{-1/2}) %*% Sig_s[,,k] %*% diag(diag(Sig_s[,,k])^{-1/2})}
    
    fit  = eigen(R_t)
    U_c_t = fit$vectors[,1:r0]
    
    fit  = eigen(Sig_t)
    U_t = fit$vectors[,1:r0]
    
    Yt = scale(Xt)
    YS = array(0,c(n1,p,K))
    for (k in 1:K) {YS[,,k] = scale(XS[,,k])}
    
    #===========================cor==============================
    #============================================================
    
    # Target-PCA 
    #------------------------------------------
    #------------------------------------------
    sam.cor = cor(Xt)
    U_single = eigen(sam.cor)$vectors[,1:r0]
    out.single = fit.evaluate(U_single,U_c_t)
    
    
    # Trans-ECM-PCA 
    #------------------------------------------
    #------------------------------------------
    R_t_hat = sam.cor
    # R_t_hat = est_sparseCov(Xt, method = 'cv', operator = 'scad', corr = T)
    R_s_hat = R_s
    pre = precompute_c_independent(Yt, YS, K, R_t_hat, R_s_hat)
    
    M.range = c(seq(0,0.2,0.05),seq(0.3,0.7,0.1))
    A.range = seq(0.3,0.9,0.3)
    
    cv.loss = matrix(0, length(A.range), length(M.range))
    R_ECM = list()
    
    for (ir in seq_along(A.range)) {
      for (im in seq_along(M.range)) {
        
        A = A.range[ir]
        M = M.range[im]
        
        R_Trans_ECM_hat = fit.trans.ECM(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre)
        l = (ir - 1)*length(M.range) + im
        R_ECM[[l]] = R_Trans_ECM_hat
        
        #----------------
        U_trans_ecm = eigen(R_Trans_ECM_hat)$vectors[,1:r0]
        cv.loss[ir,im] = fit.evaluate(U_trans_ecm,U_c_t)
      }
    }
    
    min.id = which(cv.loss == min(cv.loss),arr.ind = TRUE)
    l.opt = (min.id[1,1] - 1)*length(M.range) + min.id[1,2]
    R_Trans_ECM_hat = R_ECM[[l.opt]]
    out.Trans.ECM = min(cv.loss)
    
    # Trans-SECM-PCA 
    #------------------------------------------
    #------------------------------------------
    
    lambda.range = seq(0.01,0.1,0.01) 
    cv.loss = rep(0,length(lambda.range))
    R_Trans_SECM_hat = R_Trans_ECM_hat
    R_SECM = list()
    
    for (m in 1:length(lambda.range)) {
      
      lambda = lambda.range[m]
      
      for (j in 1:(p - 1)) {
        for (l in (j + 1):p) {
          
          if(abs(R_Trans_SECM_hat[j,l]) < lambda){
            
            R_Trans_SECM_hat[j,l] = R_Trans_SECM_hat[l,j] = 0
            
          }
          
        }
      }
      
      R_SECM[[m]] = R_Trans_SECM_hat
      
      #--------------
      U_trans_secm = eigen(R_Trans_SECM_hat)$vectors[,1:r0]
      cv.loss[m] = fit.evaluate(U_trans_secm,U_c_t)
    }
    
    min.id = which.min(cv.loss)
    R_Trans_SECM_hat = R_SECM[[min.id]]
    out.Trans.SECM = min(cv.loss)
    
    # Pool-PCA 
    #------------------------------------------
    #------------------------------------------ 
    X.A = c()
    for (k in 1:K) {X.A = rbind(X.A,XS[,,k])}
    sum.X = rbind(Xt,X.A)
    
    pool.cor = cor(sum.X)
    U_pool = eigen(pool.cor)$vectors[,1:r0]
    out.pool = fit.evaluate(U_pool,U_c_t)
    
    out.cor[ii,] = c(out.single,out.pool,out.Trans.ECM,out.Trans.SECM)
    
    #===========================cov==============================
    #============================================================
    
    # Target-PCA 
    #------------------------------------------
    #------------------------------------------
    sam.cov = cov(Xt)
    U_single = eigen(sam.cov)$vectors[,1:r0]
    out.single = fit.evaluate(U_single,U_t)
    
    
    # Trans-ECM-PCA 
    #------------------------------------------
    #------------------------------------------
    diag.elements = diag(Sig_t)
    cov_Trans_ECM_hat = diag(diag.elements^{1/2}) %*% R_Trans_ECM_hat %*% diag(diag.elements^{1/2})
    U_trans_ecm = eigen(cov_Trans_ECM_hat)$vectors[,1:r0]
    out.Trans.ECM = fit.evaluate(U_trans_ecm,U_t)
    
    
    # Trans-SECM-PCA
    #------------------------------------------
    #------------------------------------------
    cov_Trans_SECM_hat = diag(diag.elements^{1/2}) %*% R_Trans_SECM_hat %*% diag(diag.elements^{1/2})
    U_trans_secm = eigen(cov_Trans_SECM_hat)$vectors[,1:r0]
    out.Trans.SECM = fit.evaluate(U_trans_secm,U_t)
    
    
    # Pool-PCA 
    #------------------------------------------
    #------------------------------------------ 
    pool.cov = cov(sum.X)
    U_pool = eigen(pool.cov)$vectors[,1:r0]
    out.pool = fit.evaluate(U_pool,U_t)
    
    # GB
    #-----------------
    #-----------------
    
    Sigma_list = list()
    Sigma_list[[1]] = cov(Xt)
    nk_list = c(n0)
    
    for(k in 1:K){
      
      Sigma_list[[k+1]] = cov(XS[,,k])
      nk_list[k + 1] = n1
    }
    
    #---------
    rs = 0.9*r0;
    tau_star = 0.6*rs
    GB_re = GB_Kmeans(Sigma_list, nk_list, rk_list,rs, p, tau = tau_star, 
                      Ps_initial = 1, n_select = 1, Ti = 50)
    U_GB = GB_re$Pksk_finetuned
    out.GB = fit.evaluate(U_GB,U_t)
    
    
    #--------------------------
    out.cov[ii,] = c(out.single,out.pool,out.Trans.ECM,
                     out.Trans.SECM,out.GB)
    
  }
  
  return(list(out.cor = out.cor,
              out.cov = out.cov))
}


d = 0.1
num.info = 0
a0 = fit.compare.PCA(num.info,d)

num.info = 1
a1 = fit.compare.PCA(num.info,d)

# num.info = 2
# a2 = fit.compare.PCA(num.info,d)

num.info = 3
a3 = fit.compare.PCA(num.info,d)

# num.info = 4
# a4 = fit.compare.PCA(num.info,d)

num.info = 5
a5 = fit.compare.PCA(num.info,d)
#--------------------------------
d = 0.2
num.info = 0
b0 = fit.compare.PCA(num.info,d)

num.info = 1
b1 = fit.compare.PCA(num.info,d)

# num.info = 2
# b2 = fit.compare.PCA(num.info,d)

num.info = 3
b3 = fit.compare.PCA(num.info,d)

# num.info = 4
# b4 = fit.compare.PCA(num.info,d)

num.info = 5
b5 = fit.compare.PCA(num.info,d)
#--------------------------------
d = 0.3
num.info = 0
c0 = fit.compare.PCA(num.info,d)

num.info = 1
c1 = fit.compare.PCA(num.info,d)

# num.info = 2
# c2 = fit.compare.PCA(num.info,d)

num.info = 3
c3 = fit.compare.PCA(num.info,d)

# num.info = 4
# c4 = fit.compare.PCA(num.info,d)

num.info = 5
c5 = fit.compare.PCA(num.info,d)
#--------------------------------


