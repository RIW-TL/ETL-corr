
# parameter set
#---------------------------------------------

fit.GGM.new = function(num.info,d,iter){
  
  n0 = 100;n1 = 500;p = 50;K = 5
  
  out.cor = matrix(0,iter,4)
  out.cov = matrix(0,iter,5)
  colnames(out.cor) = c("Trans-ECM","Trans-SECM","Target","Pool")
  colnames(out.cov) = c("Trans-ECM","Trans-SECM","Target","Pool","Trans-CLIME")
  
  for (ii in 1:iter) {
  
    # data generation
    #---------------------------------------------
    data = dat.ggm.g(n0,n1,p,num.info,K,d)
    Xt = data$Xt
    XS = data$XS
    
    Yt <- scale(Xt)
    YS <- array(0, c(n1, p, K))
    for (k in 1:K) {
      YS[,,k] <- scale(XS[,,k])
    }
    
    Ome_t = data$Ome_t
    Sig_t = solve(Ome_t)
    Ome_s = data$Ome_s
    
    diag.elements = diag(Sig_t)
    R_t = diag(diag.elements^{-1/2}) %*% Sig_t %*% diag(diag.elements^{-1/2})
    Ome_c_t = solve(R_t)
    
    R_s = array(0,c(p,p,K))
    for (k in 1:K) {
      
      Sig_s_k = solve(Ome_s[,,k])
      diag.elements = diag(Sig_s_k)
      R_s[,,k] = diag(diag.elements^{-1/2}) %*% Sig_s_k %*% diag(diag.elements^{-1/2})
    }
    
    
    
    #================================corr=====================================
    #=========================================================================
    
    # Pooled-Glasso
    #----------------------------------------------
    #----------------------------------------------
    X.A = c()
    for (k in 1:K) {X.A = rbind(X.A,XS[,,k])}
    sum.X = rbind(Xt,X.A)
    
    pool.cor = cor(sum.X)
    # lam = seq(0.001,0.005,0.001)
    # fit = CVglasso(Yt,S = pool.cor, K = 5, trace = "none",
    #                lam = lam)
    # Ome_pool_c_hat = fit$Omega
    # norm(Ome_pool_c_hat - Ome_c_t,"F")^2/p
    
    
    lam.pool = seq(0.001,0.015,0.002)
    fit = glasso.cv(pool.cor,rho.range = lam.pool,Ome_c_t)
    Ome_pool_c_hat = fit$Omega 
    out.cor.pool = norm(Ome_pool_c_hat - Ome_c_t,"F")^2/p
    
    #   fit = CVglasso(Yt, pool.cor, K = 5, trace = "none", lam = lam.pool)
    #   Ome_pool_c_hat = fit$Omega
    # norm(Ome_pool_c_hat - Ome_c_t,"F")^2/p
    
    # Ome_single_c_hat = glasso(pool.cor, rho = 0.00005)$wi
    # norm(Ome_single_c_hat - Ome_c_t,"F")^2/p
    
    # Target-GlASSO
    #------------------------------------------
    #------------------------------------------
    sam.cor = cor(Xt)
    lam.single = seq(0.04,0.1,0.02)
    fit = CVglasso(Yt, S = sam.cor, K = 5, trace = "none", lam = lam.single) 
    Ome_single_c_hat = fit$Omega
    out.cor.single = norm(Ome_single_c_hat - Ome_c_t,"F")^2/p
    # 
    # Ome_single_c_hat = glasso(sam.cor, rho = 0.1)$wi
    # norm(Ome_single_c_hat - Ome_c_t,"F")^2/p
    
    
    # Trans-ECM-Glasso 
    #------------------------------------------
    #------------------------------------------
    R_t_hat = solve(Ome_single_c_hat)
    # R_t_hat = cor(Xt)
    # R_t_hat = R_t
    # R_t_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    # norm(R_t_hat - R_t,"F")^2/p
    R_s_hat = R_s
    
    pre = precompute_c_independent(Yt, YS, K, R_t_hat, R_s_hat)
    
    if(num.info == 0){
      
      rho.range = seq(0.1,0.3,0.1)
      M.range = seq(0,0.4,0.2)
      lam.ECM = seq(0.03,0.07,0.01)
      
    }else{
      
      rho.range = seq(0.1,0.3,0.1)
      M.range = seq(0.7,1,0.1)
      lam.ECM = seq(0.001,0.005,0.001)
    }
    
    
    O = matrix(0,length(rho.range),length(M.range))
    R_ECM = list()
    
    for (ir in 1:length(rho.range)) {
      for (im in 1:length(M.range)) {
        
        rho = rho.range[ir]
        rho.M = M.range[im] 
        
        R_Trans_ECM_hat = fit.trans.ECM.GGM(Yt, YS, K, rho.M, rho, R_t_hat, R_s_hat, pre)
        # O[ir,im] = norm(R_Trans_ECM_hat - R_t,"F")^2/p
        
        l = (ir - 1)*length(M.range) + im
        R_ECM[[l]] = R_Trans_ECM_hat
        #---------------
        
        # fit = clime.cv(R_Trans_ECM_hat,lam.ECM,Ome_c_t) 
        
        # Ome_ECM_c_hat = glasso(R_Trans_ECM_hat, rho = 0.07)$wi
        # norm(Ome_ECM_c_hat - Ome_c_t,"F")^2/p
        
        fit = CVglasso(Yt, S = R_Trans_ECM_hat, K = 5, trace = "none", lam = lam.ECM)
        Ome_ECM_c_hat = fit$Omega
        O[ir,im] = norm(Ome_ECM_c_hat - Ome_c_t,"F")^2/p
      }
    }
    
    min.id = which(O == min(O),arr.ind = TRUE)
    l.opt = (min.id[1,1] - 1)*length(M.range) + min.id[1,2]
    R_Trans_ECM_hat = R_ECM[[l.opt]]
    
    out.cor.ECM = min(O)
    
    # fit = CVglasso(Xt, S = R_Trans_ECM_hat, K = 5, trace = "none", lam = lam.ECM)
    # Ome_ECM_c_hat = fit$Omega
    # Ome_ECM_c_hat = glasso(R_Trans_ECM_hat,rho = 0.08)$wi
    # norm(Ome_ECM_c_hat - Ome_c_t,"F")^2/p
    
    # Trans-SECM-Glasso 
    #------------------------------------------
    
    lambda.range = seq(0.01,0.05,0.01) 
    O.SECM = rep(0,length(lambda.range))
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
      fit = CVglasso(Yt, S = R_Trans_SECM_hat, K = 5, trace = "none", lam = lam.ECM)
      # fit = clime.cv(R_Trans_ECM_hat,lam.ECM,Ome_c_t) 
      Ome_SECM_c_hat = fit$Omega
      O.SECM[m] = norm(Ome_SECM_c_hat - Ome_c_t,"F")^2/p
    }
    
    min.id = which.min(O.SECM)
    R_Trans_SECM_hat = R_SECM[[min.id]]
    out.cor.SECM = min(O.SECM)
    
    # fit = CVglasso(Yt, S = R_Trans_SECM_hat, K = 5, trace = "none", lam = lam.ECM)
    # Ome_SECM_c_hat = fit$Omega
    # Ome_SECM_hat = glasso(R_Trans_SECM_hat,rho = 0.0005)$wi
    # norm(Ome_SECM_hat - Ome_t,"F")^2/p
    
    #======================================cov===================================
    #============================================================================
    
    # Target-GlASSO
    #------------------------------
    #------------------------------
    sam.cov = cov(Xt)
    Ome_single_hat = CVglasso(Xt,S = sam.cov, K = 5 ,trace = "none",
                              lam = seq(0.05,0.12,0.02))$Omega
    # rho.range = seq(0.04,0.1,0.02)
    # fit = clime.cv(sam.cov,rho.range,Ome_t) 
    # Ome_single_hat = fit$Omega
    out.cov.single = norm(Ome_single_hat - Ome_t,"F")^2/p
    
    # Ome_single_hat = glasso(sam.cov, rho = 0.05)$wi
    # norm(Ome_single_hat - Ome_t,"F")^2/p
    
    
    # Trans-ECM-Glasso 
    #------------------------------------------
    #------------------------------------------
    diag.elements = diag(Sig_t)
    cov_Trans_ECM_hat = diag(diag.elements^{1/2}) %*% R_Trans_ECM_hat %*% diag(diag.elements^{1/2})
    Ome_ECM_hat = CVglasso(Xt,S = cov_Trans_ECM_hat, K = 5 ,trace = "none",
                           lam = lam.ECM)$Omega
    # # fit = clime.cv(cov_Trans_ECM_hat,lam.ECM,Ome_t) 
    # Ome_ECM_hat = fit$Omega
    out.cov.ECM = norm(Ome_ECM_hat - Ome_t,"F")^2/p
    
    # Ome_ECM_hat = glasso(cov_Trans_ECM_hat, rho = 0.08)$wi
    # norm(Ome_ECM_hat - Ome_t,"F")^2/p
    
    
    # Trans-SECM-Glasso 
    #------------------------------------------
    #------------------------------------------
    cov_Trans_SECM_hat = diag(diag.elements^{1/2}) %*% R_Trans_SECM_hat %*% diag(diag.elements^{1/2})
    Ome_SECM_hat = CVglasso(Xt,S = cov_Trans_SECM_hat, K = 5 ,trace = "none",
                            lam = lam.ECM)$Omega
    # fit = clime.cv(cov_Trans_SECM_hat,lam.ECM,Ome_t)
    # Ome_SECM_hat = fit$Omega
    out.cov.SECM = norm(Ome_SECM_hat - Ome_t,"F")^2/p
    
    # Pooled-Glasso
    #----------------------------------------------
    #----------------------------------------------
    X.A = c()
    for (k in 1:K) {X.A = rbind(X.A,XS[,,k])}
    sum.X = rbind(Xt,X.A)
    pool.cov = cov(sum.X)
    
    
    
    if(num.info == 0){
      
      lam.pool = seq(0.4,0.8,0.05)
    }else{
      
      lam.pool = seq(0.02,0.2,0.02)
    }
    fit = glasso.cv(pool.cov,rho.range = lam.pool,Ome_t)
    Ome_pool_hat = fit$Omega
    out.cov.pool = norm(Ome_pool_hat - Ome_t,"F")^2/p
    
    # Ome_pool_hat = glasso(pool.cov, rho = 0.8)$wi
    # norm(Ome_pool_hat - Ome_t,"F")^2/p
    
    # Trans-CLIME 
    #------------------------------------------
    #------------------------------------------
    n.vec = c(n0,rep(n1,K))
    n00 = round(n.vec[1]*4/5) 
    const = 0.5
    Theta.re0 = Myfastclime.s(X = Xt[1:n00,], Bmat = diag(1,p), lambda = const*2*sqrt(log(p)/n00))
    Theta.init = Theta.re0$Theta.hat
    Omega.tl1 = Trans.CLIME(X = Xt[1:n00,], X.A, const = const, 
                            X.til = Yt[(n00 + 1):n0,], Theta.cl = Theta.init)
    ind2 = (n0 - n00 + 1):n0
    Theta.re0 = Myfastclime.s(X = Xt[ind2,], Bmat = diag(1,p), lambda = const*2*sqrt(log(p)/length(ind2)))
    Theta.init = Theta.re0$Theta.hat
    Omega.tl2 = Trans.CLIME(X = Xt[ind2,], X.A, const = const,
                            X.til = Yt[1:(n0 - n00),], Theta.cl = Theta.init)
    Ome_trans_clime_hat = (Omega.tl1 + Omega.tl2)/2
    out.cov.trans.clime = norm(Ome_trans_clime_hat - Ome_t,"F")^2/p
    
    
    #-------------------------------------------------
    out.cor[ii,] = c(out.cor.ECM,out.cor.SECM,out.cor.single,out.cor.pool)
    out.cov[ii,] = c(out.cov.ECM,out.cov.SECM,out.cov.single,out.cov.pool,
                     out.cov.trans.clime)
  }
  
  return(list(out.cor = colMeans(out.cor),
              out.cov = colMeans(out.cov)))
}

d = 0.1;num.info = 0;iter = 50
a0 = fit.GGM.new(num.info,d,iter)

num.info = 1
a1 = fit.GGM.new(num.info,d,iter)

num.info = 3
a3 = fit.GGM.new(num.info,d,iter)

num.info = 5
a5 = fit.GGM.new(num.info,d,iter)
#------------------------------------
d = 0.2;num.info = 0
b0 = fit.GGM.new(num.info,d,iter)

num.info = 1
b1 = fit.GGM.new(num.info,d,iter)

num.info = 3
b3 = fit.GGM.new(num.info,d,iter)

num.info = 5
b5 = fit.GGM.new(num.info,d,iter)
#------------------------------------
d = 0.3;num.info = 0
c0 = fit.GGM.new(num.info,d,iter)

num.info = 1
c1 = fit.GGM.new(num.info,d,iter)

num.info = 3
c3 = fit.GGM.new(num.info,d,iter)

num.info = 5
c5 = fit.GGM.new(num.info,d,iter)
#------------------------------------
d = 0.4;num.info = 0
d0 = fit.GGM.new(num.info,d,iter)

num.info = 1
d1 = fit.GGM.new(num.info,d,iter)

num.info = 3
d3 = fit.GGM.new(num.info,d,iter)

num.info = 5
d5 = fit.GGM.new(num.info,d,iter)

