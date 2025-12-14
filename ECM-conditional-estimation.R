
case2.fun.conditional = function(n0,n1,p,K,num.info,d,case,iter){
  
  cor.F = matrix(0,iter,5)
  colnames(cor.F) = c("Trans-ECM","Trans-SECM","Pool","Hard","sample")
  
  for (ii in 1:iter) {
    
    # data generation
    #-------------------------------------------
    data = dat.g(n0,n1,p,num.info,K,d,case,s = 10)
    Xt = data$Xt
    XS = data$XS
    
    Sig_t = data$Sig_t
    Sig_s = data$Sig_s
    
    diag.elements = diag(Sig_t)
    R_t = diag(diag.elements^{-1/2}) %*% Sig_t %*% diag(diag.elements^{-1/2})
    
    R_s = array(0,c(p,p,K))
    for (k in 1:K) {
      
      diag.elements = diag(Sig_s[,,k])
      R_s[,,k] = diag(diag.elements^{-1/2}) %*% Sig_s[,,k] %*% diag(diag.elements^{-1/2})
    }
    
    # pool
    #-----------------
    X.A = c()
    for (k in 1:K) {X.A = rbind(X.A, XS[,,k])}
    sum.X = rbind(Xt,X.A)
    R_hat_pool = cor(sum.X)
    # norm(R_hat_pool - R_t,"F")^2/p
    
    # Target-Hard
    #------------------
    R_hard_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    
    # Target-sample
    #------------------
    R_sample_hat = cor(Xt)
    
    # Trans-ECM 
    #------------------
    
    # initial estimators
    # R_t_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    R_t_hat = R_sample_hat 
    R_s_hat = R_s
    
    # computation of Trans-ECM
    A.range = seq(0.3,0.9,0.3)
    M.range = c(0.1,seq(0.5,2,0.5))
    R_Trans_ECM_hat_F = fit.trans.ECM.cv(Xt, XS, K, R_t, 
                                         R_t_hat, R_s_hat, A.range, M.range)
    
    # Trans-SECM 
    #-----------------
    R_Trans_SECM_hat_F = fit.trans.SECM(R_Trans_ECM_hat_F,Xt,R_t)
    
   
    
    #-------------------------------------
    cor.F[ii,] = c(
      norm(R_Trans_ECM_hat_F - R_t,"F")^2/p,
      norm(R_Trans_SECM_hat_F - R_t,"F")^2/p,
      norm(R_hat_pool - R_t,"F")^2/p,
      norm(R_hard_hat - R_t,"F")^2/p,
      norm(R_sample_hat - R_t,"F")^2/p)
  }
  
  cor.F = colMeans(cor.F)
  return(cor.F)
}

n0 = 100;n1 = 500;p = 100;K = 5;case = 2;iter = 100
num.info = 0;d = 0.08
b0 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)

num.info = 1
b1 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)

num.info = 2
b2 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)

num.info = 3
b3 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)

num.info = 4
b4 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)

num.info = 5
b5 = case2.fun.conditional(n0,n1,p,K,num.info,d,case,iter)
#--------------------------------------------------------


rbind(b0$cor.F,b1$cor.F,b2$cor.F,b3$cor.F,b4$cor.F,b5$cor.F)
rbind(b0$cor.2,b1$cor.2,b2$cor.2,b3$cor.2,b4$cor.2,b5$cor.2)

