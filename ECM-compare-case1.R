
case2.fun.compare = function(n0,n1,p,K,num.info,d,case,iter){
  
  cor.F = cov.F = cor.2 = cov.2 = matrix(0,iter,5)
  var.re = matrix(0,iter,3)
  
  colnames(cor.F) = colnames(cov.F) = colnames(cor.2) = colnames(cov.2) = 
    c("Trans-ECM","Trans-SECM","Pool","Hard","sample")
  colnames(var.re) = c("Trans","pool","sample")
  
  for (ii in 1:iter) {
    
    # data generation
    #-------------------------------------------
    data = dat.g(n0,n1,p,num.info,K,d,case,s = 100)
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
    
    #==================================cor============================
    #---------------------------------------------------------------
    
    # Target-sam
    #------------------
    R_sample_hat = cor(Xt)
    
    # Target-Hard
    #------------------
    # R_hard_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    R_hard_hat = R_sample_hat
    lambda.range = seq(0.02,0.5,0.02);
    R_hat = list();
    error_2 = error_F = c()
    
    for (l in seq_along(lambda.range)) {
      
      lambda = lambda.range[l]
      R_hard_hat[abs(R_hard_hat) < lambda] = 0
      R_hat[[l]] = R_hard_hat
      error_2[l] = norm(R_hard_hat - R_t,"2")
      error_F[l] = norm(R_hard_hat - R_t,"F")
    }
    R_hard_hat_2 = R_hat[[which.min(error_2)]]
    R_hard_hat_F = R_hat[[which.min(error_F)]]
    
    # Trans-ECM 
    #------------------
    
    # initial estimators
    # R_t_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    R_t_hat = R_hard_hat_F 
    R_s_hat = R_s
    
    # computation of Trans-ECM
    A.range = seq(0.3,0.9,0.3)
    M.range = c(0.05,seq(0.1,0.5,0.1))
    fit = fit.trans.ECM.cv(Xt, XS, K, R_t, R_t_hat, R_s_hat,A.range, M.range)
    R_Trans_ECM_hat_2 = fit$out.2
    R_Trans_ECM_hat_F = fit$out.F
    
    # Trans-SECM 
    #-----------------
    fit = fit.trans.SECM(R_Trans_ECM_hat_2,Xt,R_t)
    R_Trans_SECM_hat_2 = fit$out.2
    
    fit = fit.trans.SECM(R_Trans_ECM_hat_F,Xt,R_t)
    R_Trans_SECM_hat_F = fit$out.F
    
    # pool
    #-----------------
    X.A = c()
    for (k in 1:K) {X.A = rbind(X.A, XS[,,k])}
    sum.X = rbind(Xt,X.A)
    R_hat_pool = cor(sum.X)
    
    #-------------------------------------
    cor.F[ii,] = c(
      norm(R_Trans_ECM_hat_F - R_t,"F")^2/p,
      norm(R_Trans_SECM_hat_F - R_t,"F")^2/p,
      norm(R_hat_pool - R_t,"F")^2/p,
      norm(R_hard_hat_F - R_t,"F")^2/p,
      norm(R_sample_hat - R_t,"F")^2/p)
    
    cor.2[ii,] = c(
      norm(R_Trans_ECM_hat_2 - R_t,"2"),
      norm(R_Trans_SECM_hat_2 - R_t,"2"),
      norm(R_hat_pool - R_t,"2"),
      norm(R_hard_hat_2 - R_t,"2"),
      norm(R_sample_hat - R_t,"2"))
    
    
    #==========================var==============================
    #===========================================================
    
    # Target-only
    #------------------------------
    target.var.hat = diag(cov(Xt))
    
    # Pool
    #------------------------------
    pool.var.hat = diag(cov(sum.X))
    
    # Trans
    #-------------------------------
    A.range = seq(0.3,0.9,0.3)
    M.range = seq(0.2,0.5,0.1)
    trans.var.hat = fit.2rd.moment.trans.cv(Xt,XS,K,Sig_t,A.range,M.range)
    
    
    #-------------------------------
    SSE.target = sum((target.var.hat - diag(Sig_t))^2)
    SSE.pool = sum((diag(cov(X.A)) - diag(Sig_t))^2)
    SSE.trans = sum((trans.var.hat - diag(Sig_t))^2)
    
    var.re[ii,] = c(SSE.trans,SSE.pool,SSE.target)
    
    #==========================cov==============================
    #===========================================================
    
    # Target-sample
    #------------------
    cov_sample_hat = cov(Xt)
    
    # Target-Hard
    #------------------
    # R_hard_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)
    cov_hard_hat = cov_sample_hat
    lambda.range = seq(0.02,0.5,0.02);
    cov_hat = list();
    error_2 = error_F = c()
    
    for (l in seq_along(lambda.range)) {
      
      lambda = lambda.range[l]
      cov_hard_hat[abs(cov_hard_hat) < lambda] = 0
      cov_hat[[l]] = cov_hard_hat
      error_2[l] = norm(cov_hard_hat - Sig_t,"2")
      error_F[l] = norm(cov_hard_hat - Sig_t,"F")
    }
    cov_hard_hat_2 = cov_hat[[which.min(error_2)]]
    cov_hard_hat_F = cov_hat[[which.min(error_F)]]
    
    
    # Trans-ECM 
    #------------------
    cov_Trans_ECM_hat_F = diag(trans.var.hat^{1/2}) %*% R_Trans_ECM_hat_F %*% diag(trans.var.hat^{1/2})
    cov_Trans_ECM_hat_2 = diag(trans.var.hat^{1/2}) %*% R_Trans_ECM_hat_2 %*% diag(trans.var.hat^{1/2})
    
    # Trans-SECM 
    #-----------------
    cov_Trans_SECM_hat_F = diag(trans.var.hat^{1/2}) %*% R_Trans_SECM_hat_F %*% diag(trans.var.hat^{1/2})
    cov_Trans_SECM_hat_2 = diag(trans.var.hat^{1/2}) %*% R_Trans_SECM_hat_2 %*% diag(trans.var.hat^{1/2})
    
    # pool
    #-----------------
    cov_hat_pool = cov(sum.X)
    
    #------------------------------------------
    
    cov.F[ii,] = c(
      norm(cov_Trans_ECM_hat_F - Sig_t,"F")^2/p,
      norm(cov_Trans_SECM_hat_F - Sig_t,"F")^2/p,
      norm(cov_hat_pool - Sig_t,"F")^2/p,
      norm(cov_hard_hat_F - Sig_t,"F")^2/p,
      norm(cov_sample_hat - Sig_t,"F")^2/p)
    
    cov.2[ii,] = c(
      norm(cov_Trans_ECM_hat_2 - Sig_t,"2"),
      norm(cov_Trans_SECM_hat_2 - Sig_t,"2"),
      norm(cov_hat_pool - Sig_t,"2"),
      norm(cov_hard_hat_2 - Sig_t,"2"),
      norm(cov_sample_hat - Sig_t,"2"))
  }
  
  return(list(cor.2 = colMeans(cor.2),
              cov.2 = colMeans(cov.2),
              cor.F = colMeans(cor.F),
              cov.F = colMeans(cov.F),
              var.re = colMeans(var.re)))
  
}

n0 = 100;n1 = 500;p = 100;K = 5;case = 1;iter = 50
num.info = 0;d = 0.04
a0 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 1
a1 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 2
a2 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 3
a3 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 4
a4 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 5
a5 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

#-----------------------------------------------------

n0 = 100;n1 = 500;p = 100;K = 5
num.info = 0;d = 0.08
b0 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 1
b1 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 2
b2 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 3
b3 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 4
b4 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 5
b5 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
#--------------------------------------------------------

n0 = 100;n1 = 500;p = 100;K = 5
num.info = 0;d = 0.12
c0 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 1
c1 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 2
c2 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 3
c3 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 4
c4 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

num.info = 5
c5 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

#--------------------------------------------
# n0 = 100;n1 = 500;p = 100;K = 5
# num.info = 0;d = 0.2
# d0 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
# 
# num.info = 1
# d1 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
# 
# num.info = 2
# d2 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
# 
# num.info = 3
# d3 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
# 
# num.info = 4
# d4 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)
# 
# num.info = 5
# d5 = case2.fun.compare(n0,n1,p,K,num.info,d,case,iter)

rbind(a0$cor.F,a1$cor.F,a2$cor.F,a3$cor.F,a4$cor.F,a5$cor.F)
rbind(b0$cor.F,b1$cor.F,b2$cor.F,b3$cor.F,b4$cor.F,b5$cor.F)
rbind(c0$cor.F,c1$cor.F,c2$cor.F,c3$cor.F,c4$cor.F,c5$cor.F)

rbind(a0$cor.2,a1$cor.2,a2$cor.2,a3$cor.2,a4$cor.2,a5$cor.2)
rbind(b0$cor.2,b1$cor.2,b2$cor.2,b3$cor.2,b4$cor.2,b5$cor.2)
rbind(c0$cor.2,c1$cor.2,c2$cor.2,c3$cor.2,c4$cor.2,c5$cor.2)
#----------------------------------------------------------------
rbind(a0$cov.F,a1$cov.F,a2$cov.F,a3$cov.F,a4$cov.F,a5$cov.F)
rbind(b0$cov.F,b1$cov.F,b2$cov.F,b3$cov.F,b4$cov.F,b5$cov.F)
rbind(c0$cov.F,c1$cov.F,c2$cov.F,c3$cov.F,c4$cov.F,c5$cov.F)

rbind(a0$cov.2,a1$cov.2,a2$cov.2,a3$cov.2,a4$cov.2,a5$cov.2)
rbind(b0$cov.2,b1$cov.2,b2$cov.2,b3$cov.2,b4$cov.2,b5$cov.2)
rbind(c0$cov.2,c1$cov.2,c2$cov.2,c3$cov.2,c4$cov.2,c5$cov.2)
#----------------------------------------------------------------
rbind(a0$var.re,a1$var.re,a2$var.re,a3$var.re,a4$var.re,a5$var.re)
rbind(b0$var.re,b1$var.re,b2$var.re,b3$var.re,b4$var.re,b5$var.re)
rbind(c0$var.re,c1$var.re,c2$var.re,c3$var.re,c4$var.re,c5$var.re)
