library(mvtnorm)
library(ks)
library(Matrix)
library(sparseCov)
library(moments)
#----------------------------------------------------

fit.def.sym = function(Sig){
  
  p = dim(Sig)[2]
  eigen_decomp = eigen(Sig)
  eigenvalues = eigen_decomp$values
  
  if(min(eigenvalues) >= 0.1){
    
    Sig = Sig
  }else{
    
    c = 0.1 - min(eigenvalues)
    Sig = Sig + c*diag(1,p)
  }
  
  Sig[abs(Sig)< 10^{-10}] = 0
  
  Sig = (Sig + t(Sig))/2
  return(Sig)
}

# 注意normal and t distributions
dat.g = function(n0,n1,p,num.info,K,d,case,s){
  
  ## parameter generation
  ##---------------------------------------
  
  # target covariance matrix
  if(case == 1){
    
    Sig_t = diag(1,p)
    
    for (j in 1:(p - 1)) {
      for (l in (j + 1):p) {
        
        Sig_t[j,l] = 0.7^{abs(j - l)}
        Sig_t[l,j] = Sig_t[j,l]
      }
    }
    
  }else{
    
    # target precision matrix
    B = matrix(0, nrow = p, ncol = p)
    for (i in 1:p) {
      for (j in i:p) {
        
        if(i == j){
          
          B[i, j] = 1 + runif(1, 0, 0.2)
          
        }else {
          
          value = runif(1, 0, 0.2)
          B[i, j] <- value
          B[j, i] <- value
        }
      }
    }
    
    B = (B + t(B))/2
    
    # sparse covariance matrix
    if(s < p){
      
      for(j in 1:p){
        
        B[,j]<-B[,j]*(abs(B[,j])>=quantile(abs(B[,j]),1-s/p))
        B[j,]<-B[j,]*(abs(B[j,])>=quantile(abs(B[j,]),1-s/p))
      }
    }
    
    
    Sig_t = fit.def.sym(B) 
  }
  
  # ensure the positive definiteness and symmetry
  Sig_t = fit.def.sym(Sig_t) 
  
  # source covariance matrices
  Sig_s = array(0,c(p,p,K))
  
  fit = eigen(Sig_t)
  U = fit$vectors
  Lambda.T = fit$values
  
  for (k in 1:K) {
    
    if(k <= num.info){
      
      # Lambda.S = Lambda.T + runif(p,min = 0,max = d)
      Lambda.S = Lambda.T
      
      US1 = U + matrix(runif(p^2,0,d),p)
      qr_result = qr(US1)
      US1 = qr.Q(qr_result)
      # US1 = U
      
      Sig_s[,,k] = US1 %*% diag(Lambda.S) %*% t(US1)
    }else{
      
      Lambda.S = Lambda.T + runif(p,min = 1,max = 2)
      
      US2 = U + matrix(runif(p^2,0.5,1),p)
      qr_result = qr(US2)
      US2 = qr.Q(qr_result)
      
      Sig_s[,,k] = US2 %*% diag(Lambda.S) %*% t(US2)
    }
    
    Sig_s[,,k] = (Sig_s[,,k] + t(Sig_s[,,k]))/2
    
  }
  
  ## data generation
  ##-----------------------------------------------
  # Xt = rmvnorm(n0,rep(0,p),Sig_t)
  # Xt = rmvt(n0,Sig_t,df = 5)
  
  # XS = array(0,c(n1,p,K))
  # for (k in 1:K) {XS[,,k] = rmvnorm(n1,rep(0,p),Sig_s[,,k])}
  # for (k in 1:K) {XS[,,k] = rmvt(n1,Sig_s[,,k],df = 5)}

  Zt = matrix(0,n0,p)
  
  Zt[,1] = rchisq(n0,df = 5)
  Zt[,2] = rt(n0,df = 5)
  Zt[,3] = runif(n0,min = -1,max = 1)
  
  for(j in 4:p){Zt[,j] = rnorm(n0)}
  
  eig <- eigen(Sig_t)
  Q <- eig$vectors
  Lambda <- diag(eig$values)
  Lambda_sqrt <- diag(sqrt(eig$values))
  Sig_t_sqrt <- Q %*% Lambda_sqrt %*% t(Q)
  
  Xt = Zt %*% Sig_t_sqrt
  
 #----------------------------------
  XS = array(0,c(n1,p,K))
  ZS = array(0,c(n1,p,K))
  
  for (k in 1:K) {
    
    ZS[,1,k] = rchisq(n1,df = 5)
    ZS[,2,k] = rt(n1,df = 5)
    ZS[,3,k] = runif(n1,min = -1,max = 1)
    
    for(j in 4:p){ZS[,j,k] = rnorm(n1)}
    
    eig <- eigen(Sig_s[,,k])
    Q <- eig$vectors
    Lambda <- diag(eig$values)
    Lambda_sqrt <- diag(sqrt(eig$values))
    Sig_k_sqrt <- Q %*% Lambda_sqrt %*% t(Q)
    
    XS[,,k] = ZS[,,k] %*% Sig_k_sqrt
  }

  
  return(list(Xt = Xt, XS = XS, 
              Sig_t = Sig_t, Sig_s = Sig_s))
}

#---------------------------------------------------------------------
# 辅助函数：预计算所有与c无关的中间变量

# 注意计算one or two-dimensional density 
precompute_c_independent = function(Yt, YS, K, R_t_hat, R_s_hat) {
  
  p <- dim(Yt)[2]
  n1 <- dim(YS)[1]
  
  # 存储预计算结果的列表
  pre <- list()
  
  # 预计算目标域的交叉项（与c无关）
  pre$Yt_cross <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in 1:p) {
      pre$Yt_cross[j, l] <- sum(Yt[, j] * Yt[, l])
    }
  }
  
  # 预计算每个(j,l,k)的误差项、KDE结果和权重
  pre$errors <- array(list(), dim = c(p, p, K))  # 存储error_0, error_1, eta
  pre$weights <- array(0, dim = c(n1, K, p, p))  # 存储weight[,k]
  # pre$A <- array(0, dim = c(p, p, K))  # 存储A（与c无关）
  # pre$h_jl <- array(0, dim = c(p, p, K))  # 存储h_jl
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      for (k in 1:K) {
        
        # 误差项和eta（与c无关）
        eta <- (R_t_hat[j, l] - R_s_hat[j, l, k]) * YS[, l, k]
        error_0 <- YS[, j, k] - R_t_hat[j, l] * YS[, l, k]
        error_1 <- YS[, j, k] - R_s_hat[j, l, k] * YS[, l, k]
        pre$errors[[j, l, k]] <- list(eta = eta, error_0 = error_0, error_1 = error_1)
        
        # KDE和权重计算（与c无关）
        kde.fit <- kde(x = error_1)
        fenzi <- (predict(kde.fit, x = error_0) + predict(kde.fit, x = -error_0)) / 2
        # # fenzi <- predict(kde.fit, x = error_0) 
        fenmu <- predict(kde.fit, x = error_1)
        
        #-------------------------------------------------------
        # fenzi = dnorm(error_0)
        # 
        # # 计算条件密度估计
        # if(j <= 3){
        #   
        # fit_joint <- kde(x = YS[, c(j,l),k])
        # fenmu.fenzi <- predict(fit_joint, x = YS[, c(j,l), k])
        # }else{
        #   
        #   fenmu.fenzi = dmvnorm(x = YS[,c(j,l),k], 
        #                         mean = rep(0,2), sigma = R_s_hat[c(j,l),c(j,l),k],
        #                         log = FALSE)
        # }
        # 
        # # fenmu.fenzi = dmvt(x = YS[,c(j,l),k],sigma = R_s_hat[c(j,l),c(j,l),k],
        # #                    df = 5,log = FALSE)
        # fit_marginal <- kde(x = YS[, l, k])        # 边缘密度估计（X2）
        # fenmu.fenmu <- predict(fit_marginal, x = YS[, l, k])
        # fenmu = fenmu.fenzi/fenmu.fenmu
        #-------------------------------------------------------
        
        pre$weights[, k, j, l] <- fenzi / fenmu  # 存储权重
      }
    }
  }
  
  return(pre)
}

# 优化后的主函数：仅计算与c相关的部分
fit.trans.ECM = function(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre) {
  
  p <- dim(Yt)[2]
  R_trans <- diag(1, p)
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      fenzi.sum <- fenmu.sum <- rep(0, K)
      
      for (k in 1:K) {
        # 直接复用预计算的结果
        err <- pre$errors[[j, l, k]]
        weight <- pre$weights[, k, j, l]
        # A <- pre$A[j, l, k]
        # rho <- 0.7
        # A <- quantile(abs(err$error_1), rho)
        # M <- quantile(abs(err$eta), rho.M)
        
        # boxplot(err$eta,main = "eta")
        # boxplot(err$error_1,main = "error1")
        # boxplot(weight,main = "weight")
        # l = l + 1
        
        # 仅与c相关的计算：M和样本选择
        # M <- c / max(h_jl, 0.001)
        # A = 0.5;M = 0.05
        select.id1 <- which(abs(err$error_1) <= A)
        select.id2 <- which(abs(err$eta) <= M)
        id <- intersect(select.id1, select.id2)
        
        # boxplot(err$eta[id],main = "eta[id]")
        # boxplot(err$error_1[id],main = "error1[id]")
        # boxplot(weight[id],main = "weight[id]")
        
        # 计算分子和分母（复用预计算的权重和YS）
        fenzi.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, j, k])
        fenmu.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, l, k])
      }
      
      # 复用预计算的目标域交叉项
      fenzi <- sum(fenzi.sum) + pre$Yt_cross[j, l]
      fenmu <- sum(fenmu.sum) + pre$Yt_cross[l, l]
      R_trans[j, l] <- fenzi / fenmu
      R_trans[l, j] <- R_trans[j, l]
    }
  }
  
  return(R_trans)
}

# 优化后的交叉验证函数：整合预计算步骤
# num.info
fit.trans.ECM.cv = function(Xt, XS, K, R_t, 
                            R_t_hat, R_s_hat, A.range, M.range) {
  # 1. 基础计算（标准化和初始估计）
  n0 <- dim(Xt)[1]
  p <- dim(Xt)[2]
  n1 <- dim(XS[,,1])[1]
  Yt <- scale(Xt)
  YS <- array(0, c(n1, p, K))
  for (k in 1:K) {
    YS[,,k] <- scale(XS[,,k])
  }
  
  # 2. 预计算所有与c无关的变量（核心优化）
  pre <- precompute_c_independent(Yt, YS, K, R_t_hat, R_s_hat)
  
  cv.loss.F = matrix(0, length(A.range), length(M.range))
  out = list()
  
  for (ir in seq_along(A.range)) {
    for (im in seq_along(M.range)) {
      
      A = A.range[ir]
      M = M.range[im]
      
      l = (ir - 1)*length(M.range) + im
      out[[l]] = fit.trans.ECM(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre)
      cv.loss.F[ir,im] <- norm(out[[l]] - R_t, "F")
    }
  }
  
  # ic.opt.F <- which.min(cv.loss.F)
  # ic.opt.2 <- which.min(cv.loss.2)
  
  min.id = which(cv.loss.F == min(cv.loss.F),arr.ind = TRUE)
  l.opt = (min.id[1,1] - 1)*length(M.range) + min.id[1,2]
  out.F = out[[l.opt]]
  
  # min.id = which(cv.loss.2 == min(cv.loss.2),arr.ind = TRUE)
  # l.opt = (min.id[1,1] - 1)*length(M.range) + min.id[1,2]
  # out.2 = out[[l.opt]]
  
  return(out.F)
  
}
#---------------------------------------------------------------------
fit.trans.SECM = function(R_trans,Xt,R_t){
  
  lambda.range = seq(0.01,0.1,0.02) # can be change in practice
  R_strans = R_trans
  
  cv.SECM = list()
  cv.loss.2 = cv.loss.F = rep(0,length(lambda.range))
  
  for (m in 1:length(lambda.range)) {
    
    lambda = lambda.range[m]
    
    for (j in 1:(p - 1)) {
      for (l in (j + 1):p) {
        
        if(abs(R_strans[j,l]) < lambda){R_strans[j,l] = R_strans[l,j] = 0}
        
      }
    }
    
    cv.SECM[[m]] = R_strans
    # cv.loss.2[m] = norm(cv.SECM[[m]] - R_t,"2")
    cv.loss.F[m] = norm(cv.SECM[[m]] - R_t,"F")
  }
  
  ic.opt.F <- which.min(cv.loss.F)
  # ic.opt.2 <- which.min(cv.loss.2)
  out.F = cv.SECM[[ic.opt.F]]
  
  return(out.F)
}

fit.1st.moment.trans = function(Xt,XS,K,A,M){
  
  n0 = nrow(Xt)
  n1 = nrow(XS[,,1])
  p = ncol(Xt)
  
  # initial estimators
  mean.t.hat = colMeans(Xt)
  mean.s.hat = matrix(0,p,K)
  for (k in 1:K) {mean.s.hat[,k] = colMeans(XS[,,k])}
  
  # Transfer-enhanced mean estimation
  mean.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = mean.t.hat[j] - mean.s.hat[j,]
    kp.id = which(abs(eta) <= M)
    # kp.id = 1:num.info
    
    if(length(kp.id) == 0){
      
      mean.trans[j] = colMeans(Xt)[j]
      
    }else{
      
      if(length(kp.id) == 1){
        
        XS.use = array(0,c(n1,p,1))
        XS.use[,,1] = XS[,,1]
      }else{
        
        XS.use = XS[,,kp.id]
        
        fenzi.sum = fenmu.sum = rep(0,length(kp.id))
        
        for (k in 1:length(kp.id)) {
          
          # the estimation of density functions
          error_0 = XS.use[,j,k] - mean.t.hat[j]
          error_1 = XS.use[,j,k] - mean.s.hat[j,k]
          kde.fit = kde(x = error_1)
          
          # the estimation of importance weights
          fenzi = (predict(kde.fit,x = error_0) +
                     predict(kde.fit,x = -error_0))/2
          fenmu = predict(kde.fit,x = error_1)
          weight = fenzi/fenmu
          
          # sample selection
          id = which(abs(error_1) <= A)
          
          # weights estimator (mean estimation)
          fenzi.sum[k] = sum(XS.use[id,j,k]*weight[id])
          fenmu.sum[k] = sum(weight[id])
        }
        
        mean.trans[j] = (sum(Xt[,j]) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
      }
    }
    
  }
  
  return(mean.trans)
}

fit.2rd.moment.trans = function(Xt,XS,K,A,M,m2.t.hat,m2.s.hat){
  
  n0 = nrow(Xt)
  n1 = nrow(XS[,,1])
  p = ncol(Xt)
  
  # transfer-enhanced (2rd) moment estimation
  m2.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = m2.t.hat[j] - m2.s.hat[j,]
    kp.id = which(abs(eta) <= M)
    # kp.id = 1:num.info
    
    if(length(kp.id) == 0){
      
      m2.trans[j] = m2.t.hat[j]
      
    }else{
      
      if(length(kp.id) == 1){
        
        XS.use = array(0,c(n1,p,1))
        XS.use[,,1] = XS[,,1]
      }else{
        
        XS.use = XS[,,kp.id]
      }
      
      fenzi.sum = fenmu.sum = rep(0,length(kp.id))
      
      for (k in 1:length(kp.id)) {
        
        # the estimation of density functions
        error_0 = (XS.use[,j,k])^2 - m2.t.hat[j]
        error_1 = (XS.use[,j,k])^2 - m2.s.hat[j,k]
        kde.fit = kde(x = error_1)
        
        # the estimation of importance weights
        fenzi = (predict(kde.fit,x = error_0) +
                   predict(kde.fit,x = -error_0))/2
        fenmu = predict(kde.fit,x = error_1)
        weight = fenzi/fenmu
        
        # sample selection
        id = which(abs(error_1) <= A)
        
        # boxplot(weight)
        # boxplot(error_1)
        # boxplot(weight[id])
        # boxplot(error_1[id])
        
        # weights estimator (mean estimation)
        fenzi.sum[k] = sum(((XS.use[id,j,k])^2)*weight[id])
        fenmu.sum[k] = sum(weight[id])
      }
      
      m2.trans[j] = (sum((Xt[,j])^2) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
    }
    
  }
  
  return(m2.trans)
}

fit.2rd.moment.trans.cv = function(Xt,XS,K,Sig_t, A.range, M.range){
  
  # initial estimates
  m2.t.hat = moment(Xt,order = 2)
  
  m2.s.hat = matrix(0,p,K)
  for (k in 1:K) {m2.s.hat[,k] = moment(XS[,,k],order = 2)}
  # for (k in 1:K) {m2.s.hat[,k] = diag(Sig_s[,,k])}
  
  cv.loss = matrix(0,length(A.range),length(M.range))
  
  for (a in 1:length(A.range)) {
    for (m in 1:length(M.range)) {
      
      A = A.range[a]
      M = M.range[m]
      
      m2.trans.hat = fit.2rd.moment.trans(Xt,XS,K,A,M,m2.t.hat,m2.s.hat)
      
      cv.loss[a,m] = sum((m2.trans.hat - diag(Sig_t))^2)
    }
  }
  
  id = which(cv.loss == min(cv.loss), arr.ind = TRUE)
  A.opt = A.range[id[1,1]]
  M.opt = M.range[id[1,2]]
  
  m2.trans.hat = fit.2rd.moment.trans(Xt,XS,K,A.opt,M.opt,m2.t.hat,m2.s.hat)
  return(m2.trans.hat)
}


