library(mvtnorm)
library(ks)
library(Matrix)
library(sparseCov)
library(moments)
#----------------------------------------------------

# symmetrization and positive definitization
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

## data generation 
# num.info: the size of informative sources
# d: informative level
# s: sparsity level
dat.g = function(n0,n1,p,num.info,K,d,case,s){
  
  ## parameter generation
  ##---------------------------------------
  
  # target covariance matrix (two settings)
  if(case == 1){
    
    Sig_t = diag(1,p)
    
    for (j in 1:(p - 1)) {
      for (l in (j + 1):p) {
        
        Sig_t[j,l] = 0.7^{abs(j - l)}
        Sig_t[l,j] = Sig_t[j,l]
      }
    }
    
    Sig_t = fit.def.sym(Sig_t) 
    
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
    
    # generate a sparse covariance matrix. Each column remains s elements. 
    if(s < p){
      
      for(j in 1:p){
        
        B[,j]<-B[,j]*(abs(B[,j])>=quantile(abs(B[,j]),1-s/p))
        B[j,]<-B[j,]*(abs(B[j,])>=quantile(abs(B[j,]),1-s/p))
      }
    }
    
    
    Sig_t = fit.def.sym(B) 
  }
  
  # source covariance matrices
  Sig_s = array(0,c(p,p,K))
  
  fit = eigen(Sig_t)
  U = fit$vectors
  Lambda.T = fit$values
  
  for (k in 1:K) {
    
    if(k <= num.info){
      
      Lambda.S = Lambda.T + runif(p,min = 0,max = 0.2)
      
      US1 = U + matrix(runif(p^2,0,d/p),p)
      qr_result = qr(US1)
      US1 = qr.Q(qr_result)
      
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
  
  # Gaussian 
  Xt = rmvnorm(n0,rep(0,p),Sig_t)
  XS = array(0,c(n1,p,K))
  for (k in 1:K) {XS[,,k] = rmvnorm(n1,rep(0,p),Sig_s[,,k])}
  
  # # non-Guassian
  # Zt = matrix(0,n0,p)
  # 
  # Zt[,1] = rchisq(n0,df = 5)
  # Zt[,2] = rt(n0,df = 5)
  # Zt[,3] = runif(n0,min = -1,max = 1)
  # 
  # for(j in 4:p){Zt[,j] = rnorm(n0)}
  # 
  # eig <- eigen(Sig_t)
  # Q <- eig$vectors
  # Lambda <- diag(eig$values)
  # Lambda_sqrt <- diag(sqrt(eig$values))
  # Sig_t_sqrt <- Q %*% Lambda_sqrt %*% t(Q)
  # 
  # Xt = Zt %*% Sig_t_sqrt
  
  # XS = array(0,c(n1,p,K))
  # ZS = array(0,c(n1,p,K))
  # 
  # for (k in 1:K) {
  #   
  #   ZS[,1,k] = rchisq(n1,df = 5)
  #   ZS[,2,k] = rt(n1,df = 5)
  #   ZS[,3,k] = runif(n1,min = -1,max = 1)
  #   
  #   for(j in 4:p){ZS[,j,k] = rnorm(n1)}
  #   
  #   eig <- eigen(Sig_s[,,k])
  #   Q <- eig$vectors
  #   Lambda <- diag(eig$values)
  #   Lambda_sqrt <- diag(sqrt(eig$values))
  #   Sig_k_sqrt <- Q %*% Lambda_sqrt %*% t(Q)
  #   
  #   XS[,,k] = ZS[,,k] %*% Sig_k_sqrt
  # }

  return(list(Xt = Xt, XS = XS, 
              Sig_t = Sig_t, Sig_s = Sig_s))
}

# Main function in implementing ETL-corr
pre_compute = function(Xt, XS, R_t_hat, R_s_hat) {
  
  p = dim(Xt)[2]
  n1 = dim(XS)[1]
  K = dim(XS)[3]
  
  # Step 1: data normaliation  
  #---------------------------------------
  Yt <- scale(Xt)
  YS <- array(0, c(n1, p, K))
  for (k in 1:K) {
    YS[,,k] <- scale(XS[,,k])
  }
  
  # List for storing precomputed results
  pre <- list()
  
  pre$Yt_cross <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in 1:p) {
      pre$Yt_cross[j, l] <- sum(Yt[, j] * Yt[, l])
    }
  }
  
  # Precompute the error term, KDE result and weight for each (j, l, k)
  pre$errors <- array(list(), dim = c(p, p, K))  
  pre$weights <- array(0, dim = c(n1, K, p, p)) 
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      for (k in 1:K) {
        
        eta <- (R_t_hat[j, l] - R_s_hat[j, l, k]) * YS[, l, k]
        error_0 <- YS[, j, k] - R_t_hat[j, l] * YS[, l, k]
        error_1 <- YS[, j, k] - R_s_hat[j, l, k] * YS[, l, k]
        pre$errors[[j, l, k]] <- list(eta = eta, error_0 = error_0, error_1 = error_1)
        
        # kernel density estimation (one-dimensional density ratio)
        kde.fit <- kde(x = error_1)
        fenzi <- (predict(kde.fit, x = error_0) + predict(kde.fit, x = -error_0)) / 2
        fenmu <- predict(kde.fit, x = error_1)
        
        # # kernel density estimation (two-dimensional density ratio)
        # fenzi = dnorm(error_0)
        # 
        # fit_joint <- kde(x = YS[, c(j,l),k])
        # fenmu.fenzi <- predict(fit_joint, x = YS[, c(j,l), k])
        # 
        # fit_marginal <- kde(x = YS[, l, k])
        # fenmu.fenmu <- predict(fit_marginal, x = YS[, l, k])
        # fenmu = fenmu.fenzi/fenmu.fenmu

        
        pre$weights[, k, j, l] <- fenzi / fenmu  
      }
    }
  }
  
  return(pre)
}

fit.ETL.corr = function(Xt, XS, A, M, R_t_hat, R_s_hat, pre) {
  
  n0 <- dim(Xt)[1]
  p <- dim(Xt)[2]
  n1 <- dim(XS[,,1])[1]
  K = dim(XS)[3]
  
  # Step 1: data normaliation  
  #---------------------------------------
  Yt <- scale(Xt)
  YS <- array(0, c(n1, p, K))
  for (k in 1:K) {
    YS[,,k] <- scale(XS[,,k])
  }
  
  # Step 2: kernel density estimation: pre
  # Step 3: weights estimation: pre

  R_trans <- diag(1, p)
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      fenzi.sum <- fenmu.sum <- rep(0, K)
      
      for (k in 1:K) {

        err <- pre$errors[[j, l, k]]
        weight <- pre$weights[, k, j, l]
        
        # sample selection
        select.id1 <- which(abs(err$error_1) <= A)
        select.id2 <- which(abs(err$eta) <= M)
        id <- intersect(select.id1, select.id2)
        
        fenzi.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, j, k])
        fenmu.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, l, k])
      }
      
      # ETL-corr estimator construction
      fenzi <- sum(fenzi.sum) + pre$Yt_cross[j, l]
      fenmu <- sum(fenmu.sum) + pre$Yt_cross[l, l]
      R_trans[j, l] <- fenzi / fenmu
      R_trans[l, j] <- R_trans[j, l]
    }
  }

  R_trans = fit.def.sym(R_trans)
  return(R_trans)
}

# ETL-Scor via hard-thresholding
fit.ETL.Scor = function(R_trans,lambda){
  
  R_strans = R_trans
    
    for (j in 1:(p - 1)) {
      for (l in (j + 1):p) {
        
        if(abs(R_strans[j,l]) < lambda){R_strans[j,l] = R_strans[l,j] = 0}
        
      }
    }
  
  return(R_strans)
  }

# Tranfer-enhanced variance estimation
fit.1st.moment.trans = function(Xt,XS,A1,M1){
  
  n0 = nrow(Xt)
  n1 = nrow(XS[,,1])
  p = ncol(Xt)
  K = dim(XS)[3]
  
  # initial estimators
  mean.t.hat = colMeans(Xt)
  mean.s.hat = matrix(0,p,K)
  for (k in 1:K) {mean.s.hat[,k] = colMeans(XS[,,k])}
  
  # Transfer-enhanced mean estimation
  mean.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = mean.t.hat[j] - mean.s.hat[j,]
    kp.id = which(abs(eta) <= M1)
    
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
          id = which(abs(error_1) <= A1)
          
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

fit.2rd.moment.trans = function(Xt,XS,A2,M2){
  
  n0 = nrow(Xt)
  n1 = nrow(XS[,,1])
  p = ncol(Xt)
  K = dim(XS)[3]
  
  m2.t.hat = moment(Xt,order = 2)
  m2.s.hat = matrix(0,p,K)

  for (k in 1:K) {m2.s.hat[,k] = moment(XS[,,k],order = 2)}
  
  # transfer-enhanced (2rd) moment estimation
  m2.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = m2.t.hat[j] - m2.s.hat[j,]
    kp.id = which(abs(eta) <= M2)
    
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
        id = which(abs(error_1) <= A2)
        
        # weights estimator (mean estimation)
        fenzi.sum[k] = sum(((XS.use[id,j,k])^2)*weight[id])
        fenmu.sum[k] = sum(weight[id])
      }
      
      m2.trans[j] = (sum((Xt[,j])^2) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
    }
    
  }
  
  return(m2.trans)
}

fit.ETL.var = function(Xt,XS,A1,M1,A2,M2){
  
  m1.trans.hat = fit.1st.moment.trans(Xt,XS,A1,M1)
  m2.trans.hat = fit.2rd.moment.trans(Xt,XS,A2,M2)
  
  ETL.var.hat = m2.trans.hat - m1.trans.hat^2
  
  return(ETL.var.hat)
}




