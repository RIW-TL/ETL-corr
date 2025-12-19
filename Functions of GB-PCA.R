library(mvtnorm)
library(Matrix)
#---------------------------------------------------

dat.g.PCA = function(n0,n1,p,num.info,K,d){
  
  # rank
  r0 = 10;rk = 6
  rk_list = c(r0)
  for (k in 1:K) {rk_list[k + 1] = rk}
  
  # target covariance matrix
  U1 = matrix(rnorm(p*p),p,p)
  qr_result_U = qr(U1)
  U2 = qr.Q(qr_result_U)
  
  sigma0 = diag(c(rep(10,r0),rep(0.5,p - r0)))
  Sig_t = as.matrix(U2 %*% sigma0 %*% t(U2))
  
  # source covariance matrix array
  Sig_s = array(0,c(p,p,K))
  
  for (k in 1:K) {
    
    if (k <= num.info) {
      
      U1k = U1 + matrix(runif(p*p,0,d/p),p,p)
      qr_result_Uk = qr(U1k)
      U2k = qr.Q(qr_result_Uk)
      sigmak = sigma0
      
    }else{
      
      U1k = U1 + matrix(runif(p*p,0.5,1),p,p)
      qr_result_Uk = qr(U1k)
      U2k = qr.Q(qr_result_Uk)
      sigmak = diag(c(rep(10,rk),rep(0.5,p - rk)))
    }
    
    Sig_s[,,k]<- as.matrix(U2k %*% sigmak %*% t(U2k))
  }
  
  
  ## data generation
  XS = array(0,c(n1,p,K))
  mu_0 = mu_1 = rep(0,p)
  Xt = rmvnorm(n0,mu_0,Sig_t)
  XS = array(0,c(n1,p,K))
  for (k in 1:K) {
    
    XS[,,k] = rmvnorm(n1,mu_1,Sig_s[,,k])
  }
  
  return(list(Xt = Xt,XS = XS,
              Sig_t = Sig_t,Sig_s = Sig_s,
              rk_list = rk_list))
}

fit.evaluate = function(U,U_t){
  
  index = norm(U %*% t(U)  - U_t %*% t(U_t),"F")
  return(index)
}

#--------------------------------------------------------------
EP = function(Pt_list, nk_list) {
  
  K_a <- length(Pt_list)
  EP <- nk_list[1] * Pt_list[[1]]
  if(K_a>1){
    for (k in 2:K_a) {
      EP <- EP + nk_list[k] * Pt_list[[k]]
    }
  }
  return(EP)
}

# Pt_list=PA
# Step 1 of the Oracle algorithm: construct the Grassmannian barycenter--P_s
GB = function(Pt_list, nk_list, rs) {
  
  if(rs == 0){
    
    U = matrix(0,p,1)
    Ps = U %*% t(U)
  }else{
    
    EPA <- EP(Pt_list, nk_list)
    UG <- eigen(EPA/sum(nk_list))$vectors[,1:rs]
    Ps <- UG %*% t(UG)
  }
  
  return(Ps)
}

# Step 2 of the Oracle algorithm: finetune
finetune = function(Sigma, Ps, rp, p) {
  
  if(rp == 0){
    Uf = matrix(0,p,1)
  }else{
    
    Ps_perp <- diag(p) - Ps
    Uf <- eigen(Ps_perp %*% Sigma %*% Ps_perp)$vectors[, 1:rp] # rp = r0 - rs
  }
  
  return(Uf %*% t(Uf) + Ps)
}

# Better initial estimator
Ps_initial_1 = function(Pks, nk_list, K_a, n_select, rs) {
  
  PA <- list(Pks[[1]])
  k_selected <- c(1) # 一定将target的P0和对应指标放入进迭代初始值的估计之中
  
  dis0k <- numeric(K_a - 1) # K_a <- length(Sigma_list)
  for (k in 1:(K_a - 1)) {
    
    dis0k[k] <- sum(diag(Pks[[1]] %*% Pks[[k+1]])) # 计算trace
  }
  
  result <- sort(dis0k, decreasing = TRUE) # 降序排列
  tau0 <- result[n_select]
  
  for (k in 1:(K_a - 1)) {
    if (dis0k[k] > tau0) {
      
      PA[[length(PA) + 1]] <- Pks[[k+1]]
      k_selected <- c(k_selected, k+1)
    }
  }
  
  Ps <- GB(PA, nk_list[k_selected], rs)
  return(Ps)
}

GB_Kmeans = function(Sigma_list, nk_list, rk_list,rs, p, tau, 
                     Ps_initial = 1, n_select, Ti) { 
  
  
  K_a <- length(Sigma_list) # K_a表示source的个数 + 1
  Pks <- list()             # Pks表示所有的投影矩阵 Pk = Uk %*% t(Uk)
  
  for (k in 1:K_a) {
    
    Uk <- eigen(Sigma_list[[k]])$vectors[, 1:rk_list[k]] # r_list表示所有的秩
    Pks[[k]] <- Uk %*% t(Uk)
  }
  
  ## provide initial estimator
  # if (Ps_initial == 0) {
  #   
  #   Ps <- GB(Pks, nk_list, rs)
  # } else if (Ps_initial == 1) {
  #   
  #   Ps <- Ps_initial_1(Pks, nk_list, K_a, n_select, rs)
  # }
  
  uks = eigen(Sigma_list[[1]])$vectors[,1:rs]
  Ps = uks %*% t(uks)
  
  ## 迭代循环得到k_selected
  for (t in 1:Ti) {
    
    k_selected <- c(1)
    for (k in 2:length(Pks)) {
      
      if (sum(diag(Ps %*% Pks[[k]])) > tau) {
        k_selected <- c(k_selected, k)
      }
    }
    
    # k_selected包含target one
    Ps <- GB(Pks[k_selected], nk_list[k_selected], rs)
  }
  
  rp = rk_list[[1]] - rs
  Pksk_finetuned = finetune(Sigma_list[[1]], Ps, rp, p)
  
  list = list(k_selected = k_selected,Pksk_finetuned = Pksk_finetuned)
}

#---------------------------------------------------------------