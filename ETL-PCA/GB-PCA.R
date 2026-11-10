#--------------------------------------------------------------
#输入一个矩阵列表 并将所有列表中的矩阵加权求和 输出一个矩阵
EP <- function(Pt_list, nk_list) {
  
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
#oracle算法的第一步：构造出Grassmannian barycenter Ps
GB <- function(Pt_list, nk_list, rs) {
  
  if(rs == 0){
    
    U = matrix(0,p,1)
    Ps = U %*% t(U)
  }else{
    
    EPA <- EP(Pt_list, nk_list)
    UG <- eigen(EPA / sum(nk_list))$vectors[, 1:rs]
    Ps <- UG %*% t(UG)
  }
  
  return(Ps)
}

#oracle算法的第二步finetune
finetune <- function(Sigma, Ps, rp, p) {
  
  if(rp == 0){
    Uf = matrix(0,p,1)
  }else{
    
    Ps_perp <- diag(p) - Ps
    Uf <- eigen(Ps_perp %*% Sigma %*% Ps_perp)$vectors[, 1:rp] # rp = r0 - rs
  }
  
  return(Uf %*% t(Uf) + Ps)
}

GB_Kmeans = function(Sigma_list, nk_list, rk_list,rs, p, tau, Ti) { 
  
  
  K_a <- length(Sigma_list) # K_a表示source的个数 + 1
  Pks <- list()             # Pks表示所有的投影矩阵 Pk = Uk %*% t(Uk)
  
  for (k in 1:K_a) {
    
    Uk <- eigen(Sigma_list[[k]])$vectors[, 1:rk_list[[k]]] # r_list表示所有的秩
    Pks[[k]] <- Uk %*% t(Uk)
  }
  
  uks = eigen(Sigma_list[[1]])$vectors[, 1:rs]
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

fit_GB = function(X0,XS,r0, rk_list,
                  rs_ratio = c(0.25,0.5,0.75)){
  
  p = ncol(X0)
  
  Sigma_list = list()
  nk_list = c()
  Sigma_list[[1]] = cov(X0);nk_list[1] = nrow(X0)
  
  for(k in 1:K){
    
    Sigma_list[[k+1]] = cov(XS[[k]])
    nk_list[k+1] = nrow(XS[[k]])
  }
  
  rs = floor(rs_ratio*r0);rp = r0 - rs
  GB_re = GB_Kmeans(Sigma_list, nk_list, rk_list,rs, p, tau = rs - 0.5, Ti = 50)
  P_GB = GB_re$Pksk_finetuned
  
  return(P_GB)
}
