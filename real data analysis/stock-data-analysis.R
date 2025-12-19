library(mvtnorm)
library(ks)
library(Matrix)
library(glasso)
library(fastclime)
library(lavaSearch2) 
library(clime)
library(sparseCov)
library(CVglasso)
library(moments)

source("Functions of Trans-CLIME.R")
source("Main-functions.R")
source("stock-data-process.R")

Dist = function(Ome_hat,X.test){
  
  p = dim(Ome_hat)[1]
  eigens = eigen(Ome_hat)$values
  value = sum(diag(cov(X.test) %*% Ome_hat))/p - 
    sum(log(eigens[eigens > 0]))/p
  return(value)
}
##################################################################

Yt = scale(Xt,center = TRUE,scale = TRUE)
YS = list();K = 1
YS[[1]] = scale(XS,center = TRUE,scale = TRUE)


iter = 100
PE.out = matrix(0,iter,5)
colnames(PE.out) = c("Target-Glasso","ETL-Glasso","Pool-Glasso","CLIME","Trans-CLIME")

for (ii in 1:iter) {
  
  test.index = sample(1:nrow(Yt),floor(nrow(Yt)*0.3))
  Yt.test = Yt[test.index,]
  Yt.train = Yt[-test.index,]
  
  ## Target-Glasso
  ##------------------------------------------------------------
  #-------------------------------------------------------------
  sam.cov = cov(Yt.train)
  
  rho.list = seq(0.005,0.1,0.005)
  fit = CVglasso(Yt.train, sam.cov, K = 5, trace = "none", lam = rho.list)
  Ome_target_hat = fit$Ome_hat
  PE.target.Glasso = Dist(Ome_target_hat,Yt.test)
  
  # ETL-Glasso 
  ##------------------------------------------------------------
  #-------------------------------------------------------------
  # Step 1: ETL-corr
  #------------------
  R_t_hat = cor(Yt.train)
  R_s_hat = array(0,c(p,p,K))
  for (k in 1:K) {R_s_hat[,,k] = cor(YS[[k]])}
  
  pre = pre_compute(Yt.train, YS, R_t_hat, R_s_hat)
  ETL_corr_hat = fit.ETL.corr(Yt.train, YS, A = 0.5, M = 0.4, R_t_hat, R_s_hat,pre)
  
  # Step 2: ETL-var
  #-------------------------------
  ETL.var.hat = fit.ETL.var(Yt.train,YS,A1 = 0.5,M1 = 0.3,A2 = 0.3,M2 = 0.5)
  
  # Step 3: ETL-cov 
  #------------------
  ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})
  
  # Step 4: precision matrix estimation
  #----------------------------------
  lam = seq(0.002,0.01,0.002)
  fit = CVglasso(Xt, S = ETL_cov_hat, K = 5, trace = "none", lam = lam)
  ETL_glasso_hat = fit$Omega 
  PE.ETL.Glasso = Dist(ETL_glasso_hat,Yt.test)
  
  
  ## Pool-Glasso
  #-------------------------------------------------
  #-------------------------------------------------
  sum.X = Yt.train
  for (k in 1:K) {
    
    sum.X = rbind(sum.X,YS[[k]])
  }
  sam.pool = cov(sum.X)
  
  rho.list = seq(0.01,0.1,0.01)
  fit = CVglasso(Yt.train, sam.pool, K = 5, trace = "none", lam = rho.list)
  Pool_Glasso_hat = fit$Omega
  PE.pool.Glasso = Dist(Pool_Glasso_hat,Yt.test)
  
  
  # CLIME
  #------------------------------------------
  #------------------------------------------
  const = 0.5
  Theta.re0 = Myfastclime.s(X = Yt.train, Bmat = diag(1,p),
                            lambda = 2*const*sqrt(log(p)/nrow(Yt.train)))
  Theta.init.0 = Theta.re0$Theta.hat
  Theta.initial = (Theta.init.0 + t(Theta.init.0))/2
  PE.clime = Dist(Theta.initial,Yt.test)  
  
  # Trans-CLIME
  #------------------------------------------
  #------------------------------------------
  n0 = nrow(Yt.train)
  X.A = sum.X[-c(1:n0),]
  const = 0.5
  n = round(n0*4/5)
  Omega.tl1 = Trans.CLIME(X = Yt.train[1:n,], X.A = X.A, 
                          const = const, agg = T,
                          X.til = Yt.train[(n + 1):n0,], 
                          Theta.cl = Theta.initial)
  ind  = (n0 - n + 1): n0
  Omega.tl2 = Trans.CLIME(X = Yt.train[ind,], X.A = X.A, 
                          const = const, agg = T,
                          X.til = Yt.train[-ind,], 
                          Theta.cl = Theta.initial)
  
  Omega.tl = (Omega.tl1 + Omega.tl2)/2
  Omega.trans.clime.hat = (Omega.tl + t(Omega.tl))/2
  PE.trans.clime = Dist(Omega.trans.clime.hat,Yt.test)  

  
  #--------------------------------------------
  PE.out[ii,] = c(PE.target.Glasso,PE.ETL.Glasso,PE.pool.Glasso,PE.clime,PE.trans.clime)
}


#-------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)

out1 = as.data.frame(PE.out)
colnames(out1) = c("Target-Glasso","ETL-Glasso","Pool-Glasso","CLIME","Trans-CLIME")

out_long <- out1 %>%
  mutate(Trial = 1:n()) %>%
  pivot_longer(-Trial, names_to = "Method", values_to = "Value") %>%
  mutate(Method = factor(Method, levels = breaks.new))

ggplot(out_long, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0.15,0.45) +
  labs(title = "",
       x = "",
       y = "Prediction Error") + 
  theme_bw() +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),'cm')) + 
  theme(legend.position = "top")
