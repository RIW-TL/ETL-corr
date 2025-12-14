
Yt = scale(Xt,center = TRUE,scale = TRUE)
YS = list();K = 1
YS[[1]] = scale(XS,center = TRUE,scale = TRUE)


iter = 50
PE.out = matrix(0,iter,6)
colnames(PE.out) = c("Target-Glasso","Trans-ECM-Glasso","Trans-SECM-Glasso",
                     "Pooled-Glasso","CLIME","Trans-CLIME")

for (ii in 1:iter) {
  
  test.index = sample(1:nrow(Yt),floor(nrow(Yt)*0.3))
  Yt.test = Yt[test.index,]
  Yt.train = Yt[-test.index,]
  
  ## Target-glasso
  ##------------------------------------------------------------
  #-------------------------------------------------------------
  sam.cov = cov(Yt.train)
  
  rho.list = seq(0.005,0.1,0.005)
  fit = glasso.cv(sam.cov,Yt.test,rho.list);fit$error
  # fit = CVglasso(Yt.train, sam.cov, K = 5, trace = "none", lam = rho.list)
  Ome_single_hat = fit$Ome_hat
  PE.single = Dist(Ome_single_hat,Yt.test)
  
  # Trans-ECM-Glasso 
  ##------------------------------------------------------------
  #-------------------------------------------------------------
  m2.t.hat = moment(Yt.train,order = 2);m2.t.hat[1]
  
  m2.s.hat = matrix(0,p,K)
  for (k in 1:K) {m2.s.hat[,k] = moment(YS[[k]],order = 2)}
  trans.var.hat = fit.2rd.moment.trans(Yt.train,YS,K, M = 1, A = 0.5, m2.t.hat, m2.s.hat)
  
  # correlation matrix estimation
  R_t_hat = est_sparseCov(Yt.train, method = 'cv', operator = 'hard', corr = T)
  K = 1
  R_s_hat = array(0,c(p,p,K))
  R_s_hat[,,1] = est_sparseCov(YS[[1]], method = 'cv', operator = 'scad', corr = T)
  
  pre = precompute_c_independent(Yt.train, YS, K = 1, R_t_hat, R_s_hat)
  
  A.range = seq(0.3,1.2,0.3)
  M.range = seq(0,4,1)
  rho.list = seq(0.005,0.1,0.005)
  
  O = matrix(0,length(A.range),length(M.range))
  R_ECM = Ome_ECM = list()
  
  for (ir in 1:length(A.range)) {
    for (im in 1:length(M.range)) {
      
      A = A.range[ir]
      M = M.range[im]
      
      # A = 0.9;M = 3
      R_Trans_ECM_hat = fit.trans.ECM.GGM(Yt.train, YS, K, M, A, R_t_hat, R_s_hat, pre)
      
      # min(eigen(R_Trans_ECM_hat)$values)
      
      R_Trans_ECM_hat = fit.def.sym(R_Trans_ECM_hat)
      cov_Trans_ECM_hat = diag(trans.var.hat^{1/2}) %*% R_Trans_ECM_hat %*% diag(trans.var.hat^{1/2})
      
      l = (ir - 1)*length(M.range) + im
      R_ECM[[l]] = R_Trans_ECM_hat
      fit = glasso.cv(cov_Trans_ECM_hat,Yt.test,rho.list);
      # fit$error
      Ome_ECM[[l]] = fit$Ome_hat
      # fit = CVglasso(Yt.train, R_Trans_ECM_hat, K = 5, trace = "none", lam = rho.list)
      Ome_ECM_hat = fit$Ome_hat
      O[ir,im] = Dist(Ome_ECM_hat,Yt.test)
    }
  }
  
  min.id = which(O == min(O),arr.ind = TRUE)
  l.opt = (min.id[1,1] - 1)*length(M.range) + min.id[1,2]
  R_Trans_ECM_hat = R_ECM[[l.opt]]
  Ome_ECM_hat = Ome_ECM[[l.opt]]
  PE.ECM = min(O)
  
  # Trans-SECM-Glasso 
  #------------------------------------------
  lambda.range = seq(0.01,0.1,0.01) 
  O.SECM = rep(0,length(lambda.range))
  Ome_SECM_hat = list()
  R_Trans_SECM_hat = R_Trans_ECM_hat
  
  for (m in 1:length(lambda.range)) {
    
    lambda = lambda.range[m]
    R_Trans_SECM_hat[abs(R_Trans_SECM_hat) < lambda] = 0
    R_Trans_SECM_hat = fit.def.sym(R_Trans_SECM_hat)
    cov_Trans_SECM_hat = diag(trans.var.hat^{1/2}) %*% R_Trans_SECM_hat %*% diag(trans.var.hat^{1/2})
    
    #--------------
    # fit = CVglasso(Xt.train, S = cov_Trans_SECM_hat, K = 5, trace = "none", lam = lam.ECM)
    # Ome_SECM_hat[[m]] = fit$Omega
    fit = glasso.cv(cov_Trans_SECM_hat,Yt.test,rho.list)
    Ome_SECM_hat[[m]] = fit$Ome_hat
    O.SECM[m] = Dist(Ome_SECM_hat[[m]],Yt.test)
  }
  
  min.id = which.min(O.SECM)
  Ome_SECM_hat = Ome_SECM_hat[[min.id]]
  PE.SECM = min(O.SECM)
  
  ## pool-glasso
  #-------------------------------------------------
  #-------------------------------------------------
  rho.list = seq(0.01,0.1,0.01)
  sum.X = Yt.train
  for (k in 1:K) {
    
    sum.X = rbind(sum.X,YS[[k]])
  }
  sam.pool = cov(sum.X)
  
  fit = glasso.cv(sam.pool,Yt.test,rho.list);fit$error
  Ome_pool_hat = fit$Ome_hat
  # fit = CVglasso(Yt.train, sam.pool, K = 5, trace = "none", lam = rho.list)
  # Ome_pool_hat = fit$Omega
  
  PE.pool = Dist(Ome_pool_hat,Yt.test)
  
  
  # CLIME
  #---------------------------------------
  const.list = seq(0.2,1,0.2)
  Theta.init = list();clime.error = c()
  
  for (s in 1:length(const.list)) {
    
    const = const.list[s]
    Theta.re0 = Myfastclime.s(X = Yt.train, Bmat = diag(1,p),
                              lambda = 2*const*sqrt(log(p)/nrow(Yt.train)))
    Theta.init.0 = Theta.re0$Theta.hat
    Theta.init[[s]] = (Theta.init.0 + t(Theta.init.0))/2
    clime.error[s] = Dist(Theta.init[[s]],Yt.test)
  }
  
  min.id = which.min(clime.error)
  Theta.initial = Theta.init[[which.min(clime.error)]]
  
  PE.clime = min(clime.error)
  
  ## Trans-CLIME
  ##------------------------------------------
  n0 = nrow(Yt.train)
  X.A = sum.X[-c(1:n0),]
  const.list = seq(0.2,1,0.2)
  error = c();Omega.trans.clime.hat = list()
  for (s in 1:length(const.list)) {
    
    const = const.list[s]
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
    Omega.trans.clime.hat[[s]] = (Omega.tl + t(Omega.tl))/2
    error[s] = Dist(Omega.trans.clime.hat[[s]],Yt.test)
  }
  
  min.id = which.min(error)
  Omega.trans.clime.hat = Omega.trans.clime.hat[[min.id]]
  
  PE.trans.clime = min(error)
  
  #--------------------------------------------
  PE.out[ii,] = c(PE.single,PE.ECM,PE.SECM,PE.pool,PE.clime,PE.trans.clime)
}


#-------------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)

out1 = as.data.frame(PE.out)[,-3]*2
colnames(out1) = c("Target-Glasso","ETL-Glasso","Pool-Glasso","CLIME","Trans-CLIME")
breaks.new = c("ETL-Glasso","Pool-Glasso","Target-Glasso",
               "Trans-CLIME", "CLIME")
out_long <- out1 %>%
  mutate(Trial = 1:n()) %>%
  pivot_longer(-Trial, names_to = "Method", values_to = "Value") %>%
  mutate(Method = factor(Method, levels = breaks.new))

ggplot(out_long, aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA) +
  # scale_colour_manual("",breaks = breaks.new)+
  ylim(0.15,0.45) +
  # theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +  # ????? x ??????????????????
  labs(title = "",
       x = "",
       y = "Prediction Error") + 
  theme_bw() +
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),'cm')) + 
  theme(legend.position = "top")
