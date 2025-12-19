library(mvtnorm)
library(ks)
library(Matrix)
#--------------------------------------------
source("Functions of GB-PCA.R")
source("Main-functions.R")
##################################################################

# settings (as an example)
#--------------------------------------------
n0 = 100;n1 = 500;p = 50;K = 5
num.info = 2;d = 10

# data generation
#-----------------------------------------
data = dat.g.PCA(n0,n1,p,num.info,K,d)

Xt = data$Xt
XS = data$XS
rk_list = data$rk_list
r0 = rk_list[[1]]

Sig_t = data$Sig_t
Sig_s = data$Sig_s
fit  = eigen(Sig_t)
U_t = fit$vectors[,1:r0]

R_t = diag(diag(Sig_t)^{-1/2}) %*% Sig_t %*% diag(diag(Sig_t)^{-1/2})
R_s = array(0,c(p,p,K))
for (k in 1:K) {R_s[,,k] = diag(diag(Sig_s[,,k])^{-1/2}) %*% Sig_s[,,k] %*% diag(diag(Sig_s[,,k])^{-1/2})}

###################################################

# Target-PCA 
#------------------------------------------
#------------------------------------------
sam.cov = cov(Xt)
Target_PCA_hat = eigen(sam.cov)$vectors[,1:r0]
error_Target_PCA = fit.evaluate(Target_PCA_hat,U_t);error_Target_PCA


# Pool-PCA 
#------------------------------------------
#------------------------------------------ 
X.A = c()
for (k in 1:K) {X.A = rbind(X.A,XS[,,k])}
sum.X = rbind(Xt,X.A)
pool.cov = cov(sum.X)

Pool_PCA_hat = eigen(pool.cov)$vectors[,1:r0]
error_Pool_PCA = fit.evaluate(Pool_PCA_hat,U_t)

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
error_GB = fit.evaluate(U_GB,U_t);error_GB


# ETL-PCA
#-------------------------------
#-------------------------------

# Step 1: ETL-corr
#------------------
R_t_hat = cor(Xt)
R_s_hat = array(0,c(p,p,K))
for (k in 1:K) {R_s_hat[,,k] = cor(XS[,,k])}

pre = pre_compute(Xt, XS, R_t_hat, R_s_hat)
ETL_corr_hat = fit.ETL.corr(Xt, XS, A = 0.5, M = 0.2, R_t_hat, R_s_hat,pre)

#----------------------------------------------------------
# A.range = seq(0.5,2,0.5)
# M.range = seq(0.1,0.5,0.1)
# cv.loss = matrix(0, length(A.range), length(M.range))
# for (ir in seq_along(A.range)) {
#   for (im in seq_along(M.range)) {
#     
#     A = A.range[ir]
#     M = M.range[im]
#     
#     ETL_corr_hat = fit.ETL.corr(Xt, XS, A, M, R_t_hat, R_s_hat,pre)
#     # ETL.var.hat = fit.ETL.var(Xt,XS,A1 = 0.5,M1 = 0.5,A2 = 0.5,M2 = 0.5)
#     ETL.var.hat = diag(Sig_t)
#     ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})
#     
#     # Step 4: PCA
#     #----------------------------------
#     ETL_PCA_hat = eigen(ETL_cov_hat)$vectors[,1:r0]
#     cv.loss[ir,im] = fit.evaluate(ETL_PCA_hat,U_t)
#   }
# }
# cv.loss
#-------------------------------------------------------

# Step 2: ETL-var
#-------------------------------
ETL.var.hat = fit.ETL.var(Xt,XS,A1 = 0.5,M1 = 0.3,A2 = 0.3,M2 = 0.5)

# Step 3: ETL-cov 
#------------------
ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})

# Step 4: PCA
#----------------------------------
ETL_PCA_hat = eigen(ETL_cov_hat)$vectors[,1:r0]
error_ETL_PCA = fit.evaluate(ETL_PCA_hat,U_t);error_ETL_PCA


##################################################################
error.PCA = matrix(0,1,4)
colnames(error.PCA) = c("ETL-PCA","Target-PCA","Pool-PCA","GB")
error.PCA[1,] = c(error_ETL_PCA,error_Target_PCA,error_Pool_PCA,error_GB);error.PCA


