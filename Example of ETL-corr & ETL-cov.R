
source("Main-functions.R")
#--------------------------------------------

# settings (as an example)
#--------------------------------------------
n0 = 100;n1 = 500;p = 100;K = 5;case = 1
num.info = 2;d = 8

# data generation
#-------------------------------------------
data = dat.g(n0,n1,p,num.info,K,d,case,s = 10)
Xt = data$Xt 
XS = data$XS

Sig_t = data$Sig_t
# Sig_s = data$Sig_s

diag.elements = diag(Sig_t)
R_t = diag(diag.elements^{-1/2}) %*% Sig_t %*% diag(diag.elements^{-1/2})

R_s = array(0,c(p,p,K))
for (k in 1:K) {
  
  diag.elements = diag(Sig_s[,,k])
  R_s[,,k] = diag(diag.elements^{-1/2}) %*% Sig_s[,,k] %*% diag(diag.elements^{-1/2})
}

#==================================cor============================
#---------------------------------------------------------------

# Target-corr
#------------------
Target_corr_hat = cor(Xt)

# Target-Scor
#------------------
Target_Scor_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)

# Pool-corr
#-----------------
X.A = c()
for (k in 1:K) {X.A = rbind(X.A, XS[,,k])}
sum.X = rbind(Xt,X.A)
Pool_corr_hat = cor(sum.X)

# ETL-corr
#------------------
R_t_hat = Target_corr_hat
R_s_hat = array(0,c(p,p,K))
for (k in 1:K) {R_s_hat[,,k] = cor(XS[,,k])}

pre <- pre_compute(Xt, XS, R_t_hat, R_s_hat)
ETL_corr_hat = fit.ETL.corr(Xt, XS, A = 0.5, M = 0.4, R_t_hat, R_s_hat,pre)

#---------------------------------------------------------
# A.range = seq(0.3,0.9,0.3)
# M.range = c(0.05,seq(0.1,0.5,0.1))
# cv.loss.F = matrix(0, length(A.range), length(M.range))
# 
# for (ir in seq_along(A.range)) {
#   for (im in seq_along(M.range)) {
#     
#     A = A.range[ir]
#     M = M.range[im]
#     
#     l = (ir - 1)*length(M.range) + im
#     out = fit.ETL.corr(Xt, XS, K, A, M, R_t_hat, R_s_hat,pre)
#     cv.loss.F[ir,im] <- norm(out - R_t, "F")^2/p
#   }
# }
# cv.loss.F
#---------------------------------------------------------


# ETL-Scor 
#-----------------
ETL_Scor_hat = fit.ETL.Scor(R_trans = ETL_corr_hat,lambda = 0.1)


#-----------------------------------------------------
# lambda.range = seq(0.01,0.1,0.02)
# cv.loss.F = rep(0,length(lambda.range))
# 
# for (m in 1:length(lambda.range)) {
#   
#   lambda = lambda.range[m]
#   cv.SECM = fit.ETL.Scor(R_trans = ETL_corr_hat,lambda = 0.1)
#   
#   cv.loss.F[m] = norm(cv.SECM - R_t,"F")^2/p
# }
# cv.loss.F
#----------------------------------------------------


error.F = matrix(0,1,5)
colnames(error.F) = c("ETL-corr","ETL-Scor","Pool-corr",
                      "Target-corr","Target-Scor")
error.F[1,] = c(norm(ETL_corr_hat - R_t,"F")^2/p,
                norm(ETL_Scor_hat - R_t,"F")^2/p,
                norm(Pool_corr_hat - R_t,"F")^2/p,
                norm(Target_corr_hat - R_t,"F")^2/p,
                norm(Target_Scor_hat - R_t,"F")^2/p);error.F


#==========================var==============================
#===========================================================

# Taget-var
#------------------------------
target.var.hat = diag(cov(Xt))

# Pool
#------------------------------
pool.var.hat = diag(cov(sum.X))

# ETL-var
#-------------------------------
ETL.var.hat = fit.ETL.var(Xt,XS,A1 = 0.5,M1 = 0.3,A2 = 0.3,M2 = 0.5)

#-------------------------------------------------------
# A.range = seq(0.3,0.9,0.3)
# M.range = seq(0.2,1,0.1)
# cv.loss = matrix(0,length(A.range),length(M.range))
# for (a in 1:length(A.range)) {
#   for (m in 1:length(M.range)) {
#     
#     A2 = A.range[a]
#     M2 = M.range[m]
#     
#     m.trans.hat = fit.2rd.moment.trans(Xt,XS,K,A2,M2)
#     
#     # cv.loss[a,m] = sum((m.trans.hat)^2)
#     cv.loss[a,m] = sum((m.trans.hat - diag(Sig_t))^2)
#   }
# }
# cv.loss
#--------------------------------------------------------


#-------------------------------
error.var = matrix(0,1,3)
colnames(error.var) = c("Target-var","Pool-var","ETL-var")
error.var[1,] = c(sum((target.var.hat - diag(Sig_t))^2),
              sum((pool.var.hat - diag(Sig_t))^2),
              sum((ETL.var.hat - diag(Sig_t))^2));error.var

#==========================cov==============================
#===========================================================

# Target-cov
#------------------
target_cov_hat = cov(Xt)

# Target-Scov
#------------------
target_Scov_hat = est_sparseCov(Xt, method = 'cv', operator = 'hard', corr = T)


# ETL-cov 
#------------------
ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})

# ETL-Scov
#-----------------
ETL_Scov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_Scor_hat %*% diag(ETL.var.hat^{1/2})

# pool
#-----------------
pool_cov_hat = cov(sum.X)

#------------------------------------------

error.cov.F = matrix(0,1,5)
colnames(error.cov.F) = c("ETL-cov","ETL-Scov","Pool-cov","Target-Scov","Target-cov")
error.cov.F[1,] = c(norm(ETL_cov_hat - Sig_t,"F")^2/p,
                norm(ETL_Scov_hat - Sig_t,"F")^2/p,
                norm(pool_cov_hat - Sig_t,"F")^2/p,
                norm(target_Scov_hat - Sig_t,"F")^2/p,
                norm(target_cov_hat - Sig_t,"F")^2/p);error.cov.F

