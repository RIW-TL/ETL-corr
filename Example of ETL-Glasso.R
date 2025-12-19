library(mvtnorm)
library(ks)
library(Matrix)
library(glasso)
library(fastclime)
library(lavaSearch2) 
library(CVglasso)
library(clime)
#--------------------------------------------

source("Functions of Trans-CLIME.R")
source("Main-functions.R")
##################################################################

# settings (as an example)
#--------------------------------------------
n0 = 100;n1 = 500;p = 50;K = 5
num.info = 2;d = 10

# data generation
#-------------------------------------------
fit = dat.g(n0,n1,p,num.info,K,d,case = 2,s = 10)
Ome_t = fit$Sig_t
Ome_s = fit$Sig_s

Xt = rmvnorm(n0,rep(0,p), solve(Ome_t))
XS = array(0,c(n1,p,K))
for (k in 1:K) {XS[,,k] = rmvnorm(n1,rep(0,p),solve(Ome_s[,,k]))}

##################################################################

# Target-Glasso
#----------------------------------------------
#----------------------------------------------
sam.cov = cov(Xt)
lam = seq(0.05,0.12,0.02)
target_glasso_hat = CVglasso(Xt,S = sam.cov, K = 5 ,trace = "none",lam = lam)$Omega
# norm(target_glasso_hat - Ome_t,"F")^2/p


# Pool-Glasso
#----------------------------------------------
#----------------------------------------------
X.A = c()
for (k in 1:K) {X.A = rbind(X.A,XS[,,k])}
sum.X = rbind(Xt,X.A)
pool.cov = cov(sum.X)

lam = seq(0.01,0.05,0.01)
pool_glasso_hat = CVglasso(Xt,S = pool.cov, K = 5 ,trace = "none",lam = 0.03)$Omega
# norm(pool_glasso_hat - Ome_t,"F")^2/p


# Trans-CLIME 
#------------------------------------------
#------------------------------------------
n.vec = c(n0,rep(n1,K))
n00 = round(n.vec[1]*4/5) 
const = 0.5
Theta.re0 = Myfastclime.s(X = Xt[1:n00,], Bmat = diag(1,p), lambda = const*2*sqrt(log(p)/n00))
Theta.init = Theta.re0$Theta.hat
Omega.tl1 = Trans.CLIME(X = Xt[1:n00,], X.A, const = const, 
                        X.til = Xt[(n00 + 1):n0,], Theta.cl = Theta.init)
ind2 = (n0 - n00 + 1):n0
Theta.re0 = Myfastclime.s(X = Xt[ind2,], Bmat = diag(1,p), lambda = const*2*sqrt(log(p)/length(ind2)))
Theta.init = Theta.re0$Theta.hat
Omega.tl2 = Trans.CLIME(X = Xt[ind2,], X.A, const = const,
                        X.til = Xt[1:(n0 - n00),], Theta.cl = Theta.init)
trans_clime_hat = (Omega.tl1 + Omega.tl2)/2
# norm(trans_clime_hat - Ome_t,"F")^2/p


# ETL-Glasso
#-------------------------------
#-------------------------------

# Step 1: ETL-corr
#------------------
R_t_hat = cor(Xt)
R_s_hat = array(0,c(p,p,K))
for (k in 1:K) {R_s_hat[,,k] = cor(XS[,,k])}

pre = pre_compute(Xt, XS, R_t_hat, R_s_hat)
ETL_corr_hat = fit.ETL.corr(Xt, XS, A = 0.5, M = 0.4, R_t_hat, R_s_hat,pre)

# Step 2: ETL-var
#-------------------------------
ETL.var.hat = fit.ETL.var(Xt,XS,A1 = 0.5,M1 = 0.3,A2 = 0.3,M2 = 0.5)

# Step 3: ETL-cov 
#------------------
ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})

# Step 4: precision matrix estimation
#----------------------------------
lam = seq(0.002,0.01,0.002)
fit = CVglasso(Xt, S = ETL_cov_hat, K = 5, trace = "none", lam = lam)
ETL_glasso_hat = fit$Omega 
# norm(ETL_glasso_hat - Ome_t,"F")^2/p

##################################################################
error.glasso = matrix(0,1,4)
colnames(error.glasso) = c("ETL-Glasso","Target-Glasso","Pool-Glasso","Trans-CLIME")
error.glasso[1,] = c(norm(ETL_glasso_hat - Ome_t,"F")^2/p,
                     norm(target_glasso_hat - Ome_t,"F")^2/p,
                     norm(pool_glasso_hat - Ome_t,"F")^2/p,
                     norm(trans_clime_hat - Ome_t,"F")^2/p);error.glasso



