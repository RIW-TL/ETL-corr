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

## Target-Glasso
##------------------------------------------------------------
#-------------------------------------------------------------
sam.cov = cov(Yt)

rho.list = seq(0.005,0.1,0.005)
fit = CVglasso(Yt, sam.cov, K = 5, trace = "none", lam = rho.list)
Ome_target_hat = fit$Ome_hat

# ETL-Glasso 
##------------------------------------------------------------
#-------------------------------------------------------------
# Step 1: ETL-corr
#------------------
R_t_hat = cor(Yt)
R_s_hat = array(0,c(p,p,K))
for (k in 1:K) {R_s_hat[,,k] = cor(YS[[k]])}

pre = pre_compute(Yt, YS, R_t_hat, R_s_hat)
ETL_corr_hat = fit.ETL.corr(Yt, YS, A = 0.5, M = 0.4, R_t_hat, R_s_hat,pre)

# Step 2: ETL-var
#-------------------------------
ETL.var.hat = fit.ETL.var(Yt,YS,A1 = 0.5,M1 = 0.3,A2 = 0.3,M2 = 0.5)

# Step 3: ETL-cov 
#------------------
ETL_cov_hat = diag(ETL.var.hat^{1/2}) %*% ETL_corr_hat %*% diag(ETL.var.hat^{1/2})

# Step 4: precision matrix estimation
#----------------------------------
lam = seq(0.002,0.01,0.002)
fit = CVglasso(Xt, S = ETL_cov_hat, K = 5, trace = "none", lam = lam)
ETL_glasso_hat = fit$Omega 


## Pool-Glasso
#-------------------------------------------------
#-------------------------------------------------
sum.X = Yt
for (k in 1:K) {
  
  sum.X = rbind(sum.X,YS[[k]])
}
sam.pool = cov(sum.X)

rho.list = seq(0.01,0.1,0.01)
fit = CVglasso(Yt, sam.pool, K = 5, trace = "none", lam = rho.list)
Pool_Glasso_hat = fit$Omega


# CLIME
#------------------------------------------
#------------------------------------------
const = 0.5
Theta.re0 = Myfastclime.s(X = Yt, Bmat = diag(1,p),
                          lambda = 2*const*sqrt(log(p)/nrow(Yt)))
Theta.init.0 = Theta.re0$Theta.hat
Theta.initial = (Theta.init.0 + t(Theta.init.0))/2

# Trans-CLIME
#------------------------------------------
#------------------------------------------
n0 = nrow(Yt)
X.A = sum.X[-c(1:n0),]
const = 0.5
n = round(n0*4/5)
Omega.tl1 = Trans.CLIME(X = Yt[1:n,], X.A = X.A, 
                        const = const, agg = T,
                        X.til = Yt[(n + 1):n0,], 
                        Theta.cl = Theta.initial)
ind  = (n0 - n + 1): n0
Omega.tl2 = Trans.CLIME(X = Yt[ind,], X.A = X.A, 
                        const = const, agg = T,
                        X.til = Yt[-ind,], 
                        Theta.cl = Theta.initial)

Omega.tl = (Omega.tl1 + Omega.tl2)/2
Omega.trans.clime.hat = (Omega.tl + t(Omega.tl))/2


#===============================================================
#===============================================================
edge = function(Ome_hat){
  
  p = ncol(Ome_hat)
  
  density.in = c()
  density.in[1] = length(which(Ome_hat[1:35,1:35] != 0))/(35^2)
  density.in[2] = length(which(Ome_hat[36:67,36:67] != 0))/(32^2)
  density.in[3] = length(which(Ome_hat[68:104,68:104] != 0))/(37^2)

  
  density.cross = c()
  density.cross[1] = length(which(Ome_hat[1:35,36:67] != 0))/(35*32)
  density.cross[2] = length(which(Ome_hat[1:35,68:104] != 0))/(35*37)
  density.cross[3] = length(which(Ome_hat[36:67,68:104] != 0))/(32*37)
  
  return(list(density.in = density.in, density.cross = density.cross))
}


edge.detect.result = matrix(0,5,6)
rownames(edge.detect.result) = c("Target-Glasso","ETL-Glasso","Pool-Glasso","CLIME","Trans-CLIME")
colnames(edge.detect.result) = c("class-1","class-2","class-3","class1-2","class1-3","class2-3")

a = edge(Ome_target_hat)
edge.detect.result[1,] = c(a$density.in, a$density.cross)

a = edge(ETL_glasso_hat)
edge.detect.result[2,] = c(a$density.in, a$density.cross)

a = edge(Pool_Glasso_hat)
edge.detect.result[3,] = c(a$density.in, a$density.cross)

a = edge(Theta.initial)
edge.detect.result[4,] = c(a$density.in, a$density.cross)

a = edge(Omega.trans.clime.hat)
edge.detect.result[5,] = c(a$density.in, a$density.cross)
