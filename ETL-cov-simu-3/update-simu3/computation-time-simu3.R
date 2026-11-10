
source("data-generate-3.R")
source("ETL-functions-3.R")
source("ETL-var.R")

time.fun = function(p,K,iter = 10,ncores = 8){
  
  time.out = matrix(0,iter,2)
  
  for (ii in 1:iter) {
    
    #==============================================
    # data generation
    #==============================================
    dat <- generate_sim3_data(eta = 0.5, K = K, p = p, seed = ii)
    
    X0 <- dat$X0
    XS <- dat$XS
    
    ##================Implementation=========================
    ##=======================================================
    
    start.time.cor = Sys.time()
    fit <- fit.ETL.corr.cv(X0, XS, M_grid = c(0.2, 0.4, 0.6,0.8), cv_folds = 5, ncores = ncores)
    ETL_corr_hat = fit$R_hat
    end.time.cor = Sys.time()
    time.cor = as.numeric(end.time.cor - start.time.cor)
    
    
    start.time.var = Sys.time()
    fit <- ETL_var_cv(
      X0 = X0,
      XS = XS,
      A_grid = c(2.5,3,3.5,4),
      M_grid = c(0.2,0.4,0.6,0.8),
      ncores = ncores
    )
    ETL_var_hat = fit$ETL_var_hat
    end.time.var = Sys.time()
    time.var = as.numeric(end.time.var - start.time.var)
    
    
    time.out[ii,] = c(time.cor,time.var)
  }
  
  return(colMeans(time.out))
}

a1 = time.fun(p = 50,K = 1);a1
a2 = time.fun(p = 50,K = 3);a2
a3 = time.fun(p = 50,K = 5);a3

b1 = time.fun(p = 100,K = 1);b1
b2 = time.fun(p = 100,K = 3);b2
b3 = time.fun(p = 100,K = 5);b3







