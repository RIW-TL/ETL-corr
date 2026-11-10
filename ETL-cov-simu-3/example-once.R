n0 = 100;
n1 = 500
p = 100
eta = 0.1
C_t = 0.1
#----------------------------------------------
    
##==============================
## data generation
##==============================
fit = dat.g.band(n0, n1, p, eta, C_t, seed = 100 + ii) 
X0 = fit$target_data
XS = fit$source_data

Sigma0 = fit$Sigma0
diag.elements = diag(Sigma0)
R_t = diag(diag.elements^(-1/2)) %*% Sigma0 %*% diag(diag.elements^(-1/2))

Sigma1 = fit$Sigma1
diag.elements = diag(Sigma1)
R_s = diag(diag.elements^(-1/2)) %*% Sigma1 %*% diag(diag.elements^(-1/2))

##=============================
## ADT
##=============================
fit.ADT = ADT(target_data = X0, source_data = XS)
sigma0.ADT.hat = fit.ADT$estimate

# norm(sigma0.ADT.hat - Sigma0,"F")^2/p

##============================
## Target-band
##============================
target_band_hat = fit_band(X0)$estimate


##=============================
## ETL-cov
##=============================

ETL.start.time = Sys.time()

R_t_initial = fit_band(X0)$estimate
R_s_initial = fit_band(XS)$estimate

# Step 1: ETL-corr estimator
#-----------------------------------------
precomputed = precompute.ETL.corr.parallel(
  X0,
  XS,
  R_t_initial,
  R_s_initial,
  c_h2 = 1,
  density_floor = 1e-12,
  ncores = parallel::detectCores() - 2)

fit.ETL.corr = fit.ETL.corr.fast(precomputed,
                                   c_A1 = 0.6, # c_A1 \in \{0.4,0.6,0.8\}
                                   c_A2 = 1,   # c_A2 \in \{0.6,0.8,1\}
                                   c_M1 = 0.5,
                                   c_M2 = 1)
trans_corr_hat = fit.ETL.corr[[1]]

ETL.corr.end.time = Sys.time()

S_num_min = fit.ETL.corr[[2]]

## loss of ETL-corr
# loss.ETL.corr = norm(fit.ETL.corr$trans_corr_hat - R_t,"F")^2/p;loss.ETL.corr

##=============================
## ETL-Scov
##=============================

fit.scor = fit.ETL.Scor.cv(trans_corr_hat, X0, R_true = R_t, 
                           lambda.range = seq(0.02,0.3,0.04),
                           rule = "theory", c_0 = 1, S_num_min)
trans_scor_hat = fit.scor$R_scor_hat

ETL.scor.end.time = Sys.time()


# # loss of ETL-corr
# loss.ETL.scor = norm(trans_scor_hat - R_t,"F")^2/p;loss.ETL.scor


# Step 2: ETL-var estimator
#-----------------------------------------
ETL.var.start.time = Sys.time()

trans.m1.hat = fit.mean.trans(X0, XS, 
                              c_m1_A = 0.8, # c_m1_A \in \{0.4,0.6,0.8\}
                              c_m1_h = 1, 
                              M = 0.3)  

trans.m2.hat = fit.m2.trans(X0, XS, 
                            c_m2_A = 0.6, # c_m1_A \in \{0.4,0.6,0.8\}
                            c_m2_h = 0.5, 
                            M = 0.6)
ETL.var.hat = trans.m2.hat - trans.m1.hat^2

## loss of ETL-var
# loss.ETL.var = sum(abs(ETL.var.hat - diag(Sigma0)))/p;loss.ETL.var

ETL.var.end.time = Sys.time()

## Step 3: correlation-variance combination
##-----------------------------------------
ETL.cov.hat = diag(ETL.var.hat^{1/2}) %*% trans_corr_hat %*% diag(ETL.var.hat^{1/2})
ETL.scov.hat = diag(ETL.var.hat^{1/2}) %*% trans_scor_hat %*% diag(ETL.var.hat^{1/2})

##=============================
## Tar-cov
##=============================
sigma0.target.hat = cov(X0)

##=============================
## Pool-cov
##=============================
sigma0.pool.hat = cov(rbind(X0,XS))

##=============================
## Target-Scov
##=============================

target_Scov_hat = est_sparseCov(X0, method = 'cv', operator = 'hard', corr = FALSE)


##=============================
## loss
##=============================

# loss of all methods on covariance matrix estimation
out.loss = c(
  norm(ETL.cov.hat - Sigma0,"F")^2/p,
  norm(ETL.scov.hat - Sigma0,"F")^2/p,
  norm(sigma0.ADT.hat - Sigma0,"F")^2/p,
  norm(sigma0.target.hat - Sigma0,"F")^2/p,
  norm(sigma0.pool.hat - Sigma0,"F")^2/p,
  norm(target_Scov_hat - Sigma0,"F")^2/p,
  norm(target_band_hat - Sigma0,"F")^2/p)
names(out.loss) <- c("ETL-cov","ETL-scov", "ADT", "Target-cov",
                     "Pool-cov","Target-scov","Target-band")
out.loss


