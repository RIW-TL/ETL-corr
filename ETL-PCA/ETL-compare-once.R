
#========================================================
# Data generation
#========================================================
k0 = 3;d = 6;ii = ii + 1
dat <- generate_sim4_data(k0 = k0, d = d, seed = ii)

X0 <- dat$X0
XS <- dat$XS
p = ncol(dat$X0)
K = length(XS)
r0 = dat$r0
U_t = dat$U_t
rk_list = dat$rk_list

Sigma_true = dat$Sig_t


#========================================================
# 1. ETL-corr
#========================================================

fit_corr <- fit.ETL.corr(
  X0 = X0,
  XS = XS,
  c_A1 = 0.5,
  c_M1 = 0.05,
  c_M2 = 0.1
)

ETL_corr_hat <- fit_corr$R_hat
H_min <- fit_corr$H_min


#========================================================
# 2. ETL-Scor
#========================================================

ETL_scor_hat <- fit.ETL.Scor(
  R_hat = ETL_corr_hat,
  H_min = H_min,
  c_T = 1
)


#========================================================
# 3. ETL-var
#========================================================

# fit <- ETL_var(
#   X0 = X0,
#   XS = XS,
#   A_grid = c(3,4,5,6),
#   M_grid = c(0.3,0.5)
# )
# 
# ETL_var_hat = fit$ETL_var_hat
ETL_var_hat = diag(Sigma_true)

#========================================================
# 4. ETL-cov
#========================================================

D_ETL <- diag(
  sqrt(ETL_var_hat)
)

ETL_cov_hat <-
  D_ETL %*%
  ETL_corr_hat %*%
  D_ETL

ETL_cov_hat <- make_pd(
  ETL_cov_hat,
  min_eig = 0.1
)


#========================================================
# 5. ETL-Scov
#========================================================

ETL_scov_hat <-
  D_ETL %*%
  ETL_scor_hat %*%
  D_ETL

ETL_scov_hat <- make_pd(
  ETL_scov_hat,
  min_eig = 0.1
)


#========================================================
# 6. Target-cov
#========================================================

Target_cov_hat <- cov(X0)


#========================================================
# 7. Pool-cov
#========================================================

X_all <- do.call(
  rbind,
  c(list(X0), XS)
)

Pool_cov_hat <- cov(X_all)


#========================================================
# 8. Target-Scov
#========================================================

Target_scov_hat <- fit.target.scor(X0, corr = FALSE)

Target_scov_hat <- as.matrix(
  Target_scov_hat
)

#========================================================
# GB-PCA
#========================================================

GB_hat_1 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.25)
GB_hat_2 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.5)
GB_hat_3 <- fit_GB(X0,XS,r0 = r0, rk_list, rs_ratio = 0.75)

#========================================================
# Compute losses
#========================================================

loss_ETL_cov <- cov_loss_fun(
  Sigma_hat = ETL_cov_hat,
  U_t = U_t,
  r = r0
)

loss_ETL_scov <- cov_loss_fun(
  Sigma_hat = ETL_scov_hat,
  U_t = U_t,
  r = r0
)

loss_Target_cov <- cov_loss_fun(
  Sigma_hat = Target_cov_hat,
  U_t = U_t,
  r = r0
)

loss_Target_scov <- cov_loss_fun(
  Sigma_hat = Target_scov_hat,
  U_t = U_t,
  r = r0
)

loss_Pool_cov <- cov_loss_fun(
  Sigma_hat = Pool_cov_hat,
  U_t = U_t,
  r = r0
)

loss_GB1_cov <- norm(GB_hat_1  - U_t %*% t(U_t),"F")
loss_GB2_cov <- norm(GB_hat_2  - U_t %*% t(U_t),"F")
loss_GB3_cov <- norm(GB_hat_3  - U_t %*% t(U_t),"F")
#========================================================
# Store covariance losses
#========================================================

result.cov <- c(
  loss_ETL_cov,
  loss_ETL_scov,
  loss_Target_cov,
  loss_Target_scov,
  loss_Pool_cov,
  loss_GB1_cov,
  loss_GB2_cov,
  loss_GB3_cov
)
names(result.cov) = c(
  "ETL-cov",
  "ETL-Scov",
  "Target-cov",
  "Target-Scov",
  "Pool-cov",
  "GB-1",
  "GB-2",
  "GB-3"
)
result.cov
fit$A_hat
fit$M_hat