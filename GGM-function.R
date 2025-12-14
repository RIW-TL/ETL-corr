library(mvtnorm)
library(ks)
library(Matrix)
library(glasso)
library(fastclime)
library(lavaSearch2) 
library(clime)
library(CVglasso)
library(clime)
#-------------------------------------

dat.ggm.g = function(n0,n1,p,num.info,K,d){
  
  ## parameter generation
  ##---------------------------------------
  
  # target covariance matrix
  B = matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in i:p) {
      
      if(i == j){
        
        B[i, j] = 1 + runif(1, 0, 0.2)
        
      }else {
        
        value = runif(1, 0, 0.2)
        B[i, j] <- value
        B[j, i] <- value
      }
    }
  }
  
  B = (B + t(B))/2
  
  # sparse covariance matrix
  s = 10
  
  for(j in 1:p){
    
    B[,j]<-B[,j]*(abs(B[,j])>=quantile(abs(B[,j]),1-s/p))
    B[j,]<-B[j,]*(abs(B[j,])>=quantile(abs(B[j,]),1-s/p))
  }  
  
  Sig_t = fit.def.sym(B) 
  
  
  # source covariance matrices
  Sig_s = array(0,c(p,p,K))
  
  fit = eigen(Sig_t)
  U = fit$vectors
  Lambda.T = fit$values
  
  for (k in 1:K) {
    
    if(k <= num.info){
      
      # Lambda.S = Lambda.T + runif(p,min = 0,max = d)
      Lambda.S = Lambda.T
      
      US1 = U + matrix(runif(p^2,0,d),p)
      qr_result = qr(US1)
      US1 = qr.Q(qr_result)
      # US1 = U
      
      Sig_s[,,k] = US1 %*% diag(Lambda.S) %*% t(US1)
    }else{
      
      Lambda.S = Lambda.T + runif(p,min = 1,max = 2)
      
      US2 = U + matrix(runif(p^2,0.5,1),p)
      qr_result = qr(US2)
      US2 = qr.Q(qr_result)
      
      Sig_s[,,k] = US2 %*% diag(Lambda.S) %*% t(US2)
    }
    
    Sig_s[,,k] = fit.def.sym(Sig = Sig_s[,,k])
  }
  
  Ome_t = Sig_t
  Ome_s = Sig_s
  
  ## data generation
  ##-----------------------------------------------
  Xt = rmvnorm(n0,rep(0,p), solve(Ome_t))
  # X.test = rmvnorm(n.test,rep(0,p), solve(Ome_t))
  
  XS = array(0,c(n1,p,K))
  for (k in 1:K) {XS[,,k] = rmvnorm(n1,rep(0,p),solve(Ome_s[,,k]))}
  
  return(list(Xt = Xt,XS = XS, 
              Ome_t = Ome_t, Ome_s = Ome_s))
}

## Trans-CLIME
Myfastclime.s = function(X, Bmat,lambda = 0.1, scale = T, n){
  
  p<-ncol(X)
  obj=rep(-1,2*p)
  obj_bar=rep(0,2*p)
  rhs_bar<-rep(1, 2*p)
  if(isSymmetric(X, tol=10^(-4))){
    Sig.hat<-X
  }else{
    Sig.hat<-cov(X)
  }
  Sig.hat0<-Sig.hat
  feasible=T
  Theta.hat<-NULL
  Sig.diag<-diag(Sig.hat)
  Sig.hat<-cov2cor(Sig.hat)
  mat=rbind(cbind(Sig.hat,-Sig.hat),
            cbind(-Sig.hat,Sig.hat))
  for(j in 1:p){
    rhs <- c(Bmat[,j],-Bmat[,j])
    out.txt<-capture.output(  fastlp.re<-fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
    if(!grepl("optimal", out.txt) ){
      feasible=F
      break
    }
    Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)])/sqrt(Sig.diag[j])/sqrt(Sig.diag))
    if(scale & sum(Theta.hat[,j]==0)==p){
      feasible=F
      break
    }else if(scale){
      Theta.hat[,j]<-  as.numeric(Bmat[j,j]/ (Sig.hat0[j,]%*%Theta.hat[,j]))*Theta.hat[,j]
      # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]
      
    }
  }
  if(!feasible){
    cat('Theta.hat not found','\n')
    Theta.hat<-solve(cov(X)+diag(lambda,p))%*%Bmat
  }
  list(Theta.hat=Theta.hat, conv=feasible)
}

Trans.CLIME = function(X,X.A, const, agg = T, X.til = NULL,Theta.cl){
  
  if(agg &is.null(X.til)){
    cat('no aggregation samples provided.','\n')
  }
  n0=nrow(X)
  nA<-nrow(X.A)
  p<-ncol(X)
  sigA.hat<-mean(apply(X.A, 2, sd))
  sig0.hat<-mean(apply(X, 2, sd))
  
  omega.l1<-mean(apply(Theta.cl,2, function(x) sum(abs(x))))
  Delta.re <- Myfastclime.s(X=diag(1,p), Bmat=diag(1,p)-t(Theta.cl)%*%cov(X.A), 
                            lambda = 2*omega.l1*sqrt(log(p)/n0), scale = F)
  if(Delta.re$conv){Delta.hat<-Delta.re$Theta.hat}else{ Delta.hat<-diag(0,p) }
  Theta.re <- Myfastclime.s(X=cov(X.A), Bmat=diag(1,p)-t(Delta.hat),
                            lambda = 2*const*sqrt(log(p)/nA))
  Theta.hat<-Theta.re$Theta.hat
  
  if(agg){
    
    Omega.hat <- Agg(Theta.init=cbind(Theta.cl, Theta.hat), X.til=X.til)
  }else{
    Omega.hat<-Theta.hat
  }
  Omega.hat
}

Agg = function(Theta.init, X.til){
  p<-ncol(X.til)
  n.til<-nrow(X.til)
  v.mat<-sapply(1:p, function(j){
    W.j<-cov(X.til%*%cbind(Theta.init[,j], Theta.init[,p+j]))
    v0=rep(0,2)
    v0[which.min(c(W.j[1,1]-2*Theta.init[j,j], W.j[2,2]-2*Theta.init[j,p+j]))]<-1
    v0
  })
  Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])
  
  Theta.hat
}

Spd.proj = function(SigA.hat, eps=NULL){
  p=ncol(SigA.hat)
  if(is.null(eps)){
    eps<-5/p
  }
  feasible=1
  SigA.t<-SigA.hat
  if(min(eigen(SigA.t)$values) <=eps ){
    feasible=2
    SigA.t<-ADMM_proj(SigA.hat, epsilon=eps)$mat
  }
  SigA.t<-lavaSearch2:::symmetrize(SigA.t, update.upper = TRUE)
  list(mat=SigA.t, conv=feasible)
}

l1proj = function(v, b){
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}

ADMM_proj = function(mat,
                     epsilon=1e-4,
                     mu=10,
                     it.max=1e3,
                     etol=1e-4,
                     etol_distance = 1e-4){
  
  
  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE) 
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L     
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)    
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))
  
}

Dist = function(Ome_hat,X.test){
  
  p = dim(Ome_hat)[1]
  eigens = eigen(Ome_hat)$values
  value = sum(diag(cov(X.test) %*% Ome_hat))/(2*p) - 
    sum(log(eigens[eigens > 0]))/(2*p)
  return(value)
}

fit.trans.ECM.GGM = function(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre) {
  
  p <- dim(Yt)[2]
  R_trans <- diag(1, p)
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      fenzi.sum <- fenmu.sum <- rep(0, K)
      
      for (k in 1:K) {
        
        # k = 5;j = 1;l = l + 1
        err <- pre$errors[[j, l, k]]
        weight <- pre$weights[, k, j, l]
        
        # boxplot(err$eta,main = "eta")
        # boxplot(err$error_1,main = "error1")
        # boxplot(weight,main = "weight")
        
        # A = 1;M = 0.02
        select.id1 <- which(abs(err$error_1) <= A)
        select.id2 <- which(abs(err$eta) <= M)
        id <- intersect(select.id1, select.id2)
        
        # boxplot(err$eta[id],main = "eta[id]")
        # boxplot(err$error_1[id],main = "error1[id]")
        # boxplot(weight[id],main = "weight[id]")
        
        fenzi.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, j, k])
        fenmu.sum[k] <- sum(YS[id, l, k] * weight[id] * YS[id, l, k])
      }
      
      fenzi <- sum(fenzi.sum) + pre$Yt_cross[j, l]
      fenmu <- sum(fenmu.sum) + pre$Yt_cross[l, l]
      R_trans[j, l] <- fenzi / fenmu
      R_trans[l, j] <- R_trans[j, l]
    }
  }
  
  return(R_trans)
}

fit.trans.ECM.GGM.cv = function(Xt, XS, K, Ome_c_t, R_t_hat, R_s_hat, 
                                A.range, M.range) {
  
  n0 <- dim(Xt)[1]
  p <- dim(Xt)[2]
  n1 <- dim(XS[,,1])[1]
  Yt <- scale(Xt)
  YS <- array(0, c(n1, p, K))
  for (k in 1:K) {
    YS[,,k] <- scale(XS[,,k])
  }
  
  
  pre <- precompute_c_independent(Yt, YS, K, R_t_hat, R_s_hat)
  
  # M.range = seq(0.02,0.1,0.02)
  cv.loss.F = rep(0, length(A.range), length(M.range))
  R_ecm_hat = R_secm_hat = list()
  
  for (a in seq_along(A.range)) {
    for (m in seq_along(M.range)) {
      
      M = M.range[m]
      A = A.range[a]
      
      # correlation matrix estimation
      #----------------------------------
      l = (a - 1)*length(M.range) + m
      R_ecm_hat[[l]] = fit.trans.ECM.GGM(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre)
      
      # precision matrix estimation
      #----------------------------------
      lam = seq(0.005,0.2,0.005)
      fit = CVglasso(Yt, S = R_ecm_hat[[l]], K = 5, trace = "none", lam = lam)
      Ome_ECM_c_hat = fit$Omega
      
      cv.loss.F[a,m] = norm(Ome_ECM_c_hat - Ome_c_t,"F")^2/p
    }
  }
  
  min.id = which.min(cv.loss.F)
  
  return(list(R_trans_ECM_hat = R_ecm_hat[[min.id]],
              F.loss = min(cv.loss)))
}

fit.def.sym = function(Sig){
  
  p = dim(Sig)[2]
  eigen_decomp = eigen(Sig)
  eigenvalues = eigen_decomp$values
  
  if(min(eigenvalues) >= 0){
    
    Sig = Sig
  }else{
    
    c = 1e-8 - min(eigenvalues)
    Sig = Sig + c*diag(1,p)
  }
  
  Sig[abs(Sig)< 10^{-10}] = 0
  
  Sig = (Sig + t(Sig))/2
  return(Sig)
}

precompute_c_independent = function(Yt, YS, K, R_t_hat, R_s_hat) {
  
  p <- dim(Yt)[2]
  n1 <- dim(YS)[1]
  
  pre <- list()
  
  pre$Yt_cross <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in 1:p) {
      pre$Yt_cross[j, l] <- sum(Yt[, j] * Yt[, l])
    }
  }
  
  
  pre$errors <- array(list(), dim = c(p, p, K))  
  pre$weights <- array(0, dim = c(n1, K, p, p))  
  # pre$A <- array(0, dim = c(p, p, K))  
  # pre$h_jl <- array(0, dim = c(p, p, K))  
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      for (k in 1:K) {
        
        eta <- (R_t_hat[j, l] - R_s_hat[j, l, k]) * YS[, l, k]
        error_0 <- YS[, j, k] - R_t_hat[j, l] * YS[, l, k]
        error_1 <- YS[, j, k] - R_s_hat[j, l, k] * YS[, l, k]
        pre$errors[[j, l, k]] <- list(eta = eta, error_0 = error_0, error_1 = error_1)
        
        kde.fit <- kde(x = error_1)
        fenzi <- (predict(kde.fit, x = error_0) + predict(kde.fit, x = -error_0)) / 2
        # fenzi <- predict(kde.fit, x = error_0) 
        fenmu <- predict(kde.fit, x = error_1)
        pre$weights[, k, j, l] <- fenzi / fenmu 
        
        # rho <- 0.7
        # pre$A[j, l, k] <- quantile(abs(error_1), rho)
        # pre$h_jl[j, l, k] <- abs(R_t_hat[j, l] - R_s_hat[j, l, k])
      }
    }
  }
  
  return(pre)
}

glasso.cv = function(Sig,rho.range,Ome_t){
  
  cv.loss = rep(0,length(rho.range))
  Ome_hat = list()
  
  for (q in 1:length(rho.range)) {
    
    p = dim(Sig)[2]
    rho = rho.range[q] 
    Ome_hat[[q]] = glasso(Sig,rho)$wi
    cv.loss[q] = norm(Ome_hat[[q]] - Ome_t,"F")^2/p
  }
  
  min.id = which.min(cv.loss)
  
  return(list(Omega = Ome_hat[[min.id]], cv.loss = cv.loss))
}
