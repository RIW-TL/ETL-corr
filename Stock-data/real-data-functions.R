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
#--------------------------------------------------


precompute_c_independent = function(Yt, YS, K, R_t_hat, R_s_hat) {
  
  p <- dim(Yt)[2]
  
  pre <- list()
  pre$Yt_cross <- matrix(0, p, p)
  for (j in 1:p) {
    for (l in 1:p) {
      pre$Yt_cross[j, l] <- sum(Yt[, j] * Yt[, l])
    }
  }
  
  
  pre$errors <- array(list(), dim = c(p, p, K))  
  pre$weights <- vector("list", K)
  
  for (k in 1:K) {
    
    pre$weights[[k]] = matrix(vector("list", p*p), nrow = p, ncol = p)
  }
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      for (k in 1:K) {
        
        n_k <- nrow(YS[[k]])
        
        eta <- (R_t_hat[j, l] - R_s_hat[j, l, k]) * YS[[k]][, l]
        error_0 <- YS[[k]][, j] - R_t_hat[j, l] * YS[[k]][, l]
        error_1 <- YS[[k]][, j] - R_s_hat[j, l, k] * YS[[k]][, l]
        pre$errors[[j, l, k]] <- list(eta = eta, error_0 = error_0, error_1 = error_1)
        pre$errors[[l, j, k]] <- list(eta = eta, error_0 = error_0, error_1 = error_1)
        
        kde.fit <- kde(x = error_1)
        fenzi <- (predict(kde.fit, x = error_0) + predict(kde.fit, x = -error_0)) / 2
        fenmu <- predict(kde.fit, x = error_1)
        
        pre$weights[[k]][[j, l]] <- fenzi / fenmu 
        pre$weights[[k]][[l, j]] <- fenzi / fenmu 
      }
    }
  }
  
  return(pre)
}

fit.trans.ECM.GGM = function(Yt, YS, K, M, A, R_t_hat, R_s_hat, pre) {
  
  p <- dim(Yt)[2]
  R_trans <- diag(1, p)
  
  for (j in 1:(p - 1)) {
    for (l in (j + 1):p) {
      
      fenzi.sum <- fenmu.sum <- SU <- rep(0, K)
      
      for (k in 1:K) {
        
        err <- pre$errors[[j, l, k]]
        weight <- pre$weights[[k]][[j, l]]
        
        
        # boxplot(err$eta,main = "eta")
        # boxplot(err$error_1,main = "error1")
        # boxplot(weight,main = "weight")
        
        # A = 0.2;M = 0
        select.id1 <- which(abs(err$error_1) <= A)
        select.id2 <- which(abs(err$eta) <= M)
        id <- intersect(select.id1, select.id2)
        
        # boxplot(err$eta[id],main = "eta[id]")
        # boxplot(err$error_1[id],main = "error1[id]")
        # boxplot(weight[id],main = "weight[id]")
        
        fenzi.sum[k] <- sum(YS[[k]][id, l] * weight[id] * YS[[k]][id, j])
        fenmu.sum[k] <- sum(YS[[k]][id, l] * weight[id] * YS[[k]][id, l])
        # SU[k] = length(id)
      }
      
      fenzi <- sum(fenzi.sum) + pre$Yt_cross[j, l]
      fenmu <- sum(fenmu.sum) + pre$Yt_cross[l, l]
      R_trans[j, l] <- fenzi / fenmu;
      
      R_trans[l, j] <- R_trans[j, l]
    }
  }
  
  return(R_trans)
}

Dist = function(Ome_hat,X.test){
  
  p = dim(Ome_hat)[1]
  eigens = eigen(Ome_hat)$values
  value = sum(diag(cov(X.test) %*% Ome_hat))/(2*p) - 
    sum(log(eigens[eigens > 0]))/(2*p)
  return(value)
}

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

detect.edge = function(Ome_hat,Ome_single_hat){
  
  Ome_hat[which(abs(Ome_hat) < 1e-6)] = 0
  Ome_single_hat[which(abs(Ome_single_hat) < 1e-6)] = 0
  
  p = ncol(Ome_hat)
  Ome_hat_nondiag = Ome_hat - diag(diag(Ome_hat))
  Ome_single_hat = Ome_single_hat - diag(diag(Ome_single_hat))
  S = length(which(Ome_hat_nondiag != 0))/(p*(p - 1))
  PR = length(which(Ome_hat_nondiag != 0 & Ome_single_hat != 0))/length(which(Ome_single_hat != 0))
  NR = length(which(Ome_hat_nondiag != 0 & Ome_single_hat == 0))/length(which(Ome_single_hat == 0))
  
  return(c(S,PR,NR))
}

degree.top.fun = function(Ome_hat,tau = 0.2){
  
  Ome_hat[which(abs(Ome_hat) < 1e-6)] = 0
  p = ncol(Ome_hat)
  
  degree = c()
  for (j in 1:p) {degree[j] = length(which(Ome_hat[,j] != 0))}
  top.rank = as.integer(rank(degree)[1:round(tau*p)])
  
  return(top.rank)
}

fit.2rd.moment.trans = function(Xt,XS,K, M, A, m2.t.hat, m2.s.hat){
  
  n0 = nrow(Xt)
  p = ncol(Xt)
  
  # transfer-enhanced (2rd) moment estimation
  m2.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = m2.t.hat[j] - m2.s.hat[j,]
    kp.id = which(abs(eta) <= M)
    
    if(length(kp.id) == 0){
      
      m2.trans[j] = m2.t.hat[j]
      
    }else{
      
      XS.use = list()
      
      if(length(kp.id) == 1){
        
        XS.use[[1]] = XS[[kp.id]]
        
      }else{
        
        XS.use = XS[kp.id]
      }
      
      fenzi.sum = fenmu.sum = rep(0,length(kp.id))
      
      for (k in 1:length(kp.id)) {
        
        # the estimation of density functions
        error_0 = (XS.use[[k]][,j])^2 - m2.t.hat[j]
        error_1 = (XS.use[[k]][,j])^2 - m2.s.hat[j,k]
        kde.fit = kde(x = error_1)
        
        # the estimation of importance weights
        fenzi = (predict(kde.fit,x = error_0) +
                   predict(kde.fit,x = -error_0))/2
        fenmu = predict(kde.fit,x = error_1)
        weight = fenzi/fenmu
        
        # sample selection
        # A = quantile(abs(error_1),rho)
        id = which(abs(error_1) <= A)
        
        # weights estimator (mean estimation)
        fenzi.sum[k] = sum(((XS.use[[k]][id,j])^2)*weight[id])
        fenmu.sum[k] = sum(weight[id])
      }
      
      m2.trans[j] = (sum((Xt[,j])^2) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
    }
  }
  
  return(m2.trans)
}

fit.def.sym = function(Sig){
  
  p = dim(Sig)[2]
  eigen_decomp = eigen(Sig)
  eigenvalues = eigen_decomp$values
  
  if(min(eigenvalues) > 0){
    
    Sig = Sig
  }else{
    
    epsilon <- 1e-8
    c <- epsilon - min(eigenvalues)
    Sig = Sig + c*diag(1,p)
  }
  
  Sig[abs(Sig)< 10^{-10}] = 0
  
  Sig = (Sig + t(Sig))/2
  return(Sig)
}

glasso.cv = function(cov_hat,Xt.test,rho.list){
  
  glasso.error = c()
  Ome_hat = list()
  for (iv in seq_along(rho.list)) {
    
    rho = rho.list[iv]
    
    Ome_hat_0 = glasso(cov_hat,rho)$wi
    Ome_hat[[iv]] = (Ome_hat_0 + t(Ome_hat_0))/2
    glasso.error[iv] = Dist(Ome_hat[[iv]],Xt.test)
  }
  
  min.id = which.min(glasso.error)
  return(list(Ome_hat = Ome_hat[[min.id]],error = glasso.error))
}
