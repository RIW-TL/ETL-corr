library(sparseCov)
library(ks)
library(moments)


##================ETL-corr========================
##================================================

############################################################
## Stage 1: Pre-compute KDE and importance weights only once
############################################################

precompute.ETL.corr.parallel = function(
    X0,
    XS,
    R_t_initial,
    R_s_initial,
    c_h2 = 1,
    density_floor = 1e-12,
    ncores = parallel::detectCores() - 2
){
  
  ############################################################
  ## Basic dimensions
  ############################################################
  
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  if (length(dim(XS)) < 3) {
    
    dim(XS) = c(dim(XS), 1)
    dim(R_s_initial) = c(dim(R_s_initial), 1)
  }
  
  K = dim(XS)[3]
  n1 = dim(XS)[1]
  
  
  ############################################################
  ## Data normalization
  ############################################################
  
  Y0 = scale(X0)
  
  YS = array(
    0,
    c(n1, p, K)
  )
  
  for (k in 1:K) {
    YS[,,k] = scale(XS[,,k])
  }
  
  
  ############################################################
  ## Initial correlation matrices
  ############################################################
  
  # R_t_initial = R_t
  # R_s_initial = R_s
  # R_t_initial = cor(X0)
  # 
  # R_s_initial = array(
  #   0,
  #   c(p, p, K)
  # )
  # 
  # for (k in 1:K) {
  #   R_s_initial[,,k] = cor(XS[,,k])
  # }
  
  
  ############################################################
  ## All pairs (j,l), j < l
  ############################################################
  
  pair.id = which(
    upper.tri(matrix(0, p, p)),
    arr.ind = TRUE
  )
  
  
  ############################################################
  ## Target cross-products
  ############################################################
  
  target.cross = crossprod(Y0)
  
  
  ############################################################
  ## KDE + importance weights
  ############################################################
  
  compute.weight = function(ii) {
      
      j = pair.id[ii,1]
      l = pair.id[ii,2]
      
      
      ########################################################
      ## Store weights for all sources
      ########################################################
      
      weight.mat = matrix(
        0,
        nrow = n1,
        ncol = K
      )
      
      
      for (k in 1:K) {
        
        
        ######################################################
        ## Extract variables
        ######################################################
        
        yl = as.numeric(
          YS[,l,k]
        )
        
        yj = as.numeric(
          YS[,j,k]
        )
        
        
        ######################################################
        ## Residuals
        ######################################################
        
        error_0 =
          yj -
          R_t_initial[j,l] * yl
        
        
        error_1 =
          yj -
          R_s_initial[j,l,k] * yl
        
        
        ######################################################
        ## Training data for joint density
        ##
        ## f(error_1, Y_l)
        ######################################################
        
        d2.fenmu = cbind(
          error_1,
          yl
        )
        
        d2.fenmu = as.matrix(
          d2.fenmu
        )
        
        storage.mode(d2.fenmu) = "double"
        
        
        ######################################################
        ## Numerator evaluation points
        ######################################################
        
        d2.fenzi.1 = cbind(
          error_0,
          yl
        )
        
        d2.fenzi.2 = cbind(
          -error_0,
          yl
        )
        
        
        d2.fenzi.1 = as.matrix(
          d2.fenzi.1
        )
        
        d2.fenzi.2 = as.matrix(
          d2.fenzi.2
        )
        
        
        storage.mode(d2.fenzi.1) = "double"
        storage.mode(d2.fenzi.2) = "double"
        
        
        ######################################################
        ## Evaluate density only at needed points
        ######################################################
        
        eval.points = rbind(
          d2.fenmu,
          d2.fenzi.1,
          d2.fenzi.2
        )
        
        eval.points = as.matrix(
          eval.points
        )
        
        storage.mode(eval.points) = "double"
        
        
        ############################################################
        ## Normal-scale/reference bandwidth
        ## Exact simplification under standardized data
        ############################################################
        
        rho.s =
          R_s_initial[j, l, k]
        
        H2 =
          n1^(-1/3) *
          matrix(
            c(
              1 - rho.s^2, 0,
              0,             1
            ),
            nrow = 2,
            ncol = 2,
            byrow = TRUE
          )
        
        ############################################################
        ## Bandwidth multiplier
        ############################################################
        
        H2 =
          c_h2^2 *
          H2
        
        ############################################################
        ## Small ridge for numerical stability
        ############################################################
        
        H2 =
          H2 +
          1e-8 * diag(2)
        
        
        ######################################################
        ## Joint KDE
        ##
        ## No 151 x 151 grid
        ## No contour
        ######################################################
        
        # kde.fit.2d = ks::kde(
        #   
        #   x = d2.fenmu,
        #   
        #   H = H2,
        #   
        #   eval.points = eval.points,
        #   
        #   binned = FALSE,
        #   
        #   compute.cont = FALSE
        # )
        
        kde.fit.2d = ks::kde(

          x = d2.fenmu,

          H = H2,

          eval.points = eval.points,

          binned = TRUE,

          bgridsize = c(81, 81),

          compute.cont = FALSE
        )
        
        
        ######################################################
        ## Extract density values
        ######################################################
        
        density.all =
          as.numeric(
            kde.fit.2d$estimate
          )
        
        
        ######################################################
        ## Denominator
        ##
        ## f(error_1, Y_l)
        ######################################################
        
        fenmu =
          density.all[
            seq_len(n1)
          ]
        
        
        ######################################################
        ## Numerator part 1
        ##
        ## f(error_0, Y_l)
        ######################################################
        
        fenzi.1 =
          density.all[
            n1 + seq_len(n1)
          ]
        
        
        ######################################################
        ## Numerator part 2
        ##
        ## f(-error_0, Y_l)
        ######################################################
        
        fenzi.2 =
          density.all[
            2 * n1 + seq_len(n1)
          ]
        
        
        ######################################################
        ## Symmetrized numerator density
        ######################################################
        
        fenzi =
          (
            fenzi.1 +
              fenzi.2
          ) / 2
        
        
        ######################################################
        ## Avoid zero denominator
        ######################################################
        
        fenmu =
          pmax(
            fenmu,
            density_floor
          )
        
        
        ######################################################
        ## Importance weights
        ##
        ## Conditional density ratio:
        ##
        ## f(error_0 | Y) / f(error_1 | Y)
        ##
        ## = [f(error_0,Y)/f(Y)] /
        ##   [f(error_1,Y)/f(Y)]
        ##
        ## = f(error_0,Y) / f(error_1,Y)
        ######################################################
        
        weight.mat[,k] =
          fenzi /
          fenmu
      }
      
      
    return(weight.mat)
  }

  if (ncores <= 1L) {
    weight.list = lapply(seq_len(nrow(pair.id)), compute.weight)
  } else {
    cl = parallel::makeCluster(ncores)

    on.exit(
      parallel::stopCluster(cl),
      add = TRUE
    )

    parallel::clusterExport(
      cl,
      c(
        "YS",
        "R_t_initial",
        "R_s_initial",
        "pair.id",
        "K",
        "n1",
        "c_h2",
        "density_floor",
        "compute.weight"
      ),
      envir = environment()
    )

    weight.list = parallel::parLapply(
      cl,
      seq_len(nrow(pair.id)),
      compute.weight
    )
  }
  
  
  ############################################################
  ## Return precomputed quantities
  ############################################################
  
  out = list(
    
    Y0 = Y0,
    
    YS = YS,
    
    R_t_initial = R_t_initial,
    
    R_s_initial = R_s_initial,
    
    target.cross = target.cross,
    
    weight.list = weight.list,
    
    pair.id = pair.id,
    
    n0 = n0,
    
    n1 = n1,
    
    p = p,
    
    K = K,
    
    c_h2 = c_h2
  )
  
  
  return(out)
}

############################################################
## Stage 2: Fast ETL-corr for given tuning parameters
############################################################

fit.ETL.corr.fast = function(
    precomputed,
    c_A1 = 0.6,
    c_A2 = 1,
    c_M1 = 0.5,
    c_M2 = 1
){
  
  Y0 = precomputed$Y0
  YS = precomputed$YS
  
  R_t_initial = precomputed$R_t_initial
  R_s_initial = precomputed$R_s_initial
  
  target.cross = precomputed$target.cross
  
  weight.list = precomputed$weight.list
  pair.id = precomputed$pair.id
  
  p = precomputed$p
  K = precomputed$K
  
  
  ##########################################################
  ## Initialize
  ##########################################################
  
  R_trans = S_num_matrix = diag(p)
  
  
  ##########################################################
  ## Loop over pairs
  ##########################################################
  
  for (ii in seq_len(nrow(pair.id))) {
    
    j = pair.id[ii,1]
    l = pair.id[ii,2]
    
    fenzi.sum = rep(0, K)
    fenmu.sum = rep(0, K)
    
    
    ########################################################
    ## Sources
    ########################################################
    
    for (k in 1:K) {
      
      yl = YS[,l,k]
      yj = YS[,j,k]
      
      
      ######################################################
      ## Quantities needed for sample selection
      ######################################################
      
      error_1 =
        yj -
        R_s_initial[j,l,k] * yl
      
      
      eta =
        (
          R_t_initial[j,l] -
            R_s_initial[j,l,k]
        ) * yl
      
      
      ######################################################
      ## Sample selection
      ######################################################
      
      A1 = quantile(
        abs(error_1),
        c_A1
      )
      
      A2 = quantile(
        abs(yl),
        c_A2
      )
      
      
      M =
        c_M1 /
        (
          c_M2 +
            abs(
              R_t_initial[j,l] -
                R_s_initial[j,l,k]
            )
        )
      
      
      id = which(
        abs(error_1) <= A1 &
          abs(yl) <= A2 &
          abs(eta) <= M
      )
      
      S_num_matrix[j,l] = length(id)
      S_num_matrix[l,j] = S_num_matrix[j,l]
      
      ######################################################
      ## Retrieve pre-computed weights
      ######################################################
      
      weight = weight.list[[ii]][,k]
      
      
      ######################################################
      ## Weighted estimator
      ######################################################
      
      fenzi.sum[k] = sum(
        weight[id] *
          yl[id] *
          yj[id]
      )
      
      
      fenmu.sum[k] = sum(
        weight[id] *
          yl[id]^2
      )
    }
    
    
    ########################################################
    ## Plus target data
    ########################################################
    
    R_trans[j,l] =
      (
        sum(fenzi.sum) +
          target.cross[l,j]
      ) /
      (
        sum(fenmu.sum) +
          target.cross[l,l]
      )
    
    
    R_trans[l,j] =
      R_trans[j,l]
  }
  
  
  S_num_min = min(abs(S_num_matrix[upper.tri(S_num_matrix)]))
  
  return(list(R_trans, S_num_min))
}


##================ETL-scor========================
##================================================
fit.ETL.Scor.cv = function(R_trans, X0, R_true, 
                           lambda.range = seq(0.02,0.2,0.02),
                           rule = c("oracle","stein","theory"),
                           c_0, S_num_min){
  
  p = dim(X0)[2]
  
  if(rule == "theory"){
    
    thre.opt = c_0*sqrt(log(p)/S_num_min)
    
    R_strans = R_trans
    
    for (j in 1:(p - 1)) {
      for (l in (j + 1):p) {
        
        if(abs(R_strans[j,l]) < thre.opt){R_strans[j,l] = R_strans[l,j] = 0}
      }
    }
    
    loss = 0
    
  }else{
    
    cv.SECM = list()
    loss = rep(0,length(lambda.range))
    
    for (m3 in 1:length(lambda.range)) {
      
      R_strans = R_trans
      lambda = lambda.range[m3]
      
      for (j in 1:(p - 1)) {
        for (l in (j + 1):p) {
          
          if(abs(R_strans[j,l]) < lambda){R_strans[j,l] = R_strans[l,j] = 0}
        }
      }
      
      cv.SECM[[m3]] = R_strans
      
      if(rule == "stein"){
        
        loss[m3] = stein.loss(cv.SECM[[m3]],X0)
      }else{
        
        loss[m3] = norm(cv.SECM[[m3]] - R_true,"F")^2/p
      }
      
    }
    
    id = which.min(loss)
    R_strans = cv.SECM[[id]]
    thre.opt = lambda.range[id]
  }
  
  list(R_scor_hat = R_strans, loss = loss, thre.opt = thre.opt)
}

stein.loss <- function(A, X0, center = TRUE, normalize = FALSE) {
  
  A <- (A + t(A)) / 2
  
  n0 <- nrow(X0)
  p  <- ncol(X0)
  
  if (nrow(A) != p || ncol(A) != p) {
    stop("Dimension of A does not match X0.")
  }
  
  # Center target data if necessary
  if (center) {
    X <- scale(X0, center = TRUE, scale = FALSE)
  } else {
    X <- X0
  }
  
  # Cholesky decomposition: A = R^T R
  R <- chol(A)
  
  # log det(A)
  logdet.A <- 2 * sum(log(diag(R)))
  
  # x_i^T A^{-1} x_i
  # Avoid explicitly computing solve(A)
  Z <- t(backsolve(R, t(X), transpose = TRUE))
  quad.term <- mean(rowSums(Z^2))
  
  # Stein / Gaussian quasi-likelihood loss
  loss <- logdet.A + quad.term
  
  # Optional per-dimension normalization
  if (normalize) {
    loss <- loss / p
  }
  
  return(loss)
}

##================ETL-var========================
##================================================
fit.mean.trans = function(X0, XS, c_m1_A, M1){
  
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  if (length(dim(XS)) < 3) {dim(XS) = c(dim(XS), 1)}
  
  K = dim(XS)[3]
  n1 = dim(XS)[1]
  
  #-----------------------------------------------
  # mean.t.hat = colMeans(X0)
  
  mean.s.hat = matrix(0,p,K)
  for (k in 1:K) {
    
    mean.s.hat[,k] = colMeans(XS[,,k])
  }
  
  mean.t.hat = rep(0,p)

  mean.s.hat = matrix(0,p,K)
  for (k in 1:K) {

    mean.s.hat[,k] = rep(0,p)
  }
  
  mean.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = mean.t.hat[j] - mean.s.hat[j,]
    kp.id = which(abs(eta) <= M1)
    
    if(length(kp.id) == 0){
      
      mean.trans[j] = colMeans(X0)[j]
    }else{
      
      if(length(kp.id) == 1){
        
        XS.use = array(0,c(n1,p,2))
        XS.use[,,kp.id] = XS[,,kp.id]
      }else{
        
        XS.use = XS[,,kp.id]
      }
      
      fenzi.sum = fenmu.sum = rep(0,length(kp.id))
      
      for (k in 1:length(kp.id)) {
        
        # the estimation of density functions
        error_0 = XS.use[,j,k] - mean.t.hat[j]
        error_1 = XS.use[,j,k] - mean.s.hat[j,k]
        # kde.fit = kde(x = error_1, h = c_m1_h * hns(error_1))
        kde.fit = kde(x = error_1)
        
        # the estimation of importance weights
        fenzi = (predict(kde.fit,x = error_0) + predict(kde.fit,x = -error_0))/2
        fenmu = predict(kde.fit,x = error_1)
        weight = fenzi/fenmu
        
        # sample selection
        m1_A = quantile(abs(error_1),c_m1_A)
        id = which(abs(error_1) <= m1_A)
        
        # weights estimator
        fenzi.sum[k] = sum(XS.use[id,j,k]*weight[id])
        fenmu.sum[k] = sum(weight[id])
      }
      
      mean.trans[j] = (sum(X0[,j]) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
      
      
    }
  }
  
  return(mean.trans)
}

fit.m2.trans = function(X0,XS,c_m2_A, M2){
  
  n0 = dim(X0)[1]
  p = dim(X0)[2]
  
  if (length(dim(XS)) < 3) {dim(XS) = c(dim(XS), 1)}
  
  K = dim(XS)[3]
  n1 = dim(XS)[1]
  
  #-----------------------------------------------
  # m2.t.hat = moment(X0,order = 2)
  
  m2.s.hat = matrix(0,p,K)
  for (k in 1:K) {
    
    m2.s.hat[,k] = moment(XS[,,k], order = 2)
  }
  #-----------------------------------------------
  
  m2.trans = rep(0,p)
  
  for (j in 1:p) {
    
    eta = m2.t.hat[j] - m2.s.hat[j,]
    kp.id = which(abs(eta) <= M2)
    
    if(length(kp.id) == 0){
      
      m2.trans[j] = m2.t.hat[j]
    }else{
      
      if(length(kp.id) == 1){
        
        XS.use = array(0,c(n1,p,2))
        XS.use[,,kp.id] = XS[,,kp.id]
      }else{
        
        XS.use = XS[,,kp.id]
      }
      
      fenzi.sum = fenmu.sum = rep(0,length(kp.id))
      
      for (k in 1:length(kp.id)) {
        
        # the estimation of density functions
        error_0 = (XS.use[,j,k])^2 - m2.t.hat[j]
        error_1 = (XS.use[,j,k])^2 - m2.s.hat[j, kp.id[k]]
        
        # c_m2_h = 0.5
        # kde.fit = kde(x = error_1, h = c_m2_h*hns(error_1))
        kde.fit = kde(x = error_1)
        
        # the estimation of importance weights
        fenzi = (predict(kde.fit, x = error_0) + predict(kde.fit, x = -error_0))/2
        fenmu = predict(kde.fit,x = error_1)
        weight = fenzi/fenmu
        
        # sample selection
        # c_m2_A = 0.8
        m2_A = quantile(abs(error_1),c_m2_A)
        id = which(abs(error_1) <= m2_A)
        
        # boxplot(error_1[id])
        # boxplot(weight[id])
        
        # weights estimator (mean estimation)
        fenzi.sum[k] = sum(((XS.use[id,j,k])^2)*weight[id])
        fenmu.sum[k] = sum(weight[id])
      }
      
      m2.trans[j] = (sum((X0[,j])^2) + sum(fenzi.sum))/(n0 + sum(fenmu.sum))
    }
    
  }
  
  return(m2.trans)
}
