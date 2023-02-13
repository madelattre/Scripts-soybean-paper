## Creating block diagonal matrices (cf help package Matrix, Copyright (C) 2016 
## Martin Maechler, ETH Zurich)

bdiag_m <- function(lmat) {
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], 
            all(vapply(lmat, dim, integer(2)) == k)) 
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  new("dgCMatrix", Dim = c(M,M),
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

## SAEM function for parameter estimation

SAEM_GRM <- function(niter, nburnin, data, predictors, paraminit,
                     GRM_mat = GRM_mat, Xi, Zi, model, Nchain=2,Nsim=1) {
  ## niter     : total number of iterations of the algorithm
  ## nburnin   : number of iterations for the burnin phase of the algorithm
  ##             (nburnin < niter)
  ## data      : name of the dataset, in the tibble format
  ##             It should at least contain the following variables
  ##             - "id": the ID of the curves (ie the ID for the plants)
  ##             - "varID": the ID of the varieties 
  ##             - "y": the response variable
  ##             - "time": the measurement times
  ## paraminit : list of initial values for the parameters to be estimated
  ##             It should contain the following elements
  ##             - "beta": vector (fixed-effects)
  ##             - "P": matrix (covariance matrix of the random effects)
  ##             - "G": matrix (genetic covariance matrix)
  ##             - "sigma2": numeric (residual variance)
  ## GRM_mat   : genetic relationship matrix
  ## predictors: numbers of the columns in data containing the regression 
  ##             variables.  that enter in the definition of the curve model 
  ##             (eg time)
  ## Xi, Zi    : arrays containing the covariate matrices for each plant in the 
  ##             fixed-effects part and in the random effects part of the model
  ##             respectively
  ## *Unlike the variables in Xi and Zi, we are not interested in the effects of
  ## the variables informed through predictors. The variables in predictors are 
  ## structural variables that enter in the definition of the curve model, for 
  ## example time.* 
  ## model     : curve model function  
  ## Nchain    : number of Monte-Carlo chains used to simulate the plant parameters
  ## Nsim      : number of iterations for the Metropolis-Hastings algorithm used
  ##             to simulate the plant parameters
  
  ## Data processing and definition of important quantities

  d       <- dim(Zi)[3]/dim(GRM_mat)[1]
  nb.phi  <- dim(Xi)[2]
  nb.beta <- length(paraminit$beta)
  id      <- matrix(data$id)
  l_id    <- unique(id)
  N       <- length(l_id) # number of plants
  varID   <- as.matrix(data$varID)
  l_var   <- select(data, varID) %>% distinct()
  y       <- data$y  
  time    <- data$time
  A       <- GRM_mat
  invA    <- solve(A)
  Ntot    <- nrow(data) # total number of observations
  Na      <- length(unique(varID)) # number of varieties
   
  ##  Store individual variables in lists
  
  yi      <- list()
  timi    <- list()
  xij     <- list()
  Xilist  <- list()
  Zilist  <- list()
  tZilist <- list()
  
  for (i in 1:N){
    yi[[i]]      <- y[id == l_id[i]]
    timi[[i]]    <- time[id == l_id[i]]
    xij[[i]]     <- as.matrix(data[id == l_id[i],predictors])
    Xilist[[i]]  <- Xi[i,,]
    Zilist[[i]]  <- Zi[i,,]
    tZilist[[i]] <- t(Zi[i,,])
  }
  
  XX  <- do.call('rbind',Xilist)
  tZZ <- do.call('cbind',tZilist)
  ZZ  <- do.call('rbind',Zilist)
  

  
  ## Initialization of the parameters and storage of the parameter estimates
  
  beta   <- matrix(NA, niter + 1, nb.beta)
  G      <- array(NA, dim = c(niter + 1, d, d))
  P      <- array(NA, dim = c(niter + 1, nb.phi, nb.phi)) 
  sigma2 <- rep(NA, niter + 1)
  
  beta[1,]  <- as.vector(paraminit$beta)
  G[1, ,]   <- paraminit$G 
  P[1, ,]   <- paraminit$P 
  sigma2[1] <- paraminit$sigma2
  
  
  ## Initialization of the sufficient statistics used to update the parameter
  ## estimations
  
  s1 <- array(0, dim = c(niter + 1,nb.phi,N))
  s2 <- array(0, dim = c(niter + 1, d, d))
  s3 <- array(0, dim = c(niter + 1, nb.phi, nb.phi))
  s4 <- rep(0, niter + 1)
  
 
  ## Definition of the step sizes for the SA step of the algorithm
  ## The step size equals 1 during the burnin phase of the algorithm and then
  ## decreases
  
  gamma <- rep(1, niter)
  for (k in (nburnin + 1):niter) {
    gamma[k] <- (k - nburnin) ^ (-2 / 3)
  }
  
  ## Initialization of u and plant parameters
  
  z_phi <- array(0, dim = c(niter + 1, nb.phi, N, Nchain))
  for (i in 1:N){
    z_phi[1, ,i ,] <- rep(Xi[i,,]%*%beta[1,],Nchain)
  }
  z_u <- array(0, dim = c(niter + 1, Na * d, Nchain))
  
  ## Initialization of intermediate quantities
  #mco         <- matrix(0, niter, N)
  mco_MC      <- array(0, dim = c(niter, N, Nchain))
  mco_P_MC    <- array(0, dim = c(niter, N, Nchain, nb.phi,nb.phi))
  mco_beta_MC <- array(0, c(niter, N, Nchain, nb.phi))
  uinvAu_MC   <- array(0, c(niter, d, d, Nchain)) 
  
  
  for (k in 1:niter) {
    ## Loop of iterations of the algorithm
    G_k        <-  G[k, , ]
    invG_k     <- solve(G_k)
    P_k        <- P[k, , ]
    invP_k     <- solve(P_k)
    S_phi      <- P_k 
    invPkbloc  <- bdiag_m(replicate(N,invP_k,simplify=FALSE))
    ZP         <- tZZ%*%invPkbloc
    sumZPZ     <- ZP%*%ZZ
    invSigma_u <- sumZPZ + (invA %x% invG_k)
    Sigma_u    <- solve(invSigma_u)
    Xbetar     <- matrix(as.vector(XX%*%as.matrix(beta[k,],ncol=1)),nrow=nb.phi*N)
    
    
    ## S-step
    
    for (mc in 1:Nchain) {
      ## Metropolis within Gibbs sampler
      
      current_u <- matrix(z_u[k, , mc])
      
      ## 1- simulate plant \varphi_i's given u
      
      for (i in 1:N) {
        m_phi       <- Xilist[[i]] %*% beta[k, ] + Zilist[[i]] %*% current_u
        current_phi <- z_phi[k, , i, mc]
        for (j in 1:Nsim) {
          phi      <- mvnfast::rmvn(1, mu = m_phi , sigma = S_phi) 
          mu_num   <- model(phi,xij[[i]]) 
          mu_denom <- model(current_phi,xij[[i]]) 
          lograte  <- sum(log(dnorm(yi[[i]], mu_num, sqrt(sigma2[k])))) - 
            sum(log(dnorm(yi[[i]], mu_denom, sqrt(sigma2[k]))))
          v <- runif(1)
          cond1 <- (log(v) <= min(0, lograte))
          if (cond1) {
            current_phi <- phi
          }
          z_phi[k + 1, , i, mc] <- current_phi
        } 
      }
      
      ## 2- simulate u given simulated phi's
      
      dphi      <- matrix(as.vector(z_phi[k + 1, , , mc])) - Xbetar
      sumZPdphi <- ZP%*%dphi
      m_u       <- Sigma_u %*% sumZPdphi
      
      # simulate u
      z_u[k + 1, , mc] <- t(mvnfast::rmvn(1, mu = m_u, sigma = Sigma_u))
      mat_u            <- matrix(z_u[k + 1, , mc], d, Na) 
      
      uinvAu_MC[k, , , mc] <- mat_u %*% invA %*% t(mat_u)
      
      
      ## 3- compute necessary quantities to update the sufficient statistics
      ## approximations
      
      for (i in 1:N) {
        mco_P_MC[k, i, mc, , ]  <- (z_phi[k + 1, , i, mc] - Zilist[[i]] %*% z_u[k + 1, , mc])%*%t(z_phi[k + 1, , i, mc] - Zilist[[i]] %*% z_u[k + 1, , mc]) #- Xi[i,,] %*% beta[k, ]
        mco_beta_MC[k, i, mc, ] <- (z_phi[k + 1, , i, mc] - Zilist[[i]] %*% z_u[k + 1, , mc])
        mco_MC[k, i, mc]        <- sum((yi[[i]] - model(z_phi[k + 1,, i, mc], xij[[i]]))^2)
      }
      
    } 
    
    mco    <- apply(mco_MC[k, , ], 1, mean)
    uinvAu <- apply(uinvAu_MC[k, , , ], c(1, 2), mean)
    
    
    ## SA-step
  
    for (i in 1:N){
      s1[k + 1,,i] <- (1 - gamma[k]) * s1[k,,i] + gamma[k] * apply(mco_beta_MC[k, i, , ],2,mean)
    }
    s2[k + 1, , ] <- (1 - gamma[k]) * s2[k, , ] + gamma[k] * uinvAu
    s3[k + 1, , ] <- (1 - gamma[k]) * s3[k, , ] + gamma[k] * apply(mco_P_MC[k, , , , ],c(3,4),sum)/Nchain 
    s4[k + 1] <- (1 - gamma[k]) * s4[k] + gamma[k] * sum(mco)
    
    ## M-step
    
    G[k + 1, , ] <- s2[k + 1, , ] / Na
    
    XinvPX <- t(XX)%*%invPkbloc%*%XX
    
    s1list <- list()
    for (i in 1:N){
      s1list[[i]] <- t(t(s1[k + 1,,i]))
    }
    
    s1mat      <- do.call('rbind',s1list)
    XinvPs1    <- t(XX)%*%invPkbloc%*%s1mat
    beta[k+1,] <- as.vector(solve(XinvPX)%*%XinvPs1)
    
    s1col   <- do.call('cbind',s1list)
    Xbeta   <- matrix(as.vector(XX%*%as.matrix(beta[k+1,],ncol=1)),nrow=nb.phi) 
    Xbetas1 <- Xbeta%*%t(s1col)
    
    P[k+1,,] <- (s3[k+1,,] - Xbetas1 - t(Xbetas1) + Xbeta%*%t(Xbeta))/N 
    
    sigma2[k + 1] <- s4[k + 1] / Ntot
    
  } ## End of the algorithm iterations loop
  
  res <- list(beta = beta,
              G = G,
              P = P,
              sigma2 = sigma2,
              z_u = z_u,
              z_phi = z_phi)
  
  return(res)
}

## SAEM function for genetic effects prediction

predict.u <- function(niter, nburnin, data, predictors, param, GRM_mat, Xi, Zi, 
                      model, Nsim = 1, Nchain = 1, u_init = NA){
  ## niter     : total number of iterations of the algorithm
  ## nburnin   : number of iterations for the burnin phase of the algorithm
  ##             (nburnin < niter)
  ## data      : name of the dataset, in the tibble format
  ##             It should at least contain the following variables
  ##             - "id": the ID of the curves (ie the ID for the plants)
  ##             - "varID": the ID of the varieties 
  ##             - "y": the response variable
  ##             - "time": the measurement times
  ## param     : list of initial values for the parameters (values obtained by 
  ##             maximum likelihood with the previous SAEM function recommended)
  ##             It should contain the following elements
  ##             - "beta": vector (fixed-effects)
  ##             - "P": matrix (covariance matrix of the random effects)
  ##             - "G": matrix (genetic covariance matrix)
  ##             - "sigma2": numeric (residual variance)
  ## GRM_mat   : genetic relationship matrix
  ## predictors: numbers of the columns in data containing the regression 
  ##             variables.  that enter in the definition of the curve model 
  ##             (eg time)
  ## Xi, Zi    : arrays containing the covariate matrices for each plant in the 
  ##             fixed-effects part and in the random effects part of the model
  ##             respectively
  ## *Unlike the variables in Xi and Zi, we are not interested in the effects of
  ## the variables informed through predictors. The variables in predictors are 
  ## structural variables that enter in the definition of the curve model, for 
  ## example time.* 
  ## model     : curve model function  
  ## Nchain    : number of Monte-Carlo chains used to simulate the plant parameters
  ## Nsim      : number of iterations for the Metropolis-Hastings algorithm used
  ##             to simulate the plant parameters
  ## u_init    : initialization of the genetic values
  
  
  ## Parameter values
  beta   <- param$beta
  P      <- param$P
  G      <- param$G
  sigma2 <- param$sigma2
  
  ## Data processing
  Na     <- dim(GRM_mat)[1]
  d      <- dim(Zi)[3]/dim(GRM_mat)[1] 
  nb.phi <- dim(Xi)[2]
  id     <- matrix(data$id)
  l_id   <- unique(id)
  y      <- data$y  
  time   <- data$time
  N      <- length(l_id)
  
  yi      <- list()
  timi    <- list()
  xij     <- list()
  Xilist  <- list()
  Zilist  <- list()
  tZilist <- list()
  
  for (i in 1:N){
    yi[[i]]      <- y[id == l_id[i]]
    timi[[i]]    <- time[id == l_id[i]]
    xij[[i]]     <- as.matrix(data[id == l_id[i],predictors])
    Xilist[[i]]  <- Xi[i,,]
    Zilist[[i]]  <- Zi[i,,]
    tZilist[[i]] <- t(Zi[i,,])
  }
  
  XX  <- do.call('rbind',Xilist)
  tZZ <- do.call('cbind',tZilist)
  ZZ  <- do.call('rbind',Zilist)
  
  ## Initialization or computation of recurrent quantities in the algorithm
  
  Xbetar     <- matrix(as.vector(XX%*%as.matrix(beta,ncol=1)),nrow=nb.phi*N)
  invP       <- solve(P)
  invG       <- solve(G)
  invA       <- solve(GRM_mat)
  invSigma_u <- (invA %x% invG)
  S_phi      <- P 
  u <- matrix(0,d*Na,niter+1)
  if (is.na(u_init)==T){
    u[,1] <- rep(0,d*Na) 
  } else {
    u[,1] <- u_init 
  }
  
  z_phi <- array(0,dim=c(niter+1,nb.phi,N,Nchain))
  
  ## Definition of the step sizes for the SA step of the algorithm
  ## The step size equals 1 during the burnin phase of the algorithm and then
  ## decreases
  
  gamma <- rep(1, niter)
  for (k in (nburnin + 1):niter) {
    gamma[k] <- (k - nburnin) ^ (-2 / 3)
  }
  
  ## To store the sufficient statistics iteration after iteration
  s1 <- array(0,dim=c(niter+1, d*Na,1))

  for (k in 1:niter){
    ## Loop of iterations of the algorithm
    
    ## S-step: draw a value \varphi from the conditional distribution of \varphi
    ## given Y and u.
  
    for (i in 1:N) {
      m_phi <- Xilist[[i]] %*% beta + Zilist[[i]] %*% matrix(u[,k])
      for (mc in 1:Nchain){
        current_phi <- z_phi[k, , i,mc]
        for (j in 1:Nsim) {
          phi      <- mvnfast::rmvn(1, mu = m_phi , sigma = S_phi) 
          mu_num   <- model(phi,xij[[i]]) 
          mu_denom <- model(current_phi,xij[[i]]) 
          lograte  <- sum(log(dnorm(yi[[i]], mu_num, sqrt(sigma2)))) - sum(log(dnorm(yi[[i]], mu_denom, sqrt(sigma2))))
          v <- runif(1)
          cond1 <- (log(v) <= min(0, lograte))
          if (cond1) {
            current_phi <- phi
          }
          z_phi[k + 1, , i,mc] <- current_phi
        } 
      } 
    }
    
    ## SA-step
    
    invPkbloc  <- bdiag_m(replicate(N,invP,simplify=FALSE))
    ZP         <- tZZ%*%invPkbloc
    sumZPZ     <- ZP%*%ZZ
    invSigma_u <- sumZPZ + (invA %x% invG)
    Sigma_u    <- solve(invSigma_u)
    sumZPdphi  <- 0
    for (mc in 1:Nchain){
      dphi      <- matrix(as.vector(z_phi[k + 1, , , mc])) - Xbetar
      sumZPdphi <- sumZPdphi + ZP%*%dphi
    }
    sumZPdphi <- sumZPdphi/Nchain
    s1[k+1,,] <- as.vector(s1[k,,] + gamma[k]*(sumZPdphi-s1[k,,]))
    
    ## M-step
    
    u[,k+1] <- as.vector(Sigma_u%*%as.matrix(s1[k+1,,],ncol=1))
  }
  
  res <- list(u=u, z_phi=z_phi)
}

