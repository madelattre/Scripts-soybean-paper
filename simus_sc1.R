##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################


## Script of simulations under Scenario 1 in section 2.4 of the article

## Source main functions and load necessary R libraries
library(tidyverse)
library(Matrix)
library(mvnfast)

functionsrepository <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the algorithm functions repository 
resultsrepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the repository in which storing the results

source(paste(functionsrepository,"/main_functions.R",sep=""))

#getwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Model functions

### For data simulation
Esp_logistic <- function(A,xmid,time){
  scal <- 10; 
  pred <- A/(1+exp(-0.1*(time-xmid)))
  return(pred)
}  
### In the right format to be used in the SAEM functions
model_logistic <- function(phi, x){
  tim   <- as.numeric(x[,1])
  A     <- phi[1]; 
  xmid  <- phi[2]; 
  ypred <- A/(1+exp(-0.1*(tim-xmid)))
  return(ypred)
}


model <- model_logistic

## Parameter values

d <- 2

mu_A        <- 50
mu_x        <- 13
beta        <- c(mu_A, mu_x)
G_true      <- diag(c(25,1.5))
P_true      <- diag(c(6,0.4))
sigma2_true <- 4
param.list  <- list(beta=c(mu_A,mu_x),sigma2=sigma2_true,P=P_true, G=G_true)


time <- seq(1,45,1)  # Observation time points 
Nobs <- length(time) # Number of observations per plant
Nmarkers <- 10000    # Number of genetic markers for the simulation of the GRM matrix


## Algorithmic settings
###  For parameter estimation
niterSAEM   <- 1500
nburninSAEM <- 1300
### For genetic effects prediction
niterPred   <- 1500
nburninPred <- 1300

## Number of simulated dataset
nbrep <- 50

N <- 2000 ## Total number of plants
for (Na in c(50,100)){ ## Numbers of varieties used for the numerical experiment
  
    ## R objects for storing the results
    param.est         <- matrix(0,nbrep,d+1+d*(d+1))
    true.u            <- matrix(0,nbrep,Na*d)
    total.SAEM.time   <- rep(0,nbrep)
    total.Pred.time   <- rep(0,nbrep)
    u.pred            <- matrix(0,nbrep,Na*d)
    u.pred.theta.true <- matrix(0,nbrep,Na*d)
    
    ## Simulation of the genomic relationship matrix (GRM)
    X       <- matrix(sample(c(-1,0,1), Nmarkers*Na, replace =T), Na, Nmarkers)
    GRM_mat <- X%*%t(X)/(Nmarkers-1) 
    
    filename <- paste(resultsrepositoryname,'/ResSimus_Sc1_Na',Na,'_N',N,'.Rdata',sep="")
    
    for (rep in 1:nbrep){
      print(rep)
      ## Loop of repetitions of the experiment

      ## Data Simulation
      
      Ntot <- N*Nobs
      data <- tibble (id = rep(1:N, each=Nobs), varID = rep(1:Na, each= Nobs*N/Na), time= rep(time, N))
      ## 1. Simulation of the zero-mean random effects
      true.u[rep,] <- u <- mvtnorm::rmvnorm(n = 1, mean=rep(0, d*Na), sigma = GRM_mat%x%G_true)
      p            <- mvtnorm::rmvnorm(n = N, mean=rep(0, 2), sigma = P_true)   
      ## 2. Simulation of the plant parameters and of the observations
      data$phi_A   <- NA
      data$phi_x   <- NA
      data$epsilon <- rnorm(Ntot, 0, sd=sqrt(sigma2_true))
      data$y       <- NA
      l_var        <- seq(1,Na)
      Zi           <- array(0,dim=c(N,d,d*Na))
      Xi           <- array(0,dim=c(N,d,d))
      
      for (i in 1:N){
        vari                          <- data$varID[(i-1)*Nobs+1]
        Xi[i,,]                       <- diag(d)
        Zi[i,,]                       <- (diag(d)%x%t(l_var==vari))
        Ziu                           <- Zi[i,,]%*%t(u)
        A_i                           <- mu_A + Ziu[1] + p[i,1]
        x_i                           <- mu_x + Ziu[2] + p[i,2]
        data$phi_A[(i-1)*Nobs+1:Nobs] <- A_i
        data$phi_x[(i-1)*Nobs+1:Nobs] <- x_i
        data$y[(i-1)*Nobs+1:Nobs]     <- Esp_logistic(A=A_i, xmid=x_i, time=data$time[(i-1)*Nobs+1:Nobs])+ data$epsilon[(i-1)*Nobs+1:Nobs]
      }
      
      ## Running the algorithms for parameter estimation and prediction
  
      ## 1. Parameter estimation
      
      ## Initialization
      paraminit <- list(beta = matrix(c(runif(1,mu_A-10,mu_A+10),
                                        runif(1,mu_x-4,mu_x+4)),d,1),
                        G=diag(rep(50,d)), 
                        P=diag(rep(50,d)), 
                        sigma2= 20)
      
      
      TimeSAEM1 <- Sys.time()
      
      ## SAEM call
      resSAEM <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, data=data, 
                          predictors=c(3), paraminit=paraminit, GRM_mat = GRM_mat, 
                          Xi=Xi, Zi=Zi, model = model)
 
      TimeSAEM2 <- Sys.time()
      
      total.SAEM.time[rep] <- TimeSAEM2-TimeSAEM1
      
      ## Store the result
      param.est[rep,] <- c(resSAEM$beta[niterSAEM+1,],resSAEM$sigma2[niterSAEM+1],
                           resSAEM$P[niterSAEM+1,1,],resSAEM$P[niterSAEM+1,2,2],
                           resSAEM$G[niterSAEM+1,1,],resSAEM$G[niterSAEM+1,2,2])
      
      # Rearrange the parameters so that the prediction function is able to use them
      
      param.est.list <- list(beta=resSAEM$beta[niterSAEM+1,],sigma2=resSAEM$sigma2[niterSAEM+1],
                             P=resSAEM$P[niterSAEM+1,,], G=resSAEM$G[niterSAEM+1,,])

      ## 2. Genetic effects predictions

      TimePred1 <- Sys.time()
  
      ## Algorithm call: prediction obtained with the estimated parameter obtained 
      ## with SAEM
      res.pred <- predict.u(niterPred, nburninPred, data=data, predictors=c(3), 
                            param=param.est.list, GRM_mat = GRM_mat, Xi=Xi, 
                            Zi=Zi, model = model, Nsim=10, Nchain=5) 
      TimePred2 <- Sys.time()
      
      total.Pred.time[rep] <- TimePred2-TimePred1
      
      ## Algorithm call: prediction obtained with the true parameter value 
      res.pred.theta.true <- predict.u(niterPred, nburninPred, data=data, 
                                       predictors=c(3), param=param.list, 
                                       GRM_mat = GRM_mat, Xi=Xi, Zi=Zi, 
                                       model = model, Nsim=10, Nchain=5) 
      
      ## Store the results
      u.pred[rep,]            <- res.pred$u[,niterPred+1]
      u.pred.theta.true[rep,] <- res.pred.theta.true$u[,niterPred+1]
      
      full.res <- list(param.est=param.est,true.u=true.u,u.pred=u.pred,
                       u.pred.theta.true=u.pred.theta.true,
                       total.SAEM.time=total.SAEM.time,total.Pred.time=total.Pred.time)
      
      save(full.res,file=filename)
            
    }
    
  }


