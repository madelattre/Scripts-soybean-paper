
##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################


####
## Study of the prediction accuracy based on simulated data
####

####
## CV0 : predictions on new environmental conditions (leave-one-environmental-out 
## scheme)  
## CV1 : predictions on unobserved varieties in any environment 
## CV2 : predictions on varieties not evaluated in one environment, but evaluated
## accross other different environments 
####


rm(list = ls())

## Source main functions and load necessary R libraries
library(tidyverse)
library(Matrix)
library(mvnfast)
library(ggpubr)

functionsrepository <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the algorithm functions repository 
resultsrepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the repository in which storing the results

source(paste(functionsrepository,"/main_functions.R",sep=""))

### Model specification
### -------------------

## Function for data simulation
Esp_logistic <- function(A,xmid,time){
  pred <- A/(1+exp(-0.1*(time-xmid))) 
  return(pred)
}  

## Function in the right format to be used in the SAEM functions
model_logistic <- function(phi, x){
  tim   <- as.numeric(x[,1])
  A     <- phi[1]; 
  xmid  <- phi[2]; 
  ypred <- A/(1+exp(-0.1*(tim-xmid)))
  return(ypred)
}


model <- model_logistic


nb.u   <- 2 # Number of genetic effects per variety
nb.phi <- 2 # Number of growth parameters per plot


### Parameter values
### ----------------

herit <- 1 #2

mu_A        <- 50
mu_x        <- 25
coef_cov_A  <- 25
coef_cov_x  <- 2 
beta        <- as.matrix(c(mu_A, mu_x, coef_cov_A, coef_cov_x))
G_true      <- diag(c(3,0.1))*herit 
P_true      <- diag(c(3,0.1))
sigma2_true <- 8
param.list  <- list(beta=c(mu_A,mu_x),sigma2=sigma2_true,P=P_true, G=G_true)


### Sample design
### -------------

time     <- seq(1,50,length.out=10) # time points
Nobs     <- length(time)            # number of observations per plant
Na       <- 225                     # total number of varieties
seq.env  <- c(1,2,3)                # environments
N        <- 2*length(seq.env)*Na    # total number of plots, two plots per combination of genotype and environment
Ntot     <- N*Nobs                  # total number of observations
env.var  <- rep(seq.env,each=N/(Na*length(seq.env))*length(time)) 
Nmarkers <- 10000
X        <- matrix(sample(c(-1,0,1), Nmarkers*Na, replace =T), Na, Nmarkers)
GRM_mat  <- X%*%t(X)/(Nmarkers-1) # genomic relationship matrix
Ntot <- N*Nobs

### Algorithmic settings
### --------------------

niterSAEM   <- 500
nburninSAEM <- 350

paraminit <- list(beta = c(40,30,20,0),
                 G=diag(rep(50,2)),
                 P=diag(rep(50,2)),
                 sigma2= 20)

niterPred   <- 500#1500
nburninPred <- 350#1300

### Experimental design
### -------------------


nbrep    <- 50 # number of simulated dataset
nb_group <- 10 # validation based on one tenth of the total dataset (ie training based on the other nine tenth)



## R objects for storing the results
param.est.cv2 <- matrix(0,nbrep,length(beta)+1+nb.u*(nb.u+1)/2+nb.phi*(nb.phi+1)/2)
param.est.cv0 <- matrix(0,nbrep,length(beta)+1+nb.u*(nb.u+1)/2+nb.phi*(nb.phi+1)/2)
param.est.cv1 <- matrix(0,nbrep,length(beta)+1+nb.u*(nb.u+1)/2+nb.phi*(nb.phi+1)/2)

true.u            <- matrix(0,nbrep,Na*nb.u)
u.pred.cv2  <- matrix(0,nbrep,Na*nb.u)
u.pred.cv0  <- matrix(0,nbrep,Na*nb.u)
u.pred.cv1  <- matrix(0,nbrep,(floor(Na/nb_group)+1)*nb.u)

rsd.y.cv2 <- matrix(0,nbrep,Nobs)
rsd.A.cv2 <- rep(0,nbrep)
rsd.x.cv2 <- rep(0,nbrep)

rsd.y.cv0 <- matrix(0,nbrep,Nobs)
rsd.A.cv0 <- rep(0,nbrep)
rsd.x.cv0 <- rep(0,nbrep)

rsd.y.cv1 <- matrix(0,nbrep,Nobs)
rsd.A.cv1 <- rep(0,nbrep)
rsd.x.cv1 <- rep(0,nbrep)

data.cv2 <- list()
data.cv0 <- list()
data.cv1 <- list()



for (rep in 1:nbrep){
  
  print(rep)
  
  set.seed(rep*1000)
  
  ## A. Data simulation
  
  
  data <- tibble (id = rep(1:N, each=Nobs), varID = rep(1:Na, each= Nobs*N/Na), time= rep(time, N),
                  env = rep(env.var,Na))
  true.u[rep,] <- u <- mvtnorm::rmvnorm(n = 1, mean=rep(0, nb.u*Na), sigma = GRM_mat%x%G_true)
  
  p            <- mvtnorm::rmvnorm(n = N, mean=rep(0, 2), sigma = P_true)   
  data$epsilon <- rnorm(Ntot, 0, sd=sqrt(sigma2_true))
  data$phi_A     <- NA
  data$phi_x     <- NA
  data$y       <- NA
  data$GE      <- NA
  l_var        <- seq(1,Na)
  Zi           <- array(0,dim=c(N,nb.phi,nb.u*Na))
  Xi           <- array(0,dim=c(N,nb.phi,length(beta)))
  
  for (i in 1:N){
    envi                          <- data$env[(i-1)*Nobs+1]
    vari                          <- data$varID[(i-1)*Nobs+1]
    Xi[i,,]                       <- matrix(c(1,0,envi,0,0,1,0,envi),nrow=nb.phi,byrow=TRUE)
    Zi[i,,]                       <- (diag(1,nb.phi)%x%t(l_var==vari)) 
    Ziu                           <- Zi[i,,]%*%t(u)
    Xibeta                        <- Xi[i,,]%*%beta
    A_i                           <- Xibeta[1] + Ziu[1] + p[i,1]
    x_i                           <- Xibeta[2] + Ziu[2] + p[i,2]
    data$phi_A[(i-1)*Nobs+1:Nobs] <- A_i
    data$phi_x[(i-1)*Nobs+1:Nobs] <- x_i
    data$y[(i-1)*Nobs+1:Nobs]     <- Esp_logistic(A=A_i, xmid=x_i, time=data$time[(i-1)*Nobs+1:Nobs])+ data$epsilon[(i-1)*Nobs+1:Nobs]
    data$GE[(i-1)*Nobs+1:Nobs]    <- paste(vari,envi,sep="_")
  }
  
  # res.ref <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, 
  #                     data=data,predictors=c(3), paraminit=paraminit,
  #                     GRM_mat = GRM_mat, Xi=Xi, Zi=Zi, model = model)
  # 
  # param.est.ref <- list(beta=res.ref$beta[niterSAEM+1,],
  #                       sigma2=res.ref$sigma2[niterSAEM+1],
  #                       P=res.ref$P[niterSAEM+1,,],
  #                       G=res.ref$G[niterSAEM+1,,])
  
  # plot_variety <- ggplot(data,
  #                        aes(x=time, y=y, group=id,
  #                            color=as.factor(env))) +
  #   geom_line() + theme_bw() +
  #   theme(legend.position = "none",
  #         axis.title.y = element_text(size=16),
  #         axis.title.x = element_text(size=16),
  #         axis.text.y = element_text(size=12),
  #         axis.text.x = element_text(size=12, angle = 45)) +
  #   ylab("Height") + xlab("Heat time")
  # 
  # 
  # plot_variety
  
  
  ## CV2 : prediction of incomplete environmental data
  ## Leave some combinations of genotypes and environments out (but all genotypes
  ## and environments should be in a training set) and predict them.
  
  l_GE <- select(data, GE) %>% distinct() 
  n_GE <- dim(l_GE)[1]
  
  var_train.cv2 <- tibble(varID=NA)
  
  while (length(var_train.cv2$varID)!=Na){ 
    
    repartition.cv2 <- sample(rep(1:nb_group, length=(n_GE)))
    groups.cv2      <- tibble(GE=l_GE$GE,group=repartition.cv2)
    datatry.cv2     <- left_join(data,groups.cv2,"GE")
  

    # All varieties must be observed in the training set
    training.cv2 <- datatry.cv2 %>% filter(group!=10)
    var_train.cv2  <- select(training.cv2, "varID") %>% distinct()
    l_id_train.cv2 <- select(training.cv2, "id") %>% distinct()
    validation.cv2 <- datatry.cv2 %>% filter(group==10)
    l_id_valid.cv2 <- select(validation.cv2, "id") %>% distinct()
  }
    
   
  res.train.cv2 <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, 
                            data=training.cv2,predictors=c(3), paraminit=paraminit,
                            GRM_mat = GRM_mat, Xi=Xi[l_id_train.cv2$id,,], 
                            Zi=Zi[l_id_train.cv2$id,,], model = model) 
  
  param.est.cv2.list <- list(beta=res.train.cv2$beta[niterSAEM+1,],
                          sigma2=res.train.cv2$sigma2[niterSAEM+1],
                          P=res.train.cv2$P[niterSAEM+1,,],
                          G=res.train.cv2$G[niterSAEM+1,,])
  
  param.est.cv2[rep,] <- c(res.train.cv2$beta[niterSAEM+1,],
                        res.train.cv2$sigma2[niterSAEM+1],
                        res.train.cv2$P[niterSAEM+1,1,1],
                        res.train.cv2$P[niterSAEM+1,2,2],
                        res.train.cv2$P[niterSAEM+1,1,2],
                        res.train.cv2$G[niterSAEM+1,1,1],
                        res.train.cv2$G[niterSAEM+1,2,2],
                        res.train.cv2$G[niterSAEM+1,1,2]
                        )
  
  res.pred.cv2 <- predict.u(niter=niterPred, nburnin=nburninPred, 
                             data=training.cv2, predictors=c(3),
                             param=param.est.cv2.list, GRM_mat=GRM_mat,
                             Xi=Xi[l_id_train.cv2$id,,], 
                             Zi=Zi[l_id_train.cv2$id,,], model=model,
                             Nsim=10, Nchain = 10)
  
  u.pred.cv2[rep,] <- res.pred.cv2$u[,niterPred+1]
  
  # Evaluation of the predictions
  
  validation.cv2$phi_A_cv2 <- NA
  validation.cv2$phi_x_cv2 <- NA
  validation.cv2$y_cv2     <- NA
  
  for (i in 1:length(l_id_valid.cv2$id)){
    envi                                        <- validation.cv2$env[(i-1)*Nobs+1]
    vari                                        <- validation.cv2$varID[(i-1)*Nobs+1]
    X                                           <- matrix(c(1,0,envi,0,0,1,0,envi),
                                                          nrow=nb.phi,byrow=TRUE)
    Z                                           <- (diag(1,nb.phi)%x%t(l_var==vari)) 
    Zu                                          <- Z%*%u.pred.cv2[rep,]
    Xbeta                                       <- X%*%param.est.cv2.list$beta
    A_i_cv2                                     <- Xbeta[1] + Ziu[1] 
    x_i_cv2                                     <- Xbeta[2] + Ziu[2]
    validation.cv2$phi_A_cv2[(i-1)*Nobs+1:Nobs] <- A_i_cv2
    validation.cv2$phi_x_cv2[(i-1)*Nobs+1:Nobs] <- x_i_cv2
    validation.cv2$y_cv2[(i-1)*Nobs+1:Nobs]     <- Esp_logistic(A=A_i_cv2, xmid=x_i_cv2, 
                                                   time=validation.cv2$time[(i-1)*Nobs+1:Nobs])
  }
  
  # Compute kind of RSD, time point by time point, can then be averaged over time
  # if necessary 
  
  time.cv2 <- unique(validation.cv2$time) 
  
  for (j in 1:length(time.cv2)){
    tj <- which(validation.cv2$time==time.cv2[j])
    rsd.y.cv2[rep,j] <- mean(abs((validation.cv2$y[tj]-validation.cv2$y_cv2[tj])/validation.cv2$y[tj]))
  }
  
  rsd.A.cv2[rep] <- mean(abs((validation.cv2$phi_A[tj]-validation.cv2$phi_A_cv2[tj])/validation.cv2$phi_A[tj]))
  rsd.x.cv2[rep] <- mean(abs((validation.cv2$phi_x[tj]-validation.cv2$phi_x_cv2[tj])/validation.cv2$phi_x[tj]))
  
  data.cv2[[rep]] <- validation.cv2 %>% select(-c(7,9,10))
  
  ## CV0 : prediction of  tested genotypes in untested environments
  ## Leave one environment out and predict the phenotypes of all genotypes in 
  ## the left-out environment.
  
  env.valid    <- sample(seq.env,1)
  training.cv0 <- data %>% filter(env!=env.valid)
  var_train.cv0  <- select(training.cv0, "varID") %>% distinct()
  l_id_train.cv0 <- select(training.cv0, "id") %>% distinct()
  
  validation.cv0 <- data %>% filter(env==env.valid)
  l_id_valid.cv0 <- select(validation.cv0, "id") %>% distinct()
  
  
  res.train.cv0 <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, 
                            data=training.cv0,predictors=c(3), paraminit=paraminit,
                            GRM_mat = GRM_mat, Xi=Xi[l_id_train.cv0$id,,], 
                            Zi=Zi[l_id_train.cv0$id,,], model = model) 
  
  param.est.cv0.list <- list(beta=res.train.cv0$beta[niterSAEM+1,],
                        sigma2=res.train.cv0$sigma2[niterSAEM+1],
                        P=res.train.cv0$P[niterSAEM+1,,],
                        G=res.train.cv0$G[niterSAEM+1,,])
  
  param.est.cv0[rep,] <- c(res.train.cv0$beta[niterSAEM+1,],
                           res.train.cv0$sigma2[niterSAEM+1],
                           res.train.cv0$P[niterSAEM+1,1,1],
                           res.train.cv0$P[niterSAEM+1,2,2],
                           res.train.cv0$P[niterSAEM+1,1,2],
                           res.train.cv0$G[niterSAEM+1,1,1],
                           res.train.cv0$G[niterSAEM+1,2,2],
                           res.train.cv0$G[niterSAEM+1,1,2])
 
  res.pred.cv0 <- predict.u(niter=niterPred, nburnin=nburninPred, 
                             data=training.cv0, predictors=c(3),
                             param=param.est.cv0.list, GRM_mat=GRM_mat,
                             Xi=Xi[l_id_train.cv0$id,,], 
                             Zi=Zi[l_id_train.cv0$id,,], model=model,
                             Nsim=10, Nchain = 10)
  
  u.pred.cv0[rep,] <- res.pred.cv0$u[,niterPred+1]
  
  # Evaluation of the predictions
  
  validation.cv0$phi_A_cv0 <- NA
  validation.cv0$phi_x_cv0 <- NA
  validation.cv0$y_cv0     <- NA
  
  for (i in 1:length(l_id_valid.cv0$id)){
    envi                                        <- validation.cv0$env[(i-1)*Nobs+1]
    vari                                        <- validation.cv0$varID[(i-1)*Nobs+1]
    X                                           <- matrix(c(1,0,envi,0,0,1,0,envi),
                                                          nrow=nb.phi,byrow=TRUE)
    Z                                           <- (diag(1,nb.phi)%x%t(l_var==vari)) 
    Zu                                          <- Z%*%u.pred.cv0[rep,]
    Xbeta                                       <- X%*%param.est.cv0.list$beta
    A_i_cv0                                     <- Xbeta[1] + Ziu[1] 
    x_i_cv0                                     <- Xbeta[2] + Ziu[2]
    validation.cv0$phi_A_cv0[(i-1)*Nobs+1:Nobs] <- A_i_cv0
    validation.cv0$phi_x_cv0[(i-1)*Nobs+1:Nobs] <- x_i_cv0
    validation.cv0$y_cv0[(i-1)*Nobs+1:Nobs]     <- Esp_logistic(A=A_i_cv0, xmid=x_i_cv0, 
                                                                time=validation.cv0$time[(i-1)*Nobs+1:Nobs])
  }
  
  # Compute kind of RSD, time point by time point, can then be averaged over time
  # if necessary 
  
  time.cv0 <- unique(validation.cv0$time) 
  
  for (j in 1:length(time.cv0)){
    tj <- which(validation.cv0$time==time.cv0[j])
    rsd.y.cv0[rep,j] <- mean(abs((validation.cv0$y[tj]-validation.cv0$y_cv0[tj])/validation.cv0$y[tj]))
  }
  
  rsd.A.cv0[rep] <- mean(abs((validation.cv0$phi_A[tj]-validation.cv0$phi_A_cv0[tj])/validation.cv0$phi_A[tj]))
  rsd.x.cv0[rep] <- mean(abs((validation.cv0$phi_x[tj]-validation.cv0$phi_x_cv0[tj])/validation.cv0$phi_x[tj]))
  
  data.cv0[[rep]] <- validation.cv0 %>% select(-c(7,9,10))
  
  ## CV1 :
  ## Leave some genotypes out and predict the phenotypes of left-out genotypes 
  ## in 3 environments.

  repartition.cv1 <- sample(rep(1:nb_group, length=(Na)))
  groups.cv1      <- tibble(varID=seq(1,Na),group=repartition.cv1)
  datatry.cv1     <- left_join(data,groups.cv1,"varID")
  
  training.cv1 <- datatry.cv1 %>% filter(group!=1)
  var_train.cv1  <- select(training.cv1, "varID") %>% distinct()
  l_id_train.cv1 <- select(training.cv1, "id") %>% distinct()
  
  GRM_mat.cv1 <- GRM_mat[var_train.cv1$varID,var_train.cv1$varID]
  GRM_mat.cv1.cross <- GRM_mat[-var_train.cv1$varID,var_train.cv1$varID]

  Nvar.train.cv1 <- length(var_train.cv1$varID)
  N.train.cv1 <- length(l_id_train.cv1$id)
  
  
  Zi.temp           <- array(0,dim=c(N.train.cv1,nb.phi,nb.u*Na))
  Xi.cv1           <- array(0,dim=c(N.train.cv1,nb.phi,length(beta)))
  
  for (i in 1:N.train.cv1){
    id <- l_id_train.cv1$id[i]
    Xi.cv1[i,,] <- Xi[id,,]
    Zi.temp[i,,] <- Zi[id,,]
  }
  
  Zi.cv1 <- Zi.temp[,,which(apply(Zi.temp,3,sum)>0)]
  
  validation.cv1 <- data %>% filter(!(id %in% l_id_train.cv1$id))
  l_id_valid.cv1 <- select(validation.cv1, "id") %>% distinct()
  l_var_valid.cv1 <- select(validation.cv1, "varID") %>% distinct()
  
  
  res.train.cv1 <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, 
                            data=training.cv1,predictors=c(3), 
                            paraminit=paraminit,
                            GRM_mat = GRM_mat.cv1, Xi=Xi.cv1, 
                            Zi=Zi.cv1, model = model) 
  
  param.est.cv1.list <- list(beta=res.train.cv1$beta[niterSAEM+1,],
                        sigma2=res.train.cv1$sigma2[niterSAEM+1],
                        P=res.train.cv1$P[niterSAEM+1,,],
                        G=res.train.cv1$G[niterSAEM+1,,])
  
  param.est.cv1[rep,] <- c(res.train.cv1$beta[niterSAEM+1,],
                           res.train.cv1$sigma2[niterSAEM+1],
                           res.train.cv1$P[niterSAEM+1,1,1],
                           res.train.cv1$P[niterSAEM+1,2,2],
                           res.train.cv1$P[niterSAEM+1,1,2],
                           res.train.cv1$G[niterSAEM+1,1,1],
                           res.train.cv1$G[niterSAEM+1,2,2],
                           res.train.cv1$G[niterSAEM+1,1,2])

  res.pred.cv1 <- predict.u(niter=niterPred, nburnin=nburninPred, 
                             data=training.cv1, predictors=c(3),
                             param=param.est.cv1.list, GRM_mat=GRM_mat.cv1,
                             Xi=Xi.cv1, Zi=Zi.cv1, model=model,
                             Nsim=10, Nchain = 10)
  
  MM <- (GRM_mat.cv1.cross %*% solve(GRM_mat.cv1)) %x% diag(nb.u)
  
  u.pred.cv1[rep,] <- MM%*%matrix(res.pred.cv1$u[,niterPred+1],ncol=1)
  
 
  # Evaluation of the predictions
  
  validation.cv1$phi_A_cv1 <- NA
  validation.cv1$phi_x_cv1 <- NA
  validation.cv1$y_cv1     <- NA
  
  for (i in 1:length(l_id_valid.cv1$id)){
    envi                                        <- validation.cv1$env[(i-1)*Nobs+1]
    vari                                        <- validation.cv1$varID[(i-1)*Nobs+1]
    X                                           <- matrix(c(1,0,envi,0,0,1,0,envi),
                                                          nrow=nb.phi,byrow=TRUE)
    Z                                           <- (diag(1,nb.phi)%x%t(l_var_valid.cv1$varID==vari)) 
    Zu                                          <- Z%*%u.pred.cv1[rep,]
    Xbeta                                       <- X%*%param.est.cv1.list$beta
    A_i_cv1                                     <- Xbeta[1] + Ziu[1] 
    x_i_cv1                                     <- Xbeta[2] + Ziu[2]
    validation.cv1$phi_A_cv1[(i-1)*Nobs+1:Nobs] <- A_i_cv1
    validation.cv1$phi_x_cv1[(i-1)*Nobs+1:Nobs] <- x_i_cv1
    validation.cv1$y_cv1[(i-1)*Nobs+1:Nobs]     <- Esp_logistic(A=A_i_cv1, xmid=x_i_cv1, 
                                                                time=validation.cv1$time[(i-1)*Nobs+1:Nobs])
  }
  
  # Compute kind of RSD, time point by time point, can then be averaged over time
  # if necessary 
  
  time.cv1 <- unique(validation.cv1$time) 
  
  for (j in 1:length(time.cv1)){
    tj <- which(validation.cv1$time==time.cv1[j])
    rsd.y.cv1[rep,j] <- mean(abs((validation.cv1$y[tj]-validation.cv1$y_cv1[tj])/validation.cv1$y[tj]))
  }
  
  rsd.A.cv1[rep] <- mean(abs((validation.cv1$phi_A[tj]-validation.cv1$phi_A_cv1[tj])/validation.cv1$phi_A[tj]))
  rsd.x.cv1[rep] <- mean(abs((validation.cv1$phi_x[tj]-validation.cv1$phi_x_cv1[tj])/validation.cv1$phi_x[tj]))
  
  data.cv1[[rep]] <- validation.cv1 %>% select(-c(7,9,10))
  
  ## Store the results
  res.pred <- list(param.est.cv2=param.est.cv2, param.est.cv0=param.est.cv0,
                   param.est.cv1=param.est.cv1, data.cv2=data.cv2, data.cv0=data.cv0,
                   data.cv1=data.cv1, rsd.A.cv2=rsd.A.cv2, rsd.A.cv0=rsd.A.cv0,
                   rsd.A.cv1=rsd.A.cv1, rsd.x.cv2=rsd.x.cv2, rsd.x.cv0=rsd.x.cv0,
                   rsd.x.cv1=rsd.x.cv1, rsd.y.cv2=rsd.y.cv2, rsd.y.cv0=rsd.y.cv0,
                   rsd.y.cv1=rsd.y.cv1, u.pred.cv2=u.pred.cv2, u.pred.cv0=u.pred.cv0,
                   u.pred.cv1=u.pred.cv1)
  
  
  
  save(res.pred,file="Simu_pred1.Rdata")
  
}
  
