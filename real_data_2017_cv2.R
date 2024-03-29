##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################


## Script for cross-validation on the 2017 soybean growth data under scenario
## CV2 : predictions on varieties not evaluated in one environment, but evaluated
## accross other different environments 

rm(list=ls())

## Load necessary R libraries 

library(Matrix)
library(mvtnorm)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(doBy)
library(readxl)
library(ggpubr)

## Source main functions and specify repositories
functionsrepository <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the algorithm functions repository 
datarepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the data repository
resultsrepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the repository in which storing the results

source(paste(functionsrepository,"/main_functions.R",sep=""))

####
## 1- Data processing
####

## Load GRM matrix and processed data 

grm <- read.csv(paste(datarepositoryname,"/grm.csv",sep=""))
A <- data.matrix(grm[,-1])
data <- read.csv(paste(datarepositoryname,"/pheno.csv",sep=""))[,-1]
data.log <- data %>% mutate(y = log(y))


## Data characteristics

Nv    <- ncol(A) ## Number of varieties
l_id  <- unique(data.log$id) ## Sequence of plant identifiers
Ntot  <- nrow(data.log) ## Total number of observations
N     <- length(l_id) ## Number of plants
l_var <- select(data.log, varID, variety) %>% distinct() 
var_A <- tibble(order = 1:Nv,variety = colnames(A))
l_var <- left_join(l_var, var_A, by="variety") %>% arrange(order) ## Identifiers and names of the varieties
rn    <- rownames(A)

## Model definition 

asym.model <- function(param, x){
  tim   <- as.numeric(x[,1])
  A     <- param[1];
  B     <- param[2];
  ypred <- A*(1-exp(-exp(B)*tim))
  return(ypred)
}

## Definition of the covariate matrices in order to consider an effect of the 
## water condition and of its interaction with the genome on the two growth 
## parameters

nb.phi <- 2
d      <- 2*nb.phi
Zi     <- array(0,dim=c(N,nb.phi,d*Nv))
Xi     <- array(0,dim=c(N,nb.phi,2*nb.phi)) 

for (i in 1:N){
  condi <- 1-as.numeric(data.log$condD[data.log$id==l_id[i]][1])
  Xi[i,1:nb.phi,1:nb.phi] <- diag(nb.phi)
  Xi[i,1:nb.phi,(nb.phi+1):(2*nb.phi)] <- diag(condi,nb.phi)
  vari <- data.log$varID[data.log$id==l_id[i]][1]
  if (condi==0){
    Zi[i,1:nb.phi,1:(nb.phi*Nv)] <- diag(nb.phi)%x%t(l_var$varID==vari)
  } else{
    Zi[i,1:nb.phi,(nb.phi*Nv+1):(d*Nv)] <- diag(nb.phi)%x%t(l_var$varID==vari)
  }
}


## Initialization and algorithmic settings 

paraminit.log <- list(beta = matrix(c(3,-5,1,0),2*nb.phi,1),G=diag(rep(5,d)),
                      P=diag(c(5,5)) , sigma2= 10)
niterSAEM     <- 400
nburninSAEM   <- 300
niterPred     <- 1000
nburninPred   <- 800

nrep_CV  <- 4  # Nb of repetitions for cross-validation
nb_group <- 5  # Nb of samples for cross-validation

## Cross-validation


data.log <- data.log %>% mutate(GE = paste(varID,condition,sep="_"))
#l_GE     <- select(data.log, GE) %>% distinct() 
#n_GE     <- dim(l_GE)[1]

for (repCV in 1:nrep_CV){
  
  repartitionC <- sample(rep(1:nb_group, length=Nv))
  g2           <- sample(seq(1:(nb_group-1)),size=1)
  repartitionD <- (repartitionC + (g2-1))%%nb_group + 1 
  groups       <- tibble(GE=c(paste(l_var$varID,"C",sep="_"),paste(l_var$varID,"D",sep="_")),
                         group=c(repartitionC,repartitionD))
  datatry      <- left_join(data.log,groups,"GE")
  
  ## Elements to save the different results
  
  param.est.cv <- list()
  u.train.cv <- list()
  u.valid.cv <- list()
  val.data <- list()
  
  for (g in 1:nb_group){
    
    
    ## Definition of the training and validation subsets
    
    validation <- datatry %>% filter(group==g)
    var_valid  <- unique(validation$variety)
    rn_valid   <- which(rn %in% var_valid)
    A_valid    <- A[rn_valid,rn_valid]
    training   <- datatry %>% filter(group!=g)
    var_train  <- unique(training$variety)
    rn_train   <- which(rn %in% var_train)
    A_train    <- A[rn_train,rn_train]
    
    
    
    
    ## Processing of the training and validation data
    
    Nvtrain <- length(unique(training$varID))
    Ntrain  <- length(unique(training$id))
    
    l_id_train  <- unique(training$id)
    l_var_train <- dplyr::select(training, "varID", "variety") %>% distinct()
    
    Nvvalid <- length(unique(validation$varID))
    Nvalid  <- length(unique(validation$id))
    
    l_id_valid <- unique(validation$id)
    l_var_valid <- dplyr::select(validation, "varID", "variety") %>% distinct()
    
    ## Estimation of the population parameters based on the training subset
    
    param.est.cv[[g]] <- res.train <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, data=training,
                                               predictors=c(4,5), paraminit=paraminit.log,
                                               GRM_mat = A_train, #Xi=Xi.train, Zi=Zi.train,
                                               Xi=Xi[which(l_id %in% l_id_train),,], 
                                               Zi=Zi[which(l_id %in% l_id_train),,],
                                               model = asym.model, Nchain=10)
    
    param.est.train <- list(beta=res.train$beta[niterSAEM+1,],
                            sigma2=res.train$sigma2[niterSAEM+1],
                            P=res.train$P[niterSAEM+1,,],
                            G=res.train$G[niterSAEM+1,,])
    

    res.pred <- predict.u(niter=niterPred, nburnin=nburninPred, data=training,
                          predictors=c(4,5),param.est.train, GRM_mat=A_train,
                          Xi=Xi[which(l_id %in% l_id_train),,], 
                          Zi=Zi[which(l_id %in% l_id_train),,],
                          model=asym.model,Nsim=10, 
                          Nchain = 10)
    
    u.train.cv[[g]] <- res.pred$u[,niterPred+1]
    
    # As we do not know the true plant growth parameters, we can only compare
    # the predicted y with the observed values
    
    validation$ypred <- NA
    for (i in 1:length(l_id_valid)){
      idi                     <- which(l_id==l_id_valid[i])
      index                   <- which(validation$id==l_id_valid[i])
      timi                    <- validation$time[index]
      X                       <- Xi[idi,,]
      Z                       <- Zi[idi,,]
      Zu                      <- Z%*%u.train.cv[[g]]
      Xbeta                   <- X%*%param.est.train$beta
      Ai                      <- Xbeta[1] + Zu[1] 
      Bi                      <- Xbeta[2] + Zu[2]
      validation$ypred[index] <- asym.model(param=c(Ai,Bi), x=matrix(timi,ncol=1))
    }
    
    val.data[[g]] <- validation
    
    ## Store the results
    
    res.CV <- list(param.est.cv,u.train.cv,u.valid.cv,val.data)
    save(res.CV,file=paste(resultsrepositoryname,"/Soybean_data_2017_cv2_",repCV,".Rdata",sep=""))
  }
}



