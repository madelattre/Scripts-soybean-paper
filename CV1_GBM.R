##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################


## Script for cross-validation on the 2017 soybean growth data under scenario
## CV1 (prediction on unobserved varieties in any environment) using the GBM 
## method (Onogi et al., 2020) 

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
library(GenomeBasedModel)

## Source main functions and specify repositories
functionsrepository <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the algorithm functions repository 
datarepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the data repository
resultsrepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the repository in which storing the results

source(paste(functionsrepository,"/main_functions.R",sep=""))

## Load GRM matrix and processed data 

grm <- read.csv(paste(datarepositoryname,"/grm.csv",sep=""))
A <- data.matrix(grm[,-1])
data <- read.csv(paste(datarepositoryname,"/pheno.csv",sep=""))[,-1]
data.log <- data %>% mutate(y = log(y))

## Model function in appropriate format for GenomeBasedModel

asym.model <- function(input, freevec, parameters){
  A <- parameters[1]
  B <- parameters[2]
  ypred <- A * (1 - exp(- exp(B) * input))
  return(ypred)
}


nrep_CV  <- 4  # Nb of repetitions for cross-validation
nb_group <- 5  # Nb of samples for cross-validation


for (repCV in 1:nrep_CV){
  valid.pred <- list()
  
  # load data to reconstitute the same training and validation samples as with 
  # our method
  load(paste(datarepositoryname,"/Soybean_data_2017_cv1_",repCV,".Rdata",sep="")) 
  for (g in 1:nb_group){
    # reconstitute g-th training and validation samples in the repCV-th repetition 
    # of cross-validation
    val.data   <- inner_join(data.log,res.CV[[4]][[g]],by=c("id","time","y","varID","condD","variety")) 
    train.data <- anti_join(data.log,val.data) 
    # prepare the training sample to use the GenomeBasedModel function
    data.train.sorted <- train.data %>% arrange(id, time)
    tmp.train         <- data.train.sorted %>% mutate(Date=time) %>% 
      select(id, Date, time) %>% spread(Date,time)
    Input     <- t(as.matrix(tmp.train[,-1]))
    Freevec   <- 0 
    tmp.train <- data.train.sorted %>% mutate(Date=time) %>% select(id, Date, y) %>% spread(Date, y)
    Y.train   <- t(tmp.train)
    tmp.train <- data.train.sorted %>% select(id, variety, condD, time, y) %>% spread(time, y)
    Geno.train <- rbind(tmp.train$id, rep(1, nrow(tmp.train))) 
    Q.train <- rbind(rep(1, nrow(tmp.train)), as.numeric(tmp.train$condD))
    
    # prepare useful vectors and dataframes to keep informations about plot 
    # identifyers and varieties
    var.C.train <- tmp.train$variety[!tmp.train$condD]
    var.D.train <- tmp.train$variety[tmp.train$condD]
    
    id.C.train <- tmp.train$id[!tmp.train$condD]
    id.D.train <- tmp.train$id[tmp.train$condD]
    
    C.train <- data.frame(id=id.C.train,variety=var.C.train)
    D.train <- data.frame(id=id.D.train,variety=var.D.train)
    
    # prepare matrix K for the training sample to be used in the GenomeBasedModel
    K.t <- rbind(cbind(A[var.C.train, var.C.train], matrix(0, length(var.C.train), length(var.D.train))),
                 cbind(matrix(0, length(var.D.train), length(var.C.train)), A[var.D.train, var.D.train]))
    K.t <- K.t + diag(0.00001, nrow(K.t)) # this is necessary to make K positive semidefinite
    K.t <- rbind(tmp.train$id, K.t)
    
    
    # apply GenomeBasedModel on the training set
    res <- GenomeBasedModel(Input = Input, Freevec = Freevec, Y = Y.train, 
                            Missing = 99999, Np = 2, Geno = Geno.train,
                            Methodcode = c(8, 8), Referencevalues = c(4, -6),
                            SdforParameters = c(4 * 0.02, 6 * 0.01),
                            SdforV = abs(mean(Y.train)) * 0.00001,
                            Modelfunction = asym.model,
                            Q = Q.train, K = K.t,
                            Nloop = 20)
    
    # prepare validation set
    
    data.valid.sorted <- val.data %>% arrange(id, time)
    tmp.valid <- data.valid.sorted %>% select(id, time, y) %>% spread(time, y)
    Y.valid <- t(tmp.valid)
    tmp.valid <- data.valid.sorted %>% select(id, variety, condD, time, y) %>% spread(time, y)
    
    var.C.valid <- tmp.valid$variety[!tmp.valid$condD]
    var.D.valid <- tmp.valid$variety[tmp.valid$condD]
    
    # compute one parameter value per variety and per condition in the training sample
    
    phi1 <- apply(res$Para[[1]], 2, mean)
    phi2 <- apply(res$Para[[2]], 2, mean)
    
    phi.var <- tibble(id=c(id.C.train,id.D.train),
                      variety=c(var.C.train,var.D.train),
                      condD = c(rep(0,length(id.C.train)),rep(1,length(id.D.train))),
                      phi1=phi1,
                      phi2=phi2) 
    
    
    phi1.var.mean <- phi.var %>%
      group_by(variety,condD) %>%
      summarise_at(vars(phi1), list(phi1.mean = mean))
    
    phi2.var.mean <- phi.var %>%
      group_by(variety,condD) %>%
      summarise_at(vars(phi2), list(phi2.mean = mean))
    
    phi1.C.train <- phi1.var.mean$phi1.mean[phi1.var.mean$condD==0]
    E.phi1.C <- mean(phi1.C.train)
    phi1.D.train <- phi1.var.mean$phi1.mean[phi1.var.mean$condD==1]
    E.phi1.D <- mean(phi1.D.train)
    
    phi2.C.train <- phi2.var.mean$phi2.mean[phi2.var.mean$condD==0]
    E.phi2.C <- mean(phi2.C.train)
    phi2.D.train <- phi2.var.mean$phi2.mean[phi2.var.mean$condD==1]
    E.phi2.D <- mean(phi2.D.train)
    
    # predict one parameter value per variety and per condition in the validation sample
    GRM.C.vt <- A[unique(var.C.valid),unique(var.C.train)]
    GRM.C.t <- A[unique(var.C.train),unique(var.C.train)]
    MM <- GRM.C.vt %*% solve(GRM.C.t)
    
    # parameter phi1 in condition C
    phi1.C.pred <- E.phi1.C +
      MM %*% (matrix(phi1.C.train,ncol=1) - matrix(E.phi1.C,ncol=1,nrow=length(phi1.C.train)))
    # parameter phi1 in condition D
    phi1.D.pred <- E.phi1.D +
      MM %*% (matrix(phi1.D.train,ncol=1) - matrix(E.phi1.D,ncol=1,nrow=length(phi1.D.train)))
    # parameter phi2 in condition C
    phi2.C.pred <- E.phi2.C +
      MM %*% (matrix(phi2.C.train,ncol=1) - matrix(E.phi2.C,ncol=1,nrow=length(phi2.C.train)))
    # parameter phi2 in condition D
    phi2.D.pred <- E.phi2.D +
      MM %*% (matrix(phi2.D.train,ncol=1) - matrix(E.phi2.D,ncol=1,nrow=length(phi2.D.train)))
    
    valid.pred[[g]] <- tibble(variety=unique(var.C.valid),phi1.C=c(phi1.C.pred),
                              phi2.C=c(phi2.C.pred), phi1.D=c(phi1.D.pred),
                              phi2.D=c(phi2.D.pred))
    
  }    
  
  save(valid.pred,file=paste(resultsrepositoryname,"/Soybean_data_2017_GBM_cv1_",repCV,".Rdata",sep=""))
}


