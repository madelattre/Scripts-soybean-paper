##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################


## Script for the analysis of the 2017 soybean growth data

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

## Load GRM matrix

load(paste(datarepositoryname,"/GRM2017.Rdata",sep=""))

## Original data
datafile <- paste(datarepositoryname,
                  "/2017_Tottori_Jul_PlantHeight_Imputed.csv",sep="")


UAV_height  <- read_csv(datafile,col_names=TRUE) %>% 
  rowid_to_column(var = "LineID")%>%
  gather(key = "Day", value = "Height", 11:20) %>%
  mutate(Date= as.Date(substring(Day, first=2, last= nchar(Day)), format = "%y%m%d")) %>%
  mutate(condition=substring(block, first=1, last= 1)) %>%
  filter(Height != "NA") %>%
  filter(!is.na(varietyID)) %>%
  filter(variety %in% colnames(A)) %>% 
  filter(Height>=0)

## Time in heat unit
heattimefile <- paste(datarepositoryname,
                      "/SoyCREST_Weather_Tottori_2017.csv",sep="")
heat.time <- read_csv(heattimefile,col_names = TRUE) %>%
  select(Date,HeatUnit)

## Processed data
UAV_final <- merge(UAV_height,heat.time,id="Date")

data <- UAV_final %>% 
  rename(c("y"="Height","id"="LineID","varID"="varietyID",
           "time"="HeatUnit")) %>% 
  mutate(condD = condition=="D")


Nv    <- ncol(A) ## Number of varieties
l_id  <- unique(data$id) ## Sequence of plant identifiers
Ntot  <- nrow(data) ## Total number of observations
N     <- length(l_id) ## Number of plants
l_var <- select(data, varID, variety) %>% distinct() 
var_A <- tibble(order = 1:Nv,variety = colnames(A))
l_var <- left_join(l_var, var_A, by="variety") %>% arrange(order) ## Identifiers and names of the varieties

data <- data %>% select(id, y, varID, time, condD, condition, variety) 

## Convert plant height into their logarithm
data.log <- data %>% mutate(y = log(y))


## Graphical representation of the growth data as a function of heating time
## change data into data.log to represent the logarithm of height as a function
## of heating time
plot_variety <- ggplot(data,
                       aes(x=time, y=y, group=id,
                           color=as.factor(condition))) +
  geom_line() + theme_bw() + facet_grid(~condition) + 
  theme(legend.position = "none",
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 45)) + 
  ylab("Height") + xlab("Heat time")

plot_variety

# plot_variety <- ggplot(data.log,
#                        aes(x=time, y=y, group=id,
#                            color=as.factor(condition))) +
#   geom_line() + theme_bw() + facet_grid(~condition) + 
#   theme(legend.position = "none",
#         axis.title.y = element_text(size=16),
#         axis.title.x = element_text(size=16),
#         axis.text.y = element_text(size=12),
#         axis.text.x = element_text(size=12, angle = 45)) + 
#   ylab("Log Height") + xlab("Heat time")
# 
# plot_variety


save(data.log,file=paste(datarepositoryname,"/FinalData2017.Rdata",sep=""))

####
## 2- Analysis of the logarithms of the heights using the asymptotic growth model 
####

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
  condi <- 1-as.numeric(data$condD[data$id==l_id[i]][1])
  Xi[i,1:nb.phi,1:nb.phi] <- diag(nb.phi)
  Xi[i,1:nb.phi,(nb.phi+1):(2*nb.phi)] <- diag(condi,nb.phi)
  vari <- data$varID[data$id==l_id[i]][1]
  if (condi==0){
    Zi[i,1:nb.phi,1:(nb.phi*Nv)] <- diag(nb.phi)%x%t(l_var$varID==vari)
  } else{
    Zi[i,1:nb.phi,(nb.phi*Nv+1):(d*Nv)] <- diag(nb.phi)%x%t(l_var$varID==vari)
  }
}


## Initialization and algorithmic settings 

paraminit.log <- list(beta = matrix(c(3,-5,1,0),2*nb.phi,1),G=diag(rep(5,d)),
                      P=diag(c(5,5)) , sigma2= 10)
niterSAEM <- 400
nburninSAEM <- 300

## Estimation of the population parameters


res.theta <- SAEM_GRM(niter=niterSAEM, nburnin=nburninSAEM, data=data.log, 
                      predictors=c(4,5), paraminit=paraminit.log, 
                      GRM_mat = A, Xi=Xi, Zi=Zi, model = asym.model,
                      Nchain=10)

                
param.est <- list(beta=res.theta$beta[niterSAEM+1,],
                  sigma2=res.theta$sigma2[niterSAEM+1],
                  P=res.theta$P[niterSAEM+1,,],
                  G=res.theta$G[niterSAEM+1,,])



## Estimation of the predicted genetic values

niterPred <- 1000
nburninPred <- 800

res.pred <- predict.u(niter=niterPred, nburnin=nburninPred, data=data.log,
                      predictors=c(4,5), param.est, GRM_mat=A, Xi, Zi,
                      model=asym.model, Nsim=10, Nchain = 10)

## Saving results

res.soybean.2017 <- list(res.pred=res.pred, res.theta=res.theta)
res.soybean.2017.file <- 'Res.soybean.2017.Rdata'
save(res.soybean.2017,file=paste(resultsrepositoryname,'/',res.soybean.2017.file,sep=''))

load(paste(resultsrepositoryname,'/',res.soybean.2017.file,sep=''))

####
## 3- Convergence graphs for the parameter estimations
####

nbeta <- dim(res.soybean.2017$res.theta$beta)[2]
niter <- dim(res.soybean.2017$res.theta$beta)[1]

convSAEM <- tibble(beta1 = res.soybean.2017$res.theta$beta[,1],
                   beta2 = res.soybean.2017$res.theta$beta[,2],
                   beta3 = res.soybean.2017$res.theta$beta[,3],
                   beta4 = res.soybean.2017$res.theta$beta[,4],
                   sigma2 = res.soybean.2017$res.theta$sigma2,
                   P11 = res.soybean.2017$res.theta$P[,1,1],
                   P22 = res.soybean.2017$res.theta$P[,2,2],
                   P12 = res.soybean.2017$res.theta$P[,1,2],
                   G11 = res.soybean.2017$res.theta$G[,1,1],
                   G22 = res.soybean.2017$res.theta$G[,2,2],
                   G33 = res.soybean.2017$res.theta$G[,3,3],
                   G44 = res.soybean.2017$res.theta$G[,4,4],
                   G12 = res.soybean.2017$res.theta$G[,1,2],
                   G13 = res.soybean.2017$res.theta$G[,1,3],
                   G14 = res.soybean.2017$res.theta$G[,1,4],
                   G23 = res.soybean.2017$res.theta$G[,2,3],
                   G24 = res.soybean.2017$res.theta$G[,2,4],
                   G34 = res.soybean.2017$res.theta$G[,3,4]
)

plotbeta1 <- ggplot(convSAEM) + geom_line(aes(y = beta1, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(mu[A]),x="Iteration")
plotbeta2 <- ggplot(convSAEM) + geom_line(aes(y = beta2, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(mu[B]),x="Iteration")
plotbeta3 <- ggplot(convSAEM) + geom_line(aes(y = beta3, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(delta[A]),x="Iteration")
plotbeta4 <- ggplot(convSAEM) + geom_line(aes(y = beta4, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(delta[B]),x="Iteration")

plotsigma2 <- ggplot(convSAEM) + geom_line(aes(y = sigma2, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(sigma^2),x="Iteration") 

plotP11 <- ggplot(convSAEM) + geom_line(aes(y = P11, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(P[A]^2),x="Iteration") + ylim(0,0.5)
plotP22 <- ggplot(convSAEM) + geom_line(aes(y = P22, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(P[B]^2),x="Iteration") + ylim(0,0.5)
plotP12 <- ggplot(convSAEM) + geom_line(aes(y = P12, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(P[A-B]),x="Iteration")

plotG11 <- ggplot(convSAEM) + geom_line(aes(y = G11, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[mu[A]]^2),x="Iteration") + ylim(0,0.5)
plotG22 <- ggplot(convSAEM) + geom_line(aes(y = G22, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[mu[B]]^2),x="Iteration") + ylim(0,0.5)
plotG33 <- ggplot(convSAEM) + geom_line(aes(y = G33, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[delta[A]]^2),x="Iteration")+ ylim(0,0.5)
plotG44 <- ggplot(convSAEM) + geom_line(aes(y = G44, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[delta[B]]^2),x="Iteration")+ ylim(0,0.5)

plotG12 <- ggplot(convSAEM) + geom_line(aes(y = G12, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(mu[A]~","~mu[B])]),x="Iteration")
plotG13 <- ggplot(convSAEM) + geom_line(aes(y = G13, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(mu[A]~","~delta[A])]),x="Iteration")
plotG14 <- ggplot(convSAEM) + geom_line(aes(y = G14, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(mu[A]~","~delta[B])]),x="Iteration")
plotG23 <- ggplot(convSAEM) + geom_line(aes(y = G23, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(mu[B]~","~delta[A])]),x="Iteration")
plotG24 <- ggplot(convSAEM) + geom_line(aes(y = G24, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(mu[B]~","~delta[B])]),x="Iteration")
plotG34 <- ggplot(convSAEM) + geom_line(aes(y = G34, x = seq(from=0,by=1,length.out = niter))) + labs(y=expression(G[(delta[A]~","~delta[B])]),x="Iteration")

plot_grid(plotbeta1,plotbeta2,plotbeta3,plotbeta4,plotsigma2, ncol = 3, nrow = 2)
plot_grid(plotP11, plotP22, plotP12, ncol = 2, nrow = 2)
plot_grid(plotG11, plotG22, plotG33, plotG44, ncol = 2, nrow = 2)
plot_grid(plotG12, plotG13, plotG14, plotG23, plotG24, plotG34, ncol = 3, nrow = 2)

####
## 4- Estimated heritability values
####

h.A <- diag(param.est$G)[1]/(diag(param.est$G)[1]+diag(param.est$P)[1])
h.dA <- diag(param.est$G)[3]/(diag(param.est$G)[3]+diag(param.est$P)[1])

h.B <- diag(param.est$G)[2]/(diag(param.est$G)[2]+diag(param.est$P)[2])
h.dB <- diag(param.est$G)[4]/(diag(param.est$G)[4]+diag(param.est$P)[2])

# ####
# ## 5- Convergence graphs for the prediction of the genetic values
# ####

par(mfrow=c(2,3))
for (j in 1:10){
  plot(res.soybean.2017$res.pred$u[j,],type='l',main=j,xlab='Iteration',ylab=paste('u[',j,']'))
}

