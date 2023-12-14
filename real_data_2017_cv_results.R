##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################

####
## Analysis of cross-validation results on the 2017 soybean growth data under scenario
####

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

load(paste(datarepositoryname,"/FinalData2017.Rdata",sep=""))
nrep_CV  <- 4  # Nb of repetitions for cross-validation
nb_group <- 5  # Nb of samples for cross-validation

####
## Results from scenario cv1
####

# Comparion between predictions and observations per condition and time point

unique.time <- unique(data.log$time)

rsd.y.cv1.cond <- array(NA,dim=c(2,nrep_CV,nb_group,length(unique.time)))

# In condition C
for (repCV in 1:nrep_CV){
  load(paste(datarepositoryname,"/Soybean_data_2017_cv1_",repCV,".Rdata",sep=""))
  for (g in 1:nb_group){
    for (j in 1:length(unique.time)){
      tj <- which((res.CV[[4]][[g]]$time==unique.time[j]))
      ij <- which(res.CV[[4]][[g]]$condD[tj]==F)
      tj <- tj[ij]
      rsd.y.cv1.cond[1,repCV,g,j] <- mean(abs((res.CV[[4]][[g]]$y[tj]-res.CV[[4]][[g]]$ypred[tj])/res.CV[[4]][[g]]$y[tj]))
    }
  }  
}

# In condition D
for (repCV in 1:nrep_CV){
  load(paste(datarepositoryname,"/Soybean_data_2017_cv1_",repCV,".Rdata",sep=""))
  for (g in 1:nb_group){
    for (j in 1:length(unique.time)){
      tj <- which((res.CV[[4]][[g]]$time==unique.time[j]))
      ij <- which(res.CV[[4]][[g]]$condD[tj]==T)
      tj <- tj[ij]
      rsd.y.cv1.cond[2,repCV,g,j] <- mean(abs((res.CV[[4]][[g]]$y[tj]-res.CV[[4]][[g]]$ypred[tj])/res.CV[[4]][[g]]$y[tj]))
    }
  }  
}

apply(rsd.y.cv1.cond,c(1,4),mean)

ttime <- rep(unique.time,each=nb_group*nrep_CV)
ccond <- rep(c("C","D"),each=length(ttime))

rsd.y.cv1.cond <- apply(rsd.y.cv1.cond, c(1,4), c)
rsd.y.cv1.cond <- provideDimnames(rsd.y.cv1.cond, sep = "_", base = list('rep','cond','time'))
rsd.y.cv1.cond <- apply(rsd.y.cv1.cond,2,c)
tib.res.cv1 <- tibble(rsd=c(rsd.y.cv1.cond),cond=ccond,time=rep(ttime,2))
tib.res.cv1 <- mutate(tib.res.cv1, cond=as.factor(cond), time=as.factor(time))

cv1.plot <- ggplot(tib.res.cv1, aes(x=time, y=rsd*100, fill=cond)) +
  ylab("Rdiff (%)") +
  ylim(2.5,13) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90, size=12),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()



####
## Results from scenario cv2
####

# Comparion between predictions and observations per condition and time point

rsd.y.cv2.cond <- array(NA,dim=c(2,nrep_CV,nb_group,length(unique.time)))

# In condition C
for (repCV in 1:nrep_CV){
  load(paste(datarepositoryname,"/Soybean_data_2017_cv2_",repCV,".Rdata",sep=""))
  for (g in 1:nb_group){
    for (j in 1:length(unique.time)){
      tj <- which((res.CV[[4]][[g]]$time==unique.time[j]))
      ij <- which(res.CV[[4]][[g]]$condD[tj]==F)
      tj <- tj[ij]
      rsd.y.cv2.cond[1,repCV,g,j] <- mean(abs((res.CV[[4]][[g]]$y[tj]-res.CV[[4]][[g]]$ypred[tj])/res.CV[[4]][[g]]$y[tj]))
    }
  }  
}

# In condition D
for (repCV in 1:nrep_CV){
  load(paste(datarepositoryname,"/Soybean_data_2017_cv2_",repCV,".Rdata",sep=""))
  for (g in 1:nb_group){
    for (j in 1:length(unique.time)){
      tj <- which((res.CV[[4]][[g]]$time==unique.time[j]))
      ij <- which(res.CV[[4]][[g]]$condD[tj]==T)
      tj <- tj[ij]
      rsd.y.cv2.cond[2,repCV,g,j] <- mean(abs((res.CV[[4]][[g]]$y[tj]-res.CV[[4]][[g]]$ypred[tj])/res.CV[[4]][[g]]$y[tj]))
    }
  }  
}

apply(rsd.y.cv2.cond,c(1,4),mean)

ttime <- rep(unique.time,each=nb_group*nrep_CV)
ccond <- rep(c("C","D"),each=length(ttime))

rsd.y.cv2.cond <- apply(rsd.y.cv2.cond, c(1,4), c)
rsd.y.cv2.cond <- provideDimnames(rsd.y.cv2.cond, sep = "_", base = list('rep','cond','time'))
rsd.y.cv2.cond <- apply(rsd.y.cv2.cond,2,c)
tib.res.cv2 <- tibble(rsd=c(rsd.y.cv2.cond),cond=ccond,time=rep(ttime,2))
tib.res.cv2 <- mutate(tib.res.cv2, cond=as.factor(cond), time=as.factor(time))

cv2.plot <- ggplot(tib.res.cv2, aes(x=time, y=rsd*100, fill=cond)) +
  ylim(2.5,13.5) +
  ylab("Rdiff (%)") +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 90,size=12),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()


figure <- ggarrange(cv1.plot, cv2.plot,
                    labels = c("CV1", "CV2"),
                    ncol = 2, nrow = 1)
figure

ggsave(paste(resultsrepositoryname,"/cv-real-data.png",sep=""))