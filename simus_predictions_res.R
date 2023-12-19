##################################################
## Modeling soybean growth: a mixed model approach.
## Delattre, M. et al.
##################################################

####
## Analysis of prediction results on simulated data
####

rm(list = ls())

## Source main functions and load necessary R libraries
library(tidyverse)
library(Matrix)
library(mvnfast)
library(ggpubr)

resultsrepositoryname <- dirname(rstudioapi::getActiveDocumentContext()$path) ## Specify the name of the repository in which storing the results

time  <- seq(1,50,length.out=10) # time points
nbrep <- 50

### First heritability value
######

## Prediction of growth parameters

load(paste(resultsrepositoryname,"/Simu_pred1.Rdata",sep=""))

summary(res.pred$rsd.A.cv0)
summary(res.pred$rsd.A.cv1)
summary(res.pred$rsd.A.cv2)

summary(res.pred$rsd.x.cv0)
summary(res.pred$rsd.x.cv1)
summary(res.pred$rsd.x.cv2)

res.pred1 <- tibble(rsdA = c(res.pred$rsd.A.cv0,res.pred$rsd.A.cv1,res.pred$rsd.A.cv2),
                    rsdx = c(res.pred$rsd.x.cv0,res.pred$rsd.x.cv1,res.pred$rsd.x.cv2),
                    cv   = rep(c("cv0","cv1","cv2"),each=nbrep)) 


plot.A.1 <- ggplot(res.pred1, aes(x=cv, y=rsdA*100)) +
  ylab("Rdiff (%)") +
  xlab("") +
  #ylim(2.5,13) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()

plot.x.1 <- ggplot(res.pred1, aes(x=cv, y=rsdx*100)) +
  ylab("Rdiff (%)") +
  xlab("") +
  #ylim(2.5,13) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()

## Prediction of height measures

res.pred1.y <- tibble(y=c(c(res.pred$rsd.y.cv0),c(res.pred$rsd.y.cv1),c(res.pred$rsd.y.cv2)),
                      t=ceiling(rep(rep(time,each=nbrep),3)),
                      cv   = rep(c("cv0","cv1","cv2"),each=nbrep*10))
res.pred1.y <- mutate(res.pred1.y,t=as.factor(t))

plot.y.1 <- ggplot(res.pred1.y, aes(x=t, y=log(y), fill=cv)) +
  ylab("log-Rdiff (%)") +
  xlab("") +
  ylim(-3,3) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()


## Second heritability value

load(paste(resultsrepositoryname,"/Simu_pred2.Rdata",sep=""))

summary(res.pred$rsd.A.cv0)
summary(res.pred$rsd.A.cv1)
summary(res.pred$rsd.A.cv2)

summary(res.pred$rsd.x.cv0)
summary(res.pred$rsd.x.cv1)
summary(res.pred$rsd.x.cv2)

res.pred2 <- tibble(rsdA = c(res.pred$rsd.A.cv0,res.pred$rsd.A.cv1,res.pred$rsd.A.cv2),
                    rsdx = c(res.pred$rsd.x.cv0,res.pred$rsd.x.cv1,res.pred$rsd.x.cv2),
                    cv   = rep(c("cv0","cv1","cv2"),each=nbrep)) 


plot.A.2 <- ggplot(res.pred2, aes(x=cv, y=rsdA*100)) +
  ylab("Rdiff (%)") +
  xlab("") +
  #ylim(2.5,13) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()

plot.x.2 <- ggplot(res.pred2, aes(x=cv, y=rsdx*100)) +
  ylab("Rdiff (%)") +
  xlab("") +
  #ylim(2.5,13) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()

figure <- ggarrange(plot.A.1, plot.x.1,
                    plot.A.2, plot.x.2,
                    labels = c("A-1", "x-1", "A-2", "x-2"),
                    ncol = 2, nrow = 2)
figure

ggsave(paste(resultsrepositoryname,"/cv-simus-Ax.png",sep=""))

## Prediction of height measures

res.pred2.y <- tibble(y=c(c(res.pred$rsd.y.cv0),c(res.pred$rsd.y.cv1),c(res.pred$rsd.y.cv2)),
                      t=ceiling(rep(rep(time,each=nbrep),3)),
                      cv   = rep(c("cv0","cv1","cv2"),each=nbrep*10))
res.pred2.y <- mutate(res.pred2.y,t=as.factor(t))

plot.y.2 <- ggplot(res.pred2.y, aes(x=t, y=log(y), fill=cv)) +
  ylab("log-Rdiff (%)") +
  xlab("") +
  ylim(-3,3) +  
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = 45, size=12, face="bold"),
        axis.text.y=element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_boxplot()



figure <- ggarrange(plot.y.1, plot.y.2,
                    labels = c("1", "2"),
                    ncol = 2, nrow = 1)
figure

ggsave(paste(resultsrepositoryname,"/cv-simus-y.png",sep=""))
