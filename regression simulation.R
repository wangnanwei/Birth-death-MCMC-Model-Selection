setwd("C:/Users/wangn/Dropbox/R/BDgraph snp interaction/code for paper")
library(parallel)
library(lattice)
library(ggplot2)
library(latex2exp)
library(expm)
library(SOIL)
library(Rlab)
library(glmnet)
library(foreach)
library(MASS)
library(doSNOW)
cl<-makeCluster(12) #change the 2 to your number of CPU cores
registerDoSNOW(cl)
#stopCluster(cl)

# the follow function is used to generate data for linear regression
datageneration<-function(n,p,rho,beta,sigma){
  S<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      S[i,j]<-rho^{abs(i-j)}
    }
  }
  X<-mvrnorm(n,rep(0,p),S) # random sample n*p covariate matrix
  #beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5))
  dm<-X%*%beta
  y<-NULL
  for(k in 1:n){
    y[k]<-rnorm(1,dm[k],sigma)
  }
  dd<-as.data.frame(cbind(X,y))
  return(dd)
}


# First sample simulation

n=1000;p=1000;rho=0.5; beta<-c(c(4,4,4,-6*sqrt(2),4/3),rnorm(5,0,1),rnorm(10,0,0.1), rep(0,p-20));sigma<-1
#beta<-c(rnorm(20,0,0.2),rep(0,p-20))
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
which(fit1$importance>0.5)
prob<-bdmcmc(X,y,1,50,50,1)
# compare BDMCMC method and SOIL method
which(fit1$importance>0.5)

which(prob>0.5)







## Binomial case 4


n=500;p=2000;rho=0;  beta<-c(4+rnorm(10),c(4,4,4,-6*sqrt(2),4/3), rep(0,p-15));sigma<-0.1
S<-matrix(0,p,p)
for (i in 1:p){
  for (j in 1:p){
    S[i,j]<-rho^{abs(i-j)}
  }
}
X<-mvrnorm(n,rep(0,p),S)
#beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5))
xb<-X%*%beta
pb<-exp(xb)/(1+exp(xb))
y<-NULL
for (i in 1:n){
  y[i]<-rbern(1, pb[i])
}
dd<-as.data.frame(cbind(X,y))



#SOIL method:
fit1<-SOIL(X,y,family="binomial",method="union",weight_type="BIC")
prob<-bdmcmcbin(X,y,1,30,20,1)
which(prob>0.5)
which(fit1$importance>0.5)





