setwd("C:/Users/postdoc/Dropbox/R/BDgraph snp interaction/code for paper")
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
datageneration<-function(n,p,rho,beta,sigma){
  S<-matrix(0,p,p)
  for (i in 1:p){
    for (j in 1:p){
      S[i,j]<-rho^{abs(i-j)}
    }
  }
  X<-mvrnorm(n,rep(0,p),S)
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

n=200;p=1000;rho=0; beta<-c(c(4,4,4,-6*sqrt(2),4/3),rnorm(5,0,1),rnorm(10,0,0.1), rep(0,p-20));sigma<-0.1
#beta<-c(rnorm(20,0,0.2),rep(0,p-20))
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
which(fit1$importance>0.5)
prob<-bdmcmc(X,y,1,50,50,1)
  which(prob>0.5)

x<-c(1:50)
y1<-fit1$importance[1:50]
y2<-prob[1:50]

com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =1, \\sigma^2 =0.01$")
)+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))





#  sample simulation2
n=100;p=1000;rho=0.9; beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5));sigma<-0.1
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:50)
y1<-fit1$importance[1:50]
y2<-prob[1:50]


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =1, \\sigma^2 =0.01$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))










#  sample simulation3
n=100;p=1000;rho=0; beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5));sigma<-5
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:50)
y1<-fit1$importance[1:50]
y2<-prob[1:50]


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =1, \\sigma^2 =0.01$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))





#  sample simulation4
n=100;p=1000;rho=0.9; beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5));sigma<-5
dd<-datageneration(n,p,rho,beta,sigma)
write.csv(dd,file="data_soil_bad.csv")
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:1000)
y1<-fit1$importance[1:1000]
y2<-prob[1:1000]


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =1, \\sigma^2 =0.01$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))












## example 2
n=150;p=14;rho=0.9; beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5));sigma<-5

S<-matrix(0,p,p)
for (i in 1:p){
  for (j in 1:p){
    S[i,j]<-rho^{abs(i-j)}
  }
}
X<-mvrnorm(n,rep(0,p),S)
X1<-0.5*X[,1]+2*X[,4]+mvrnorm(1,rep(0,n),diag(0.01,n,n))
X<-cbind(X,X1)
beta<-c(beta,0)
#beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5))
dm<-X%*%beta
y<-NULL
for(k in 1:n){
  y[k]<-rnorm(1,dm[k],sigma)
}
dd<-as.data.frame(cbind(X,y))


fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:15)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0, \\sigma^2 =0.01$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))






## example 3

n=150;p=8;rho=0.9; beta<-rep(0,p);sigma<-5
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:8)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =25$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))






## example 4

n=150;p=8;rho=0.9; beta<-rep(1,p);sigma<-5
dd<-datageneration(n,p,rho,beta,sigma)
X<-dd[,1:p];y<-dd[,p+1]
#SOIL method:
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:8)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =25$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))



















## Binomial case 1


n=80;p=7;rho=0.9; beta<-c(1,1/2,1/3,1/4,1/5,1/6,0);sigma<-5
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
prob<-bdmcmcbin(X,y,1,50,50,1)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:7)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =25$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))+ylim(0,1)








## Binomial case 2


n=5000;p=7;rho=0.9; beta<-c(1,1/2,1/3,1/4,1/5,1/6,0);sigma<-0.1
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
prob<-bdmcmcbin(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:7)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =0.01$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))+ylim(0,1)






## Binomial case 3


n=150;p=20;rho=0.9;  beta<-c(c(4,4,4,-6*sqrt(2),4/3), rep(0,p-5));sigma<-5
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
prob<-bdmcmcbin(X,y)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:20)
y1<-as.vector(fit1$importance)
y2<-prob


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =25$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))+ylim(0,1)





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
prob<-bdmcmcbin(X,y,0.5,30,20,0)
which(prob>0.5)
which(fit1$importance>0.5)
x<-c(1:50)
y1<-as.vector(fit1$importance)[1:50]
y2<-prob[1:50]


com<-as.data.frame(cbind(x,y1,y2))
theme_update(plot.title = element_text(hjust = 0.5))
ggplot()+geom_point(data=com,aes(x,y1,colour="SOIL"))+
  geom_point(data=com,aes(x,y2,colour="BDMCMC"))+
  labs(y="Importance", x="Variable index", title=TeX("$ \\rho =0.9, \\sigma^2 =25$")
  )+scale_colour_manual("compare two method",values = c("SOIL"="red", "BDMCMC"="blue"))+ylim(0,1)




## exmaple form bdmcmc
X<-dd[,-2]
y<-dd[,2]
fit1<-SOIL(X,y,family="binomial",method="union",weight_type="BIC")
































































# compute the evidence of the baysian method:
evidence<-function(p,n,X,y){
  a<-0.01
  b<-0.01
  
  x<-as.matrix(X)
  m<-diag(1,p,p)+t(x)%*%x
  n2<-solve(m)
  u1<-n2%*%(t(x)%*%y)
  a1<-a+n/2
  b1<-b+0.5*(t(y)%*%y-t(u1)%*%u1)
  evi<-((2*pi)^(n/2)*sqrt(det(m))*gamma(a1)*(b^a))/(b1^a1*gamma(a))
  return(evi)
}

