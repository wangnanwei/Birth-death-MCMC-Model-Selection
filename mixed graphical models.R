## local mixed graphical models simulation:
library(parallel)
library(lattice)
library(ggplot2)
library(latex2exp)
library(expm)
library(SOIL)
library("GGally")
library("network")
library(sna)
library(ggplot2)
library(Rlab)
library(glmnet)
library("vcd")
library("vcdExtra")
library("foreach")
library(doMC)
registerDoMC(8) 
adj<-matrix(0,50,50)

E<-combn(50,2)

k=0.05
index<-sample(1:dim(E)[2],dim(E)[2]*k)

edges<-E[,index]
for(i in 1:dim(edges)[2]){
  adj[edges[1,i],edges[2,i]]<-1
  adj[edges[2,i],edges[1,i]]<-1
  
}

temp<-c(rep("Gaussian",30),rep("Binary",15),rep("Count",5))  


net1<-network(adj,directed=FALSE)
net1 %v% "type"<-temp
col = c("Gaussian" = "steelblue", "Binary" = "gold","Count"="tomato")
ggnet2(net1,node.size=5,color="type",palette = col,label = TRUE)


## Gibbs sampling

x0<-rep(0,50)
N=10000
dd<-matrix(0,N,50)
beta0<-rnorm(50,0,0.1)
         
#beta0[:100]<--abs(beta0[81:100])
beta<-matrix(rnorm(50^2,0,1),50,50)
beta=(beta+t(beta))/10
beta<-adj*beta

for(i in 46:50){
  beta[i,]=-abs(beta[i,])
  beta[,i]=-abs(beta[,i])
}
x=x0
for(n in 1:N){
  print(n)
  x1=x
  for(i in 1:30){
    index<-which(beta[i,]!=0)
    mu=beta0[i]+ beta[i,index]%*% x[index]
    
    x[i]<-rnorm(1,mu,1)
  }
  for(i in 31:45){
    index<-which(beta[i,]!=0)
    mu=beta0[i]+ beta[i,index]%*% x[index]
    mu=exp(mu)/(1+exp(mu))
    x[i]<-rbern(1,mu)
  }
  for(i in 46:50){
    index<-which(beta[i,]!=0)
    mu=beta0[i]+ beta[i,index]%*% x[index]
    mu=exp(mu)
    x[i]<-rpois(1,mu)
  }
  dd[n,]=x
  #if (max(x)<20){
   # dd[n,]=x 
  #} else{
   # x=x1
  #}
}

dd1=dd[5001:6000,]

y=dd1[,5]
X=dd1[,-5]
fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
prob<-bdmcmc(X,y,1,50,50,1)
which(fit1$importance>0.5)
which(prob>0.5)
prob<-bdmcmc(X,y)





## Birth-death MCMC methods
bdmcmc<-function(X,y){
  dd<-as.data.frame(cbind(X,y))
  p<-dim(X)[2]
  nei<-matrix(0,501,p)
  w<-matrix(0,500,1)
  temp<-nei[1,]
  index<-which(temp>0)
  if (length(index)==0){
    fit<-lm(y~1,data=dd)  
  } else{
    fit<-lm(y~0+.,data=dd[,c(index,p+1)])  
  }
  pg0<-BIC(fit)
  len0<-length(index)
  
  for (k in 1:500){
    print(k)
    pg<-matrix(0,1,p)
    len<-rep(0,p)
    
    temp2<-foreach (i=1:p, .combine = cbind) %dopar%{
      library(expm)
      temp<-nei[k,]
      temp[i]<-abs(temp[i]-1)
      index<-which(temp>0)
      if (length(index)==0){
        fit<-lm(y~1,data=dd)  
      } else{
        fit<-lm(y~0+.,data=dd[,c(index,p+1)])  
      }
      pg[i]<-BIC(fit)
      len[i]<-length(index)
      rbind(pg[i],len[i])
    }
    pg<-temp2[1,]
    if (max(pg)>5000) {
      pg=pg/100
      pg0=pg0/100
    }
    len<-temp2[2,]
    pg1<-exp(0.5*(pg0-pg))
    prior<-0.2^temp2[2,]/0.2^len0
    
    pg1<-pg1*prior
    #pg1[pg1>1]<-1
    w[k]<-1/sum(pg1)
    pg1<-pg1/sum(pg1)
    sn<-sample(x=c(1:p),size=1,prob=pg1)
    nei[k+1,]<-nei[k,]
    nei[k+1,sn]<-abs(nei[k+1,sn]-1)
    pg0<-pg[sn]
    len0<-length(which(nei[k+1,]>0))
    print(pg0)
    print(which(nei[k+1,]>0))
    
  }
  
  result<-cbind(nei[101:500,] ,w[101:500])
  
  result[,p+1]<-result[,p+1]/sum(result[,p+1])
  prob<-NULL
  for(i in 1:p){
    index<-which(result[,i]>0)
    prob[i]<-sum(result[index,p+1])
  }
  
  return(prob)
}





pedge2<-matrix(0,50,50)
for(i in 1:50){
  y=dd1[,i]
  X=dd1[,-i]
  if (i<31){
    cv.out<-cv.glmnet(X,y,alpha=1,family="gaussian",type.measure ="mse")
    index<-coef(cv.out,s=cv.out$lambda.1se)[2:50]
  } else if (i<46) {
    cv.out<-cv.glmnet(X,y,alpha=1,family="binomial",type.measure ="class")
    index<-coef(cv.out,s=cv.out$lambda.1se)[2:50]
  } else{
    cv.out<-cv.glmnet(X,y,alpha=1,family="poisson",type.measure ="deviance")
    index<-coef(cv.out,s=cv.out$lambda.1se)[2:50]
  }
  pedge2[i,-i]<-index
}

pedge2[pedge2!=0]<-1

b<-which(pedge2==1)
a<-which(adj==1)

allt<-length(a)

allf<-dim(E)[2]-allt
TP2<-length(intersect(a,b))/allt
FP2<- length(setdiff(b,a))/allf

















pedges1<-matrix(0,50,50)
imp<-foreach(i=1:50, .combine = "rbind") %dopar% {
  library(SOIL)
  y=dd[,i]
  X=dd[,-i]
  fit1<-SOIL(X,y,family="gaussian",method="union",weight_type="BIC")
  fit1$importance
}
