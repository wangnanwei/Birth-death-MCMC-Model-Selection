library(igraph)
library(BDgraph)
library(network)
library(ggplot2)
library(sna)
library(GGally)
library(MASS)
library(glmnet)
library(doSNOW)
library(qqman)
library(foreach)
library(mgm)
library(parallel)
library(lattice)
library(latex2exp)
library(expm)
library(Rlab)
library(glmnet)
library(vcd)
library(vcdExtra)


source("bdmcmc algorithm_v6.R")

nodesize=50
g1=barabasi.game(nodesize,directed=FALSE)
adj=g1[1:nodesize]
adj=as.matrix(adj)
adj2=adj
for(i in 1:(50-1)){
  for(j in (i+1):50) adj2[j,i]<-0
}
trueedge=which(adj2>0)


temp<-c(rep("Gaussian",30),rep("Binary",15),rep("Count",5))  
net1<-network(adj,directed=FALSE)
net1 %v% "type"<-temp
col = c("Gaussian" = "green", "Binary" = "skyblue","Count"="tomato")
ggnet2(net1,node.size=7,color="type",palette = col,label = TRUE,label.size = 5)


## Gibbs sampling

x0<-rep(0,50)
N=10^5
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

dd=dd[(1+10000):10^5,]







nodesize=50

result=matrix(0,5,3)
samplesize=c(200,500,1000,2000,3000)
for(ss in 1:5){
  
  F1mgm=F1bdmcmc=F1sbdmcmc=rep(0,20)
  for(iter in 1:20){
    
    print(iter)
    k=iter
    dd1=dd[(1+(k-1)*samplesize[ss]):(k*samplesize[ss]),]
    
    
    mgmfit=mgm(data=dd1,type=c(rep("g",30),rep("c",15),rep("p",5)),level=c(rep(1,30),rep(2,15),rep(1,5)),lambda.sel="EBIC",k=2)
    edgemgm=mgmfit$interactions$indicator[[1]]
    adjmgm=matrix(0,nodesize,nodesize)
    for(i in 1:dim(edgemgm)[1]){
      adjmgm[edgemgm[i,1],edgemgm[i,2]]=1
    }
    mgmedge=which(adjmgm>0)
    TP=sum(!is.na(match(mgmedge,trueedge)))
    FP=sum(is.na(match(mgmedge,trueedge)))
    FN=sum(is.na(match(trueedge,mgmedge)))
    F1mgm[iter]=2*TP/(2*TP+FP+FN)
    NS<-bdgraph(data=dd1,method="gcgm",algorithm="bdmcmc",iter=2000,burnin=1000)
    #summary(NS)
    adjbdmcmc<-summary(NS)$selected_g
    bdmcmcedge=which(adjbdmcmc>0)
    TP=sum(!is.na(match(bdmcmcedge,trueedge)))
    FP=sum(is.na(match(bdmcmcedge,trueedge)))
    FN=sum(is.na(match(trueedge,bdmcmcedge)))
    F1bdmcmc[iter]=2*TP/(2*TP+FP+FN)
    
    cl<-makeCluster(12) #change the 2 to your number of CPU cores
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = nodesize, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    pedge<-matrix(0,nodesize,nodesize)
    prob=foreach(i = 1:nodesize, .combine = cbind, .options.snow = opts) %dopar% {
      
      y=dd1[,i]
      X=dd1[,-i]
      if (i<31){
        prob=bdmcmc(X,y,1,20,20,0)
      } else if (i<46) {
        prob<-bdmcmcbin(X,y,1,20,20,0)
      } else{
        prob<-bdmcmcpos(X,y,1,20,20,0)
      }
      
    }
    stopCluster(cl)
    for(i in 1:nodesize){
      pedge[i,-i]<-prob[,i]
      
    }
    
    pedge[which(pedge>=0.5)]=1
    pedge[which(pedge<0.5)]=0
    pedge=pedge+t(pedge)
    pedge[which(pedge>0)]=1
    for(i in 1:(nodesize-1)){
      for(j in (i+1):nodesize) pedge[j,i]<-0
    }
    sbdmcmcedge=which(pedge>0)
    TP=sum(!is.na(match(sbdmcmcedge,trueedge)))
    FP=sum(is.na(match(sbdmcmcedge,trueedge)))
    FN=sum(is.na(match(trueedge,sbdmcmcedge)))
    F1sbdmcmc[iter]=2*TP/(2*TP+FP+FN)
    
    
  }
  result[ss,]=rowMeans(rbind(F1mgm,F1bdmcmc,F1sbdmcmc))
}



df <- data.frame("sample size"=rep(samplesize, 3), "F1 score"=c(result[,1],result[,2],result[,3]),variable=c(rep("MGM",5),rep("BDMCMC",5),rep("SBDMCMC",5)))
ggplot(df,aes(x=sample.size, y=F1.score)) + geom_line(aes(colour=variable))