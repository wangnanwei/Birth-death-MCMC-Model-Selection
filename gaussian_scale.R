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

source("bdmcmc algorithm_v6.R")




nodesize=100
g1=barabasi.game(nodesize,directed=FALSE)
adj=g1[1:nodesize]
adj=as.matrix(adj)
net1<-network(adj,directed=FALSE)
col = c( "binary" = "gold")
ggnet2(net1,node.size=6,color="green",label = TRUE,label.size = 5)
K=matrix(0,nodesize,nodesize)
K=rgwish(n=1,adj=adj,b=3)
S=solve(K)

adj2=adj
for(i in 1:(nodesize-1)){
  for(j in (i+1):nodesize) adj2[j,i]<-0
}
trueedge=which(adj2>0)


result=matrix(0,5,3)
samplesize=c(200,500,1000,2000,3000)

for(ss in 1:5){
  
  F1mgm=F1bdmcmc=F1sbdmcmc=rep(0,50)
  
  for(iter in 1:50){
    print(iter)
    dd=mvrnorm(samplesize[ss],rep(0,nodesize),S)
    #dd=rmvggm(n.samples = samplesize[ss]-100,net.str = adj,cor=0.1)$dat
    
    mgmfit=mgm(data=dd,type=rep("g",nodesize),lambda.sel="EBIC",k=2)
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
    F1mgm[iter]
    NS<-bdgraph.mpl(data=dd,method="ggm",algorithm="bdmcmc",iter=5000,burnin=2000)
    #summary(NS)
    adjbdmcmc<-summary(NS)$selected_g
    bdmcmcedge=which(adjbdmcmc>0)
    TP=sum(!is.na(match(bdmcmcedge,trueedge)))
    FP=sum(is.na(match(bdmcmcedge,trueedge)))
    FN=sum(is.na(match(trueedge,bdmcmcedge)))
    F1bdmcmc[iter]=2*TP/(2*TP+FP+FN)
    F1bdmcmc[iter]
    
    cl<-makeCluster(4) #change the 2 to your number of CPU cores
    registerDoSNOW(cl)
    
    pedge<-matrix(0,nodesize,nodesize)
    prob=foreach(i = 1:nodesize, .combine = cbind) %dopar% {
      
      y=dd[,i]
      X=dd[,-i]
      prob<-bdmcmc(X,y,1,50,50,0)
    }
    stopCluster(cl)
    for(i in 1:nodesize){
      pedge[i,-i]<-prob[,i]
      
    }
    
    pedge[which(pedge>=0.5)]=1
    pedge[which(pedge<0.5)]=0
    
    for(i in 1:(nodesize-1)){
      for(j in (i+1):nodesize) pedge[j,i]<-0
    }
    sbdmcmcedge=which(pedge>0)
    TP=sum(!is.na(match(sbdmcmcedge,trueedge)))
    FP=sum(is.na(match(sbdmcmcedge,trueedge)))
    FN=sum(is.na(match(trueedge,sbdmcmcedge)))
    F1sbdmcmc[iter]=2*TP/(2*TP+FP+FN)
    F1sbdmcmc[iter]
    
  }
  
  rbind(F1mgm,F1bdmcmc,F1sbdmcmc)
  result[ss,]=rowMeans(rbind(F1mgm,F1bdmcmc,F1sbdmcmc))
  
}
#save.image(file="g.RData")
df <- data.frame("sample size"=rep(samplesize, 3), "F1 score"=c(result[,1],result[,2],result[,3]),Algorithm=c(rep("MGM",5),rep("BDMCMC",5),rep("SBDMCMC",5)))
ggplot(df,aes(x=sample.size, y=F1.score)) + geom_line(aes(colour=Algorithm),size=1)+theme(legend.position = c(.8, .3),legend.text = element_text(size = 18))



