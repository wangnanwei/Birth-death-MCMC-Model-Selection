obs<-bdgraph.sim(p=16,graph="random",n=1000,type="binary",vis=TRUE)
g0<-obs$G
dd<-obs$data


temp<-c(rep("binary",16))  


net1<-network(g0,directed=FALSE)
net1 %v% "type"<-temp
col = c( "binary" = "gold")
ggnet2(net1,node.size=5,color="gold",label = TRUE)


NS<-bdgraph.mpl(data=dd,method="dgm-binary",algorithm="bdmcmc",iter=5000,burnin=2000,save.all=TRUE)
summary(NS)
g1<-summary(NS)$selected_g


a<-which(g0==1)

allt<-length(a)
E<-combn(16,2)
allf<-dim(E)[2]-allt
b1<-which(g1==1)


TP<-length(intersect(a,b1))/allt
FP<- length(setdiff(b1,a))/allf


pedge<-matrix(0,16,16)
for(i in 1:16){
  print(i)
  y=dd[,i]
  X=dd[,-i]
  prob<-bdmcmcbin(X,y)
  pedge[i,-i]<-prob
}



pe<-(t(pedge)+pedge)/2
pe[pe>=0.5]<-1
pe[pe<0.5]<-0


allt<-length(a)

allf<-dim(E)[2]-allt
b1<-which(pe==1)


TP<-length(intersect(a,b1))/allt
FP<- length(setdiff(b1,a))/allf




pedge2<-matrix(0,16,16)
for(i in 1:16){
  y=dd[,i]
  X=dd[,-i]
    cv.out<-cv.glmnet(X,y,alpha=1,family="binomial",type.measure ="class")
    index<-coef(cv.out,s=cv.out$lambda.1se)[2:16]
    
  pedge2[i,-i]<-index
}

pedge2[pedge2!=0]<-1

b<-which(pedge2==1)


TP2<-length(intersect(a,b))/allt
FP2<- length(setdiff(b,a))/allf

















## Birth-death MCMC methods binomial case
bdmcmcbin<-function(X,y){
  dd<-as.data.frame(cbind(X,y))
  p<-dim(X)[2]
  nei<-matrix(0,501,p)
  w<-matrix(0,500,1)
  temp<-nei[1,]
  index<-which(temp>0)
  if (length(index)==0){
    fit<-glm(y~1,data=dd,family = binomial(link = "logit"))  
  } else{
    fit<-glm(y~0+.,data=dd[,c(index,p+1)],family=binomial(link = "logit"))  
  }
  pg0<-BIC(fit)
  len0<-length(index)
  
  for (k in 1:500){

    pg<-matrix(0,1,p)
    len<-rep(0,p)
    
    temp2<-foreach (i=1:p, .combine = cbind) %dopar%{
      library(expm)
      temp<-nei[k,]
      temp[i]<-abs(temp[i]-1)
      index<-which(temp>0)
      if (length(index)==0){
        fit<-glm(y~1,data=dd,family = binomial(link = "logit"))  
      } else{
        fit<-glm(y~0+.,data=dd[,c(index,p+1)],family = binomial(link = "logit"))  
      }
      pg[i]<-BIC(fit)
      len[i]<-length(index)
      rbind(pg[i],len[i])
    }
    pg<-temp2[1,]
    len<-temp2[2,]
    if (max(pg0-pg)<500){
      pg1<-exp(0.5*(pg0-pg))
    } else{
      pg1<-exp((pg0-pg)/10)
    }
    
    prior<-0.1^temp2[2,]/0.1^len0
  
    pg1<-pg1*prior
    #pg1[pg1>1]<-1
    w[k]<-1/sum(pg1)
    pg1<-pg1/sum(pg1)
    sn<-sample(x=c(1:p),size=1,prob=pg1)
    nei[k+1,]<-nei[k,]
    nei[k+1,sn]<-abs(nei[k+1,sn]-1)
    pg0<-pg[sn]
    len0<-length(which(nei[k+1,]>0))
   
    
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
