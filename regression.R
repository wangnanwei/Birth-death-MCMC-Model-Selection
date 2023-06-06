library(SOIL)
library(network)
library(ggplot2)
library(GGally)


dd_exp_all=read.table("Expression_data_all.csv", header = TRUE,sep = ",")

dd_exp_all=t(dd_exp_all)
dd_exp_all=as.data.frame(dd_exp_all)

colnames(dd_exp_all)=dd_exp_all[2,]
dd_exp_all= dd_exp_all[-c(1:2),]
dd1=dd_exp_all
dd1 <- as.data.frame(sapply(dd1, as.numeric))

dd1=log(dd1)
dd1=scale(dd1)
dd1=as.data.frame(dd1)

dd1=dd1[,which(colSums(is.na(dd1))==0)]

rownames(dd1)<-rownames(dd_exp_all)
dd_pam=read.table("Expression_data_normalized_PAM50.csv", header = TRUE,sep = ",")

dd_pam=t(dd_pam)
dd_pam=as.data.frame(dd_pam)

colnames(dd_pam)=dd_pam[2,]

dd_pam= dd_pam[-c(1:2),]
rn<-rownames(dd_pam)
dd_pam<- as.data.frame(sapply(dd_pam, as.numeric))
rownames(dd_pam)<-rn
dd_pam=log(dd_pam+1e-20)
dd_pam=as.data.frame(scale(dd_pam))
dd_mut=read.table("Mutation_file_new.csv",header = TRUE,sep=",")

dd_mut=t(dd_mut)
dd_mut=as.data.frame(dd_mut)
colnames(dd_mut)<-dd_mut[1,]
dd_mut=dd_mut[-1,]
dd_mut=as.matrix(dd_mut)
 dd_mut[which(dd_mut=="DEL")]<-1
dd_mut[which(dd_mut=="INS")]<-1
dd_mut[which(dd_mut=="SNP")]<-0
dd_mut[which(dd_mut=="0")]<-0
dd_mut=as.data.frame(dd_mut)
for( i in 1:4864){
  dd_mut[,i]=as.numeric(dd_mut[,i])
}
index=NULL
for( i in 1:4864){
if (sum(dd_mut[i])>=5) index=c(index,i)  
}
dd_mut=dd_mut[,index]


dd_subtype=read.table("Subtype_info.csv",header = TRUE,sep=",")
dd_subtype=dd_subtype[complete.cases(dd_subtype$subtype.info),]
dd_subtype=dd_subtype[!duplicated(dd_subtype$X.1),]
rn<-dd_subtype$X.1
rn=gsub("-",".",rn)
dd_subtype=dd_subtype[,-c(1,2)]
dd_subtype=as.data.frame(dd_subtype)
rownames(dd_subtype)<-rn
dd_subtype$Basal<-dd_subtype$dd_subtype=="Basal"
dd_subtype$Her2<-dd_subtype$dd_subtype=="Her2"
dd_subtype$LumA<-dd_subtype$dd_subtype=="LumA"
dd_subtype$LumB<-dd_subtype$dd_subtype=="LumB"
dd_subtype$Normal<-dd_subtype$dd_subtype=="Normal"



chooserow=match(rownames(dd_subtype),rownames(dd_mut))

dd_mut=dd_mut[chooserow[complete.cases(chooserow)],]

dd_subtype=dd_subtype[!is.na(chooserow),]
## select correlated gene expression variants:
dd<-merge(dd1,dd_subtype,by="row.names")
rownames(dd)=dd[,1]
dd=dd[,-1]

  

index=NULL
for ( i in 1: (ncol(dd)-ncol(dd_subtype))){
  fit1<-glm(dd$Basal~dd[,i],family = binomial(link="logit"))
  fit2<-glm(dd$Her2~dd[,i],family = binomial(link="logit"))
  fit3<-glm(dd$LumA~dd[,i],family = binomial(link="logit"))
  fit4<-glm(dd$LumB~dd[,i],family = binomial(link="logit"))
  fit5<-glm(dd$Normal~dd[,i],family = binomial(link="logit"))
  pvalue=c(summary(fit1)$coefficients[2,4],summary(fit2)$coefficients[2,4],
           summary(fit3)$coefficients[2,4],summary(fit4)$coefficients[2,4],
           summary(fit5)$coefficients[2,4])
  if (max(pvalue)<0.05){
    index=c(index,i)
  }
}
dd_exp_all=dd[,index]

dd_pam=dd_pam[match(row.names(dd),rownames(dd_pam)),]
dd_subtype1=dd_subtype[,-1]

## model selection
adj=matrix(0,55,55)
temp=cbind(dd_pam,dd_subtype1)
rn=colnames(temp)
for(i in 1:50){
  
  fit1<-SOIL(temp[,-i],temp[,i],family="gaussian",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  nei=adj[i,-i]
  nei[which(fit1$importance>0.5)]<-1
  adj[i,-i]=nei
}
for(i in 1:5){
  
  fit1<-SOIL(dd_pam,dd_subtype1[,i],family="binomial",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:

  adj[(ncol(dd_pam)+i),which(fit1$importance>0.5)]=1
}

adj=adj+t(adj)
adj[adj==1]<-1
net1<-network(adj,directed=FALSE)
network.vertex.names(net1)<-rn

net1 %v% "type"<-c(rep("PAM50 genes",50),rep("BC subtype",5))  

col = c("PAM50 genes" = "green", "BC subtype" = "red")
ggnet2(net1,node.size=6,color="type",palette = col,label = TRUE,label.size = 4)







## model selection with mutation
adj=matrix(0,34,34)
dd_mut1=dd_mut[,colSums(dd_mut)>3]
temp=cbind(dd_mut1,dd_subtype1)
rn=colnames(temp)
for(i in 1:29){
  
  fit1<-SOIL(temp[,-i],temp[,i],family="binomial",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  nei=adj[i,-i]
  nei[which(fit1$importance>0.5)]<-1
  adj[i,-i]=nei
}
for( i in 1:5){
  fit1<-SOIL(dd_mut1,dd_subtype1[,i],family="binomial",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  
  adj[(ncol(dd_mut1)+i),which(fit1$importance>0.5)]=1
}

adj=adj+t(adj)
adj[adj==2]<-1
net1<-network(adj,directed=FALSE)
network.vertex.names(net1)<-rn

net1 %v% "type"<-c(rep("Mutation",29),rep("BC subtype",5))  

col = c("Mutation" = "green", "BC subtype" = "red")
ggnet2(net1,node.size=6,color="type",palette = col,label = TRUE,label.size = 4)




## everything together
temp=cbind(dd_exp_all,dd_pam,dd_mut1)
temp2=cbind(temp,dd_subtype1)
adj=matrix(0,1336,1336)
for( i in 1: (ncol(dd_exp_all)+ncol(dd_pam))){
  fit1<-SOIL(temp2[,-i],temp2[,i],family="gaussian",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  nei=adj[i,-i]
  nei[which(fit1$importance>0.5)]<-1
  adj[i,-i]=nei
}
for( i in  (ncol(dd_exp_all)+ncol(dd_pam)+1): ncol(temp)){
  fit1<-SOIL(temp2[,-i],temp2[,i],family="binomial",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  nei=adj[i,-i]
  nei[which(fit1$importance>0.5)]<-1
  adj[i,-i]=nei
}

for( i in 1:5){
  fit1<-SOIL(temp,dd_subtype1[,i],family="binomial",method="union",weight_type="BIC")
  #the selected variable using SOIL method is:
  
  adj[ncol(temp)+i,which(fit1$importance>0.5)]=1
}


adj=adj+t(adj)
adj[adj==2]<-1
net1<-network(adj,directed=FALSE)
network.vertex.names(net1)<-colnames(temp2)

net1 %v% "type"<-c(rep("non-PAM50 genes",ncol(dd_exp_all)),rep("PAM50 genes",50)
                   ,rep("Mutation",29),rep("BC subtype",5))  

col = c("non-PAM50 genes" = "gold", "PAM50 genes"="blue","Mutation"="green", "BC subtype" = "red")
ggnet2(net1,node.size=6,color="type",palette = col,label = TRUE,label.size = 4)

index=NULL
for (i in 1332:1336){
  index=c(index,which(adj[i,1:ncol(dd_exp_all)]==1))
}
index=index[!duplicated(index)]

index=index[order(index)]
index2=c(index,(ncol(dd_exp_all)+1):1336)
adj2=adj[index2,index2]

net2<-network(adj2,directed=FALSE)
network.vertex.names(net2)<-colnames(temp2)[index2]

net2%v% "type"<-c(rep("non-PAM50 genes",length(index)),rep("PAM50 genes",50)
                   ,rep("Mutation",29),rep("BC subtype",5))  

col = c("non-PAM50 genes" = "gold", "PAM50 genes"="skyblue","Mutation"="green", "BC subtype" = "red")
ggnet2(net2,node.size=6,color="type",palette = col,label = TRUE,label.size = 4)
