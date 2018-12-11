
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

## Birth-death MCMC methods
bdmcmc<-function(X,y,e,N,Nburn,priortype){
  #Arguments:
  # X: n*p covariates matrix;
  #y: n*1 continous response variable;
  # priortype, e: priortype takes binary value {0,1}, e takes a vlue from [0,1], control the sparsity of the selected model.
  
  
  n<-dim(X)[1]# sample size
  p<-dim(X)[2]# number of variables
  nei<-matrix(0,N+Nburn+1,p)# matrix to store the results of MCMC iterations
  w<-matrix(0,N,1)# Posterior probability of the selected model in each MCMC iteration.
  temp<-nei[1,] 
  index<-which(temp>0) # variable "index" return the varialbes in the starting model, here it is empty as we start from empty model.
  
  if (length(index)==0){
    X1=as.matrix(rep(1,n))
    fit<-lm.fit(X1,y)  
  } else{
    X1<-as.matrix(X[,index])
    fit<-lm.fit(X1,y)  
  }# Fit a linear regression of "y" given the selected variables.
  #pg0<-BIC(fit)
  len0<-length(index) # the number of variables in the initial model "fit"
  sig=sum(fit$residuals^2)/fit$df.residual
  loglikelihood=n*log(1/(sqrt(2*pi*sig)))- sum(fit$residuals^2)/(2*sig)
  pg0<--2*loglikelihood+log(n)*(len0+1)+2*log(p+2)*(len0+1) # Extended BIC value of the model "fit"
   #pg0<--2*as.numeric(logLik(fit))+2*log(p+2)*(len0+2)+(len0+2)*log(n)
  
  
  for (k in 1:(N+Nburn)){
    print(k)# k is the iteration step of the MCMC process
    pg<-matrix(0,1,p)# initiate 1*p vector to store the BIC value of all the candiate models 
    len<-rep(0,p)# initiate 1*p vector to store the number of variables of all the candiate models
    
    temp2<-foreach (i=1:p, .combine = cbind) %dopar%{
# using parallel computing to get the BIC values of all the candiates.
     
       temp<-nei[k,] # temp:the variables in the selected model of the kth iteration.
     
       temp[i]<-abs(temp[i]-1)# remove or add the ith variable to get a candiate depending on if the ith variable is in the kth model or not.
      index<-which(temp>0) # index: return the variables in the candiate model.
      if (length(index)==0){
        X1=as.matrix(rep(1,n))
        fit<-lm.fit(X1,y)  
      } else{
        #X1<-as.matrix(X[,index])
        X1<-as.matrix(cbind(rep(1,n),X[,index]))
        fit<-lm.fit(X1,y)  
      }# Fit a linear regression of "y" given the selected variables.
      #pg0<-BIC(fit)
      len<-length(index) # the number of variables in the initial model "fit"
      sig=sum(fit$residuals^2)/fit$df.residual
      loglikelihood=n*log(1/(sqrt(2*pi*sig)))- sum(fit$residuals^2)/(2*sig)
      pg<--2*loglikelihood+log(n)*(len0+1)+2*log(p+2)*(len+1) # Extended BIC value of the model "fit"
      #pg0<--2*as.numeric(logLik(fit))+2*log(p+2)*(len0+2)+(len0+2)*log(n)
      rbind(pg,len)# return the output in the parallel computing code
    }
    pg<-temp2[1,] # BIC of all the candiates.
    
    len<-temp2[2,]# the number of variables in all the candiates.
    if (max(pg0-pg)>1000){
      pg1<-exp(0.05*(pg0-pg)) # compute part of the Bayes factor
    }else{
     pg1<-exp(0.5*(pg0-pg)) # compute part of the Bayes factor
      
   }
    #pg1<-exp(0.5*(pg0-pg))
    
    ## The prior
    ck<-len*log(exp(1)*p/len)+2*log(len+2)
    ck[is.na(ck)]<-1+2*log(2)
    if (len0==0) {
      ck0<-1+2*log(2)
    }else{
      ck0<-len0*log(exp(1)*p/len0)+2*log(len0+2)
    }
    if (priortype==0){
      prior<-e^temp2[2,]/e^len0
    } else{
      prior<-exp(-e*ck)/exp(-e*ck0)
    }
    ## the end of choose a prior
    
    pg1<-pg1*prior # Bayes factors of all the candiates compare to the kth model.
    pg1[pg1>1]<-1
    w[k]<-1/sum(pg1)
    pg1<-pg1/sum(pg1) # normalize the bayes factors to get a probability vector.
    sn<-sample(x=c(1:p),size=1,prob=pg1) # sample one varialbe to add or remove based on the probability vector "pg1"
    nei[k+1,]<-nei[k,] 
    nei[k+1,sn]<-abs(nei[k+1,sn]-1)# update the new variables in the select model, i.e. add or remove one variable.
    pg0<-pg[sn] # save the BIC value of the selecte model for the next iteration
    len0<-length(which(nei[k+1,]>0)) # save the number of variables of the selected model for the next iteration
    cat("BIC value of selected model: ", sprintf("%10.3f",pg0), "\n" )# print BIC values of the select model in the kth iteration
    print(which(nei[k+1,]>0))# print the variables in the select model in the kth iteration
    
  }
  
  result<-cbind(nei[(Nburn+1):(N+Nburn),] ,w[(Nburn+1):(N+Nburn)]) # save the MCMC result after burning Nburn sample points
  
  result[,p+1]<-result[,p+1]/sum(result[,p+1]) # the p+1 column of the matrix "result" is the waiting time of the model 
  # in the same row, normalize it to get the probability of this model
  prob<-NULL
  for(i in 1:p){
    index<-which(result[,i]>0)
    prob[i]<-sum(result[index,p+1])
  }# Bayesian averaging method to get the posterior probability of each variables belongs to the true model.
  res=list(prob,nei,result)
  return(res)
}






## Birth-death MCMC methods for logistic regression
bdmcmcbin<-function(X,y,e,N,Nburn,priortype){
  n<-dim(X)[1]
  p<-dim(X)[2]
  nei<-matrix(0,N+Nburn+1,p)
  w<-matrix(0,N,1)
  temp<-nei[1,]
  index<-which(temp>0)
  if (length(index)==0){
    fit<-glm(y~1,family = binomial(link = "logit"))  
  } else{
    X1<-as.matrix(X[,index])
    fit<-glm(y~X1,family = binomial(link = "logit"))  
  }
  #pg0<-BIC(fit)
  len0<-length(index)
  #pg0<--2*as.numeric(logLik(fit))+2*log(p+1)*(len0+1)
  pg0<--2*as.numeric(logLik(fit))+2*log(p+2)*(len0+1)+(len0+1)*log(n)
  
  
  for (k in 1:(N+Nburn)){
    print(k)
    pg<-matrix(0,1,p)
    len<-rep(0,p)
    
    temp2<-foreach (i=1:p, .combine = cbind) %dopar%{
      
      temp<-nei[k,]
      temp[i]<-abs(temp[i]-1)
      index<-which(temp>0)
      if (length(index)==0){
        fit<-glm(y~1,family = binomial(link = "logit"))  
      } else{
        X1<-as.matrix(X[,index])
        fit<-glm(y~X1,family = binomial(link = "logit")) 
      }
      #pg<-BIC(fit)
      len<-length(index)
      #pg<--2*as.numeric(logLik(fit))+2*log(p+1)*(len+1)
      pg<--2*as.numeric(logLik(fit))+2*log(p+2)*(len+1)+(len+1)*log(n)
      
      rbind(pg,len)
    }
    pg<-temp2[1,]
    if (max(pg0-pg)>1000){
      pg1<-exp(0.05*(pg0-pg)) # compute part of the Bayes factor
    }else{
      pg1<-exp(0.5*(pg0-pg)) # compute part of the Bayes factor
      
    }
    len<-temp2[2,]
    pg1<-exp(0.5*(pg0-pg))
    
    ck<-len*log(exp(1)*p/len)+2*log(len+2)
    ck[is.na(ck)]<-1+2*log(2)
    if (len0==0) {
      ck0<-1+2*log(2)
    }
    else{
      ck0<-len0*log(exp(1)*p/len0)+2*log(len0+2)
    }
    if (priortype==0){
      prior<-e^temp2[2,]/e^len0
    } else{
      prior<-exp(-e*ck)/exp(-e*ck0)
    }
    
    pg1<-pg1*prior
    #pg1[pg1>1]<-1
    w[k]<-1/sum(pg1)
    pg1<-pg1/sum(pg1)
    sn<-sample(x=c(1:p),size=1,prob=pg1)
    nei[k+1,]<-nei[k,]
    nei[k+1,sn]<-abs(nei[k+1,sn]-1)
    pg0<-pg[sn]
    len0<-length(which(nei[k+1,]>0))
    cat("BIC value of selected model: ", sprintf("%10.3f",pg0), "\n" )
    print(which(nei[k+1,]>0))
    
  }
  
  result<-cbind(nei[(Nburn+1):(N+Nburn),] ,w[(Nburn+1):(N+Nburn)])
  
  result[,p+1]<-result[,p+1]/sum(result[,p+1])
  prob<-NULL
  for(i in 1:p){
    index<-which(result[,i]>0)
    prob[i]<-sum(result[index,p+1])
  }
  
  return(prob)
}






ptm <- Sys.time()
prob<-bdmcmc(X,y,1,50,50,1)
temp<-Sys.time() - ptm
print(temp)

stopCluster(cl)
registerDoSEQ()
