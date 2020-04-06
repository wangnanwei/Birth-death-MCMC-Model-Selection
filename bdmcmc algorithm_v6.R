
printIndexes = function(ind,maxV) {
  n = maxV
  v = rep("-",  n)
  
  v[ind] = "S"
  
  v = paste(v,collapse="")
  
  #cat(v,"\n");
  v  
  
}


modelVariablesToString = function(ind, n) {
  
  if(length(ind)==0) { return("N0"); }
  
  paste(ind,collapse="-")
}




debug_all_model_evaluations = 0
model_evaluations = c()

model_evaluations_lookup = new.env()








## log-likelihood for lm objects (unweighted version)
logLik_lm_unweighted_apo <- function(object)
{
  res <- object$residuals # not resid(object) because of NA methods
  N <- length(res)
  #val <- .5* ( - N * (log(2 * pi) + 1 - log(N) +   log(  sum(res^2)  ) ) )
  sig=sum(object$residuals^2)/object$df.residual
  val=N*log(1/(sqrt(2*pi*sig)))- sum(object$residuals^2)/(2*sig)
  val
}






## Solve by cholesky decomposition (version 1)


lmFastChol_control = 
  function (B, symmetric = TRUE, tol.values = 1e-07, tol.vectors = 1e-07) 
  {
    
    A <- chol(B, pivot = TRUE)
    pivot <- attributes(A)$pivot
    rank <- attributes(A)$rank
    ok <- sort(pivot[1:rank])
    nok <- if (rank < length(pivot)) 
      pivot[(rank + 1):length(pivot)]
    else NULL
    B <- B[ok, ok]
    
    list(XTX = B, rank = rank, pivot = c(ok, nok))
  }





lmFastChol_apo = function (y, X, eigendec = FALSE, tol.solve = .Machine$double.eps ) 
{
    nvar <- ncol(X)
    nobs <- nrow(X)
    
    tol.values = 1e-07
    tol.vectors = 1e-07
    
    # eigendec = FALSE
    
    
    A <- crossprod(X)
    
    Xy <- t(crossprod(y, X))
    
    yy <- crossprod(y)
    
    
    if (eigendec)
    {  
      ris <- lmFastChol_control(A, , tol.values, tol.vectors) 
    }  else { 
      ris <- list(XTX = A, rank = nvar, pivot = 1:nvar)
    }
    
    #ris$XTX <- as(ris$XTX, "matrix")
    ok <- ris$pivot[1:ris$rank]
    #coef <- as(solve(ris$XTX, Xy[ok]), "numeric")
    coef <- solve(ris$XTX, Xy[ok] )
    
    coefficients <- rep(NA, nvar)
    coefficients[ok] <- coef
    #RSS <- yy - 2 * crossprod(coef, Xy[ok]) + t(coef) %*%  ris$XTX %*% coef
    
    
    coefficients
}



lmFastChol_logLik <- function(X,y)
{
  c1 = lmFastChol_apo(y,X)
  
  res = (  y - X %*% c1 )
  
  N <- length(res)
  dfresidual=N-length(c1)
  #val <- .5* ( - N * (log(2 * pi) + 1 - log(N) +   log(  sum(res^2)  ) ) )
  sig=sum(res^2)/dfresidual
  val=N*log(1/(sqrt(2*pi*sig)))- sum(res^2)/(2*sig)
  
  val
}




## Method 3: the lmFast function directly solves the system

lmFastWithResiduals = function(X,y) {
  
  m = solve(crossprod(X), crossprod(X,y))
  #r = y - X %*% m
  #r
  m
}



lmFast_loglik = function(X,y) {
  
  #res = lmFastWithResiduals(X,y)
  m=lmFastWithResiduals(X,y)
  res = y - X %*% m
  N <- length(res)
  dfresidual=N-length(m)
  #val <- .5* ( - N * (log(2 * pi) + 1 - log(N) +   log(  sum(res^2)  ) ) )
  sig=sum(res^2)/dfresidual
  val=N*log(1/(sqrt(2*pi*sig)))- sum(res^2)/(2*sig)
  val
  
}










bdmcmc  <-  function(X , y , e , N ,  Nburn  ,   priortype   ,  method = 2    )   {

  
  #Arguments:
  # X: n*p covariates matrix;
  #y: n*1 continous response variable;
  # priortype, e: priortype takes binary value {0,1}, e takes a vlue from [0,1], control the sparsity of the selected model.
  
  
  n = dim(X)[1]   # sample size
  p = dim(X)[2]   # number of variables
  
  nei   = matrix(0,N+Nburn+1,p)   # matrix to store the results of MCMC iterations
  w     = matrix(0,N,1)   # Posterior probability of the selected model in each MCMC iteration.
  temp  = nei[1,] 
  index = which(temp>0)    # variable "index" return the varialbes in the starting model, here it is empty as we start from empty model.
  

  
  # Fit a linear regression of "y" given the selected variables.
  
  if ( length(index)==0 )  {
    
    fit = lm( y ~ 1 )  
    
  } else{
    
    X1 = as.matrix( X[,index] )
    fit = lm( y ~ X1 )  
    
  }
  
  
  
  
  
  
  
  #pg0<-BIC(fit)
  
  
  len0 = length(index) # the number of variables in the initial model "fit"
  pg0  = -2*as.numeric(logLik(fit))+0.5*log(p)*(len0+1)+(len0+1)*log(n) # Extended BIC value of the model "fit"
  # pg0 = -2*as.numeric(logLik(fit))+2*log(p+2)*(len0+2)+(len0+2)*log(n)
  
  
  
  
  intercept_matrix = matrix( 1, ncol = 1, nrow = n )
  
  X = as.matrix(X)
  
  X = cbind(X,intercept_matrix)  # Add intercept as last column so we can avoid cbind later
  

  
  model_evaluations_lookup <<- new.env()
  
  
  
  ### Main Loop ###
  
  for (k in   1  :  (N + Nburn)  )  {
    
    cat("Iteration ",k,"\n")

    pg = matrix(0,1,p)  # initiate 1*p vector to store the BIC value of all the candiate models 
    
    len = rep(0,p)      # initiate 1*p vector to store the number of variables of all the candiate models
    
    
  
    
    
      
    
    #temp2<-foreach (i=1:p, .combine = cbind) %dopar%{
    temp2 = sapply( 1:p,  function(i) {  
      
    # using parallel computing to get the BIC values of all the candiates.
      
    temp = nei[k,] # temp:the variables in the selected model of the kth iteration.
    temp[i] = abs(temp[i]-1)# remove or add the ith variable to get a candiate depending on if the ith variable is in the kth model or not.
    
    
    index = which(temp>0) # index: variables selected in the candiate model.
    
    
    

    
    
    #cat(printIndexes(index, p ),"\n")
    
    if(debug_all_model_evaluations > 0) {
        
      
        a_ind = printIndexes(index, p )
        
        model_evaluations <<- c(  model_evaluations  ,  a_ind) 
        
        if(0) {  
            tt = model_evaluations[a_ind]
            if(is.null(tt) || is.na(tt) ) { 
              model_evaluations[a_ind] <<- 1; 
            } else { 
              model_evaluations[a_ind] <<- model_evaluations[a_ind] + 1; 
            }
        }
    }
    
    
    
    
    
    #cat(index,"\n")
    
    # fit a linear regression
    
    
    
    ## Check previous fitted models so we avoid re-evaluating
    model_logLik = NA

    a_ind = modelVariablesToString(index, p )
    
    previous_evaluation = model_evaluations_lookup[[ a_ind ]] 
    if( !is.null(  previous_evaluation  ) )  model_logLik = previous_evaluation
    
    if(0) if(!is.na(model_logLik)) { cat("resuning model: log=", model_logLik, "\n"); if(0) print(as.list(model_evaluations_lookup)) }


    if( is.na(model_logLik) ) {
          

          if (  length(index) == 0 ) {
      
            # fit = lm( y ~ 1 )  
            
            if(method==1) model_logLik = logLik_lm_unweighted_apo(  lm.fit(intercept_matrix,y) )
            
            if(method==2) model_logLik = lmFastChol_logLik( intercept_matrix , y )
            
            if(method==3) model_logLik = lmFast_loglik( intercept_matrix , y )
            
          } else {
            
            # X1 = as.matrix( X[,index] )
            # fit = lm(y~X1) 
            
            m1 = X[, c( index, p+1 ) ]   # Assume last index of matrix is intercept
            
            #m1 = cbind(intercept_matrix, X[,index] )
            if(method==1) model_logLik = logLik_lm_unweighted_apo( lm.fit(m1 , y ) )
            
            if(method==2) model_logLik = lmFastChol_logLik(m1,y)
            if(method==3) model_logLik = lmFast_loglik( m1 , y )
            
            
            
          }
          
          
          
          model_evaluations_lookup[[ a_ind ]] = model_logLik
          
    }
    
      
    
    # cat("model loglik=", model_logLik,"\n")
    
    #pg<-BIC(fit)
      
    len = length(index) # compute the number of variables
  
    #pg = -2*as.numeric(logLik(fit))+2*log(p+2)*(len+2) # return the BIC value of fit
  
    pg = -2*model_logLik+0.5*log(p)*(len+1) +(len+1)*log(n)# return the BIC value of fit
    
    
    #pg = -2*as.numeric(  runif(1,-1000,-700)  )+2*log(p+2)*(len+2) # return the BIC value of fit
      
      
      #pg<--2*as.numeric(logLik(fit))+2*log(p+2)*(len+2)+(len+2)*log(n)
      
      rbind(pg,len)# return the output in the parallel computing code
    } )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # cat("temp2=",temp2)
    # cat("finished loop1\n")
    # str(temp2)
    
    
    pg<-temp2[1,] # BIC of all the candiates.
    
    len<-temp2[2,]# the number of variables in all the candiates.
    if (max(pg0-pg)>1000){
      pg1<-exp(0.05*(pg0-pg)) # compute part of the Bayes factor
    }else{
      pg1<-exp(0.5*(pg0-pg)) # compute part of the Bayes factor
      
    }
    
    ## The prior
    ck<-len*log(exp(1)*p/len)+2*log(len+2)
    ck[is.na(ck)]<-1+2*log(2)
    if (len0==0) {
      ck0<-1+2*log(2)
    }    else{
      ck0<-len0*log(exp(1)*p/len0)+2*log(len0+2)
    }
    
    if (priortype==0) {
      prior<-e^temp2[2,]/e^len0
    } else{
      prior<-exp(-e*ck)/exp(-e*ck0)
    }
    ## the end of choose a prior
    

    
    pg1<-pg1*prior # Bayes factors of all the candiates compare to the kth model.
    #pg1[pg1>1]<-1
    w[k]<-1/sum(pg1)
    pg1<-pg1/sum(pg1) # normalize the bayes factors to get a probability vector.
    

    #print(prior)
    #print(sum(pg1))
    #print(pg1)
    
    sn<-sample(x=c(1:p),size=1,prob=pg1) # sample one varialbe to add or remove based on the probability vector "pg1"
    
    
    nei[k+1,] = nei[k,] 
    nei[k+1,sn] = abs(nei[k+1,sn]-1)# update the new variables in the select model, i.e. add or remove one variable.
    pg0 = pg[sn] # save the BIC value of the selecte model for the next iteration
    
    
    len0 = length(which(nei[k+1,]>0)) # save the number of variables of the selected model for the next iteration
   
    
    
    if(1) {
      cat("  BIC value of selected model: ", sprintf("%10.3f",pg0), "\n" )# print BIC values of the select model in the kth iteration
      vsel = which(nei[k+1,]>0)
      cat("  Variables selected: " ,  vsel ,   "\n" )
    }
    
    
    
  }
  
  
  
  
  

  
  
  
  result  =   cbind(nei[(Nburn+1):(N+Nburn),] ,w[(Nburn+1):(N+Nburn)]) # save the MCMC result after burning Nburn sample points
  
  result[,p+1] =  result[,p+1] / sum(result[,p+1]) # the p+1 column of the matrix "result" is the waiting time of the model 
  # in the same row, normalize it to get the probability of this model
  prob<-NULL
  for(i in 1:p){
    index<-which(result[,i]>0)
    prob[i]<-sum(result[index,p+1])
  }# Bayesian averaging method to get the posterior probability of each variables belongs to the true model.
  
  return(prob)
}






## Birth-death MCMC methods for logistic regression
bdmcmcbin  <-  function(X,y,e,N,Nburn,priortype, method = 1)   {
  
  glm_family = binomial(link = "logit")
  
  n<-dim(X)[1] # sample size
  
  p<-dim(X)[2] # number of variables
  
  nei<-matrix(0,N+Nburn+1,p) # matrix to store the results of MCMC iterations
  
  w<-matrix(0,N,1)  # Posterior probability of the selected model in each MCMC iteration.
  
  temp<-nei[1,]
  
  index <- which(temp>0)   # variable "index" return the varialbes in the starting model, here it is empty as we start from empty model.
  
  
  
  
  
  if (length(index)==0) {
    
    fit <- glm(y~1,glm_family )  
    
  } else{
    
    X1 <- as.matrix(X[,index])
    
    fit  <-  glm(y~X1,glm_family)  
    
  }
  
  
  
  #pg0<-BIC(fit)
  
  
  len0 = length(index)
  
  
  #pg0<--2*as.numeric(logLik(fit))+2*log(p+1)*(len0+1)
  pg0 = -2*as.numeric(logLik(fit))+0.5*log(p)*(len0+1)+(len0+1)*log(n)
  
  
  
  intercept_matrix = matrix( 1, ncol = 1, nrow = n )
  
  X = as.matrix(X)
  
  X = cbind(X,intercept_matrix)  # Add intercept as last column so we can avoid cbind later
  
  #a_X <<- X
  #a_p <<- p
  
  
  model_evaluations_lookup <<- new.env()
  
  
  best_guess_coef = rep(NA, p+1)
  
  for(i in 1:p) {    
    m1 = X[, c(i,p+1) ];  
    fit = glm.fit(m1,y, family = glm_family) 
    best_guess_coef[i] = fit$coef[1]
  }
  
  best_guess_coef[p+1] = coef( glm(y~1,glm_family ) ) [1]
  
    
  
  
  for (k in   1 : (N + Nburn)  )  {
    
    cat("Iteration ", k , "\n")
    
    pg <- matrix(0,1,p)
    len <- rep(0,p)
    
    #temp2  <-  foreach (i=1:p, .combine = cbind) %dopar%{
    temp2 = sapply( 1:p,  function(i) {  
        
      temp  =  nei[k,]
      temp[i]  =  abs(temp[i]-1)

      
      index  =  which(temp > 0)
      
      
      
      
      model_logLik = NA
      
      
      a_ind = modelVariablesToString(index, p )
      
      previous_evaluation = model_evaluations_lookup[[ a_ind ]] 
      if( !is.null(  previous_evaluation  ) )  {  model_logLik = previous_evaluation; }
      
      
      
      if( is.na(model_logLik) ) {
        
        
        if (  length(index) == 0 ) {
          
          # fit = lm( y ~ 1 )  
          
          model_logLik = logLik( glm ( y ~ 1 , family = glm_family ) )
          
        } else {

          #cat("indexes=",index," p+1=",p+1,"\n");
          X=cbind(X,rep(1,n))
          vars = c(index, p+1) # Adding intercept
          
          m1 = X[, vars ]   # Assume last index of matrix is intercept
          
          
          
          start = best_guess_coef[vars]
          start[is.na(start)] = 0 
          
          # runif(length(is.na(start)), -1,1)
          
          

          
          # X1 <- as.matrix(X[,index])
          # fit = glm(y ~ X1,family = glm_family) 
          # model_logLik = logLik(fit)
          
          
          #fit = glm.fit3(m1,y, family = glm_family)
          
          # There was no benefit found in using speedglm
          #fit = speedglm.wfit(y,m1, family=glm_family)
          
          
          if(method==1) {
            
            fit = glm.fit(m1,y, family = glm_family)

          } else if(method==2) {
            
            fit = glm.fit(m1,y, family = glm_family, start = start, control=glm.control(epsilon = 1e-4) )
            
          }
          
          model_logLik = fit$rank  - fit$aic/2
          
          
          # best_guess_coef = rep(0 , length(best_guess_coef) )
          best_guess_coef[ vars ] = fit$coef
          # cat(best_guess_coef,"\n")
          
          #  cat("model_logLik = ", model_logLik,model_logLik2, "\n")
          #  a_y <<- y;a_m1 <<- m1;a_X1 <<- X1;a_fam <<- glm_family
          
        }
        
        model_evaluations_lookup[[ a_ind ]] = model_logLik
      }
      
      
      #cat("loglik = ", model_logLik,"\n")
      
      
      #pg<-BIC(fit)
      len = length(index)

      
      if( is.na(model_logLik) ) { cat("model log lik is NA ", index, len, "\n"); stop(); }
      
            
      #pg<--2*as.numeric(logLik(fit))+2*log(p+1)*(len+1)
      pg =  -2*as.numeric(   model_logLik    )+0.5*log(p)*(len+1)+(len+1)*log(n)
      
      if( is.na(pg) ) { cat("pg is NA\n"); stop(); }
      
      rbind(pg,len)
    } )
    
    
    
    pg<-temp2[1,]
    
    #if (max(pg)>5000) {
     # pg=pg/100
    #  pg0=pg0/100
    #}
    
    len<-temp2[2,]
    #pg1<-exp(0.5*(pg0-pg))
    if (max(pg0-pg)>1000){
      pg1<-exp(0.05*(pg0-pg)) # compute part of the Bayes factor
    }else{
      pg1<-exp(0.5*(pg0-pg)) # compute part of the Bayes factor
      
    }
    
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






## Birth-death MCMC methods for poisson regression
bdmcmcpos<-function(X,y,e,N,Nburn,priortype){
  n<-dim(X)[1]
  p<-dim(X)[2]
  nei<-matrix(0,N+Nburn+1,p)
  w<-matrix(0,N,1)
  temp<-nei[1,]
  index<-which(temp>0)
  if (length(index)==0){
    fit<-glm(y~1,family =poisson(link = "log"))  
  } else{
    X1<-as.matrix(X[,index])
    fit<-glm(y~X1,family = poisson(link = "log"))  
  }
  #pg0<-BIC(fit)
  len0<-length(index)
  #pg0<--2*as.numeric(logLik(fit))+2*log(p+1)*(len0+1)
  pg0<--2*as.numeric(logLik(fit))+0.5*log(p)*(len0+1)+(len0+1)*log(n)
  
  
  for (k in 1:(N+Nburn)){
    # print(k)
    pg<-matrix(0,1,p)
    len<-rep(0,p)
    
    temp2<-sapply(1:p, function(i){
      
      temp<-nei[k,]
      temp[i]<-abs(temp[i]-1)
      index<-which(temp>0)
      if (length(index)==0){
        fit<-glm(y~1,family = poisson(link = "log"))  
      } else{
        X1<-as.matrix(X[,index])
        fit<-glm(y~X1,family =poisson(link = "log")) 
      }
      #pg<-BIC(fit)
      len<-length(index)
      #pg<--2*as.numeric(logLik(fit))+2*log(p+1)*(len+1)
      pg<--2*as.numeric(logLik(fit))+0.5*log(p)*(len+1)+(len+1)*log(n)
      
      rbind(pg,len)
    })
    
    pg<-temp2[1,]
    if (max(pg0-pg)>1000){
      pg1<-exp(0.05*(pg0-pg)) # compute part of the Bayes factor
    }else{
      pg1<-exp(0.5*(pg0-pg)) # compute part of the Bayes factor
      
    }
    len<-temp2[2,]
    #pg1<-exp(0.5*(pg0-pg))
    
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
    #cat("BIC value of selected model: ", sprintf("%10.3f",pg0), "\n" )
    #print(which(nei[k+1,]>0))
    
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

