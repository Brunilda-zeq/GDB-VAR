#The inputs are: the number of lags (p) of the model,the data (data) we are
#using for estimation, and the number of observations in the data set (n)
#The output is: OLS estimates of the paremeters (bhat), residuals
#of the estimation (resids), estimate of the variance of the error term
#(sigmahat). The estimates are the same with those one would get via
#Maximum Likelihood estimation conditional on the first p observations.
# The function works for both AR and VAR models


ols_est<-function(data,p,n){
  
  Y<-data[(p+1):n,]
  
  X<-matrix(1,nrow = n-p,ncol = 1)  #add intercepts to the model
  
  for(i in 1:p){
    
    X<-cbind(X,data[(p+1-i):(n-i),])
  }
  
  bhat<-solve(t(X) %*% X) %*% t(X) %*% Y
  resids<-Y-X%*%bhat
  sigmahat<-var(resids)
  
  output<- list(bhat=bhat, resids=resids, sigmahat=sigmahat)
  
  return(output)
  
  
}



