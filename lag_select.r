#This function gives the optimal lag lenghth for VAR. The input is the maximum
#number of lags that are checking over(pmax). The output is the optimal
#number of lags accoding to Akaike(AIC), Shwarz-Bayesian (SBIC), Hannan-Quinn
#(HQIC) information criteria

lag_select<- function(pmax){
  
  mat<-matrix(NA,nrow = 3,ncol = pmax)
  
  for (p in 1:pmax){
    
    sigmahat<-ols_est(data,p,n)$sigmahat
    
    kk<-(k^2)*p +k  #total number of parameters in all equations
    
    AIC<-log(det(sigmahat)) + 2*kk/n
    SBIC<-log(det(sigmahat)) + kk*log(n)/n
    HQIC<- log(det(sigmahat)) + 2*kk*log(log(n))/n
    IC<-rbind(AIC,SBIC,HQIC)
    mat[,p]<-IC
    
  }
  
  score<-matrix(apply(mat,1,which.min))
  colnames(score)<- "number of lags"
  rownames(score)<- c("AIC","SBIC","HQIC")
  
  return(score)
  
}



