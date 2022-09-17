

irf_fun<- function(companion,k,p,H,hor){
  
  J<-cbind(diag(k),matrix(0,nrow = k,ncol = (k*(p-1))))
  
  irf<-matrix(H,nrow=k^2,ncol = 1)
  
  for (j in 1:hor){
    
    A1<-J%*%(companion%^%j)%*%t(J)%*%H
    A2<-matrix(A1,nrow = k^2,ncol = 1)
    irf<-cbind(irf,A2)
    
  }
  
  return(irf)
  
}





