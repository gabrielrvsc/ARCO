estimate_lambda0=function(data,fn,p.fn,start=0.3,end=0.95,treated.unity=1,lag=0,Xreg=NULL){
  
  for(i in 1:length(data)){
    if(is.null(colnames(data[[i]]))){
      colnames(data[[i]])=paste("V",i,"-U",1:ncol(data[[i]]),sep="")
      cat("Variable names not informed. Automatic names supplied.")
    }
  }
  
  for(i in 1:length(data)){
    aux=length(unique(colnames(data[[i]])))
    K=ncol(data[[i]])
    if(aux<k){
      colnames(data[[i]])=paste("V",i,"-U",1:ncol(data[[i]]),sep="")
      cat("Some variables had no name. Automatic names supplied.")
    }
  }
  
  
  T=nrow(data[[1]])
  
  starting.point=floor(start*T)
  ending.point=floor(end*T)
  
  save.delta=matrix(0,T,length(data))
  for(i in starting.point:ending.point){
    m=fitARCO(data=data,fn=fn,p.fn = p.fn,treated.unity=treated.unity,lag=lag, t0=i,Xreg=Xreg,display = FALSE)
    delta=m$delta[2]
    save.delta[i,]=delta
  }
  
  delta.norm=sqrt(rowSums((save.delta)^2))
  lambda0=which(delta.norm==max(delta.norm))
  delta=delta.norm[lambda0]
  
  return(c("lambda0"=lambda0,"delta"=delta))
}
