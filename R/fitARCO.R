fitARCO=function(data,fn,p.fn,treated.unity,t0,lag=0,Xreg=NULL,display=TRUE,HACweights=1,alpha=0.05,...){
  
  if(is.null(names(date))){
    names(data)=paste("Variable",1:length(date),sep="")
    cat("The data list was unnamed. Automatic names supplied. \n")
  }
  
  
  for(i in 1:length(data)){
    if(is.null(colnames(data[[i]]))){
      colnames(data[[i]])=paste("V",i,"-U",1:ncol(data[[i]]),sep="")
      cat("Variable names not informed. Automatic names supplied. \n")
    }
  }
  
  for(i in 1:length(data)){
    aux=length(unique(colnames(data[[i]])))
    k=ncol(data[[i]])
    if(aux<k){
      colnames(data[[i]])=paste("V",i,"-U",1:ncol(data[[i]]),sep="")
      cat("Some variables had no name. Automatic names supplied. \n")
    }
  }
  
  if(length(data)==1){
    Y=matrix(data[[1]][,treated.unity],ncol = 1)
    X=data[[1]][,-treated.unity]
  }else{
    Y=Reduce("cbind",lapply(data, function(x) x[,treated.unity]))
    X=Reduce("cbind",lapply(data, function(x) x[,-treated.unity]))
  }
  Y.raw=Y 
  
  if(lag!=0){
    X=embed(X,lag+1)
    Y=tail(Y,nrow(X))
  }
  
  if(length(Xreg)!=0){##mudou
    X=cbind(X,tail(Xreg,nrow(X)))
  }
  
  if(is.vector(Y)){
    Y=matrix(Y,length(Y),1)
  }
  
  y.fit=matrix(Y[1:(t0-1-lag),],ncol=length(data))
  y.pred=matrix(Y[-c(1:(t0-1-lag)),],ncol=length(data))
  x.fit=X[1:(t0-1-lag),]
  x.pred=X[-c(1:(t0-1-lag)),]
  
  save.cf=matrix(NA,nrow(y.pred),length(data))
  save.fitted=matrix(NA,nrow(Y),length(data))
  model.list=list()
  for(i in 1:length(data)){
    model=fn(x.fit,y.fit[,i])
    model.list[[i]]=model
    contra.fact=p.fn(model,x.pred)
    save.cf[,i]=contra.fact
    save.fitted[,i]=p.fn(model,X)
  }
  
  # == intervalo == ##
  delta.aux=tail(Y.raw,nrow(save.cf))-save.cf
  delta=colMeans(delta.aux)
  aux=matrix(0,nrow(X),length(data))
  aux[t0:nrow(aux),]=1
  vhat=Y-(save.fitted+aux*delta)
  v1=matrix(vhat[1:(t0-1),],ncol=length(data))
  v2=matrix(vhat[t0:nrow(vhat),],ncol=length(data))
  
  tau01=cov(v1)*HACweights[1]
  tau02=cov(v2)*HACweights[1]
  tauk1bart=0
  tauk2bart=0
  M=length(HACweights)-1
  if(M>0){
    for(k in 1:M){
      aux=cov(v1[(1+k):nrow(v1),],v1[1:(nrow(v1)-k),])
      tauk1bart=tauk1bart+aux%*%t(aux)* HACweights[i+1] # kernel
      aux=cov(v1[(1+k):nrow(v1),],v1[1:(nrow(v1)-k),])
      tauk2bart=tauk2bart+ aux%*%t(aux) * HACweights[i+1] # bartlett kernel
    }
  }
  tauT1=tau01+tauk1bart
  tauT2=tau02+tauk2bart
  sigmahat=tauT1/(t0-1) + tauT2/ (nrow(X)-t0+1)
  w=sqrt(diag(sigmahat))
  uI=delta+(w*qnorm(1-alpha/2))/sqrt(nrow(X))
  lI=delta-(w*qnorm(1-alpha/2))/sqrt(nrow(X))
  delta.stat=cbind("LB"=lI,"delta"=delta,"UB"=uI)
  
  # == plot == #
  if(display==TRUE){
    par(mfrow=c(1,length(data)))
    for(i in 1:ncol(Y.raw)){
      plot(Y.raw[,i],type="l",ylab=paste("Y",i,sep=""),xlab="Time",...)
      lines(c(rep(NA,t0-2),Y.raw[t0-1,i],save.cf[,i]),col=2)
      abline(v=t0,col="blue",lty=2)
    }
  }
  
  ## == Names == ##
  names(model.list)=names(data)
  colnames(cf)=names(data)
  rownames(cf)=tail(rownames(Y.raw),nrow(cf))                          
                            
  return(list("cf"=save.cf,"model"=model.list,"delta"=delta.stat))
}
