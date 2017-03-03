
plot.fitArCo=function(model,ylab=NULL,main=NULL,plot=NULL,ncol=1,display.fitted=FALSE,y.min=NULL,y.max=NULL){

  t0=model$t0
  data=model$data
  fitted=model$fitted
  treated.unity=model$treated.unity
  cf=model$cf
  
  if (length(data) == 1) {
    Y = matrix(data[[1]][, treated.unity], ncol = 1)
  }else {
    Y = Reduce("cbind", lapply(data, function(x) x[, treated.unity]))
  }
  
  aux=nrow(Y)-nrow(fitted)-nrow(cf)
  if(aux!=0){
    fitted=rbind(matrix(NA,aux,ncol(fitted)),fitted)
  }
  
  ##
  
  
  
  if (is.null(plot)) {
    if(is.null(ylab)){
      ylab=ylab = paste("Y", 1:length(data), sep = "")
    }
    if(is.null(y.min)){
      y.min=apply(rbind(Y,cf),2,min,na.rm=TRUE) 
    }
    if(is.null(y.max)){
      y.max=apply(rbind(Y,cf),2,max,na.rm=TRUE)
    } 
    
  } else {
    if(is.null(ylab)){
      ylab=ylab = paste("Y", plot, sep = "")
    }
    if(is.null(y.min)){
      y.min=apply(as.matrix(rbind(Y,cf)[,plot]),2,min,na.rm=TRUE) 
    }
    if(is.null(y.max)){
      y.max=apply(as.matrix(rbind(Y,cf)[,plot]),2,max,na.rm=TRUE)
    } 
  }
  
  
  if(is.null(plot)){
    par(mfrow = c(ceiling(length(data)/ncol), ncol))
    for (i in 1:ncol(Y)) {
      plot(Y[, i], type = "l", ylab = ylab[i], xlab = "Time",main=main[i],ylim=c(y.min[i],y.max[i]))
      lines(c(rep(NA, t0 - 2), Y[t0 - 1, i], cf[, i]), col = "blue")
      abline(v = t0, col = "blue", lty = 2)
      if(display.fitted==TRUE){
        lines(fitted[,i],col="red")
        legend("topleft",legend = c("Observed","Fitted","Counter Fact."),col=c(1,2,4),
               lwd=c(1,1,1),lty=c(1,1,1),bty="n",xjust=1,seg.len = 1,y.intersp=0.5)
      }
    }
  }else{
    par(mfrow = c(ceiling(length(plot)/ncol), ncol))  
    for(i in 1:length(plot)){
      plot(Y[, plot[i]], type = "l", ylab = ylab[i], xlab = "Time",main=main[i],ylim=c(y.min[i],y.max[i]))
      lines(c(rep(NA, t0 - 2), Y[t0 - 1, plot[i]], cf[, plot[i]]), col = "blue")
      abline(v = t0, col = "blue", lty = 2)
      if(display.fitted==TRUE){
        lines(fitted[,plot[i]],col="red")
        legend("topleft",legend = c("Observed","Fitted","Counter Fact."),col=c(1,2,4),
               lwd=c(1,1,1),lty=c(1,1,1),bty="n",xjust=1,seg.len = 1,y.intersp=0.7)
      }
    }
  }

}


