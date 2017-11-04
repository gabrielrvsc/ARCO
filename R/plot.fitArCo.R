#' Plots realized values and the counterfactual estimated by the fitArCo function
#' 
#' Plots realized values and the counterfactual estimated by the fitArCo function. The plotted variables will be on the same level as supplied to the fitArCo function.  
#' 
#' @param x An ArCo object estimated using the fitArCo function.
#' @param ylab n dimensional character vector, where n is the length of the plot argument or n=q if plot=NULL.
#' @param main n dimensional character vector, where n is the length of the plot argument or n=q if plot=NULL.
#' @param plot n dimensional numeric vector where each element represents an ArCo unit. If NULL, all units will be plotted. If, for example, plot=c(1,2,5) only units 1 2 and 5 will be plotted according to the order specified by the user on the fitArCo.
#' @param ncol Number of columns when multiple plots are displayed.    
#' @param display.fitted If TRUE the fitted values of the first step estimation are also plotted (default=FALSE). 
#' @param y.min n dimensional numeric vector defining the lower bound for the y axis. n is the length of the plot argument or n=q if plot=NULL
#' @param y.max n dimensional numeric vector defining the upper bound for the y axis. n is the length of the plot argument or n=q if plot=NULL
#' @param ... Other graphical parameters to plot.  
#' @param confidence.bands TRUE to plot the counterfactual confidence bands (default=FALSE). If the ArCo was estimated without bootstrap this argument will be forced to FALSE.
#' @param alpha Significance level for the confidence bands. 
#' @export
#' @examples 
#' ##############################################
#' ## === Example based on the q=1 fitArCo === ##
#' ##############################################
#'# = First unit was treated on t=51 by adding
#'# a constant equal to one standard deviation
#' data(data.q1)
#' data=list(data.q1) # = Even if q=1 the data must be in a list
#' ## == Fitting the ArCo using linear regression == ##
#' # = creating fn and p.fn function = #
#' fn=function(X,y){
#' return(lm(y~X))
#' }
#' p.fn=function(model,newdata){
#' b=coef(model)
#' return(cbind(1,newdata) %*% b)}
#' ArCo=fitArCo(data = data,fn = fn, p.fn = p.fn, treated.unit = 1 , t0 = 51)
#' plot(ArCo)
#' @seealso \code{\link{fitArCo}}

plot.fitArCo=function(x,ylab=NULL,main=NULL,plot=NULL,ncol=1,display.fitted=FALSE,y.min=NULL,y.max=NULL,confidence.bands=FALSE,alpha=0.05,...){
  oldmar=graphics::par()$mar
  oldmfrow=graphics::par()$mfrow
  oldoma=graphics::par()$oma
  t0=x$t0
  data=x$data
  fitted=x$fitted.values
  treated.unit=x$treated.unit
  cf=x$cf
  boot.cf=x$boot.cf
  
  if(typeof(boot.cf)=="list"){
    NAboot=Reduce(sum,boot.cf)
    if(is.na(NAboot)){
      warning("NA values on the bootstrap counterfactual: unable to plot confidence bands.")
      confidence.bands=FALSE
    }
  }
  
  if(typeof(boot.cf)!="list"){
    if(confidence.bands==TRUE){
      confidence.bands=FALSE
      cat("Confidence bands set to FALSE: The model has no confidence bands. Set boot.cf to TRUE on the ArCo estimation.")
    }
  }
  
  if (length(data) == 1) {
    Y = matrix(data[[1]][, treated.unit], ncol = 1)
  }else {
    Y = Reduce("cbind", lapply(data, function(x) x[, treated.unit]))
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
    if(display.fitted==TRUE){
      graphics::par(oma = c(4, 1, 1, 1))
    }
    graphics::par(mfrow = c(ceiling(length(data)/ncol), ncol))
    for (i in 1:ncol(Y)) {
      graphics::plot(Y[, i], type = "l", ylab = ylab[i], xlab = "Time",main=main[i],ylim=c(y.min[i],y.max[i]),...)
      
      if(confidence.bands==TRUE){
        aux=Reduce("cbind",lapply(boot.cf,function(x)x[,i]))
        aux1=t(apply(aux,1,sort))
        intervals=c(round(ncol(aux)*(alpha/2)),round(ncol(aux)*(1-alpha/2)))
        xcord=t0:length(Y[,i])
        ycord1=aux1[,intervals[1]]
        ycord2=aux1[,intervals[2]]
        graphics::polygon(c(rev(xcord),xcord),c(rev(ycord1),ycord2),col='gray60',border = NA)
        graphics::lines(Y[, i])
      }
      graphics::lines(c(rep(NA, t0 - 2), Y[t0 - 1, i], cf[, i]), col = "blue")
      graphics::abline(v = t0, col = "blue", lty = 2)
      
      if(display.fitted==TRUE){
        graphics::lines(fitted[,i],col="red")
      }
    }
    
    if(display.fitted==TRUE){
      graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      graphics::legend("bottom",legend = c("Observed","Fitted","Counterfactual"),col=c(1,2,4),
                       lwd=c(2,2,2),lty=c(1,1,1),bty="n",seg.len = 1,horiz = TRUE, inset = c(0,0),xpd=TRUE,cex=1.2)
    }
  }else{
    if(display.fitted==TRUE){
      graphics::par(oma = c(4, 1, 1, 1))
    }
    graphics::par(mfrow = c(ceiling(length(plot)/ncol), ncol))  
    for(i in 1:length(plot)){
      graphics::plot(Y[, plot[i]], type = "l", ylab = ylab[i], xlab = "Time",main=main[i],ylim=c(y.min[i],y.max[i]),...)
      
      if(confidence.bands==TRUE){
        aux=Reduce("cbind",lapply(boot.cf,function(x)x[,plot[i]]))
        aux1=t(apply(aux,1,sort))
        intervals=c(round(ncol(aux)*(alpha/2)),round(ncol(aux)*(1-alpha/2)))
        xcord=t0:length(Y[,plot[i]])
        ycord1=aux1[,intervals[1]]
        ycord2=aux1[,intervals[2]]
        graphics::polygon(c(rev(xcord),xcord),c(rev(ycord1),ycord2),col='gray60',border = NA)
        graphics::lines(Y[, plot[i]])
      }
      
      graphics::lines(c(rep(NA, t0 - 2), Y[t0 - 1, plot[i]], cf[, plot[i]]), col = "blue")
      graphics::abline(v = t0, col = "blue", lty = 2)
      
      if(display.fitted==TRUE){
        graphics::lines(fitted[,plot[i]],col="red")
      }
    }
    if(display.fitted==TRUE){
      graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      graphics::legend("bottom",legend = c("Observed","Fitted","Counterfactual"),col=c(1,2,4),
                       lwd=c(2,2,2),lty=c(1,1,1),bty="n",seg.len = 1,cex=1.2,horiz = TRUE, inset = c(0,0),xpd=TRUE)
    }
  }
  graphics::par(mar=oldmar,mfrow=oldmfrow,oma=oldoma)
  
}


