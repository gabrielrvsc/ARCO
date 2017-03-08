#' Estimates intervention time on a given treated unity
#' 
#' This description may be useful to clarify the notation and understand how the arguments must be supplied to the functions. \cr \cr
#' * units: Each unity is indexed by a number between 1,...,n. They are for exemple: countries, states, municipalities, firms, etc. \cr \cr 
#' * Variables:  For each unity and for every time period t=1,...,T we observe q_i >= 1 variables. They are for example: GDP, inflation, sales, etc.\cr \cr
#' * Intervention:  The intervention took place only in the treated unity at time t0=L0*T, where L0 is in (0,1).
#' 
#' @inheritParams fitArCo
#' @param start Initial value of L0 to be tested.
#' @param end Final value of L0 to be tested.
#' @export
#' @import Matrix glmnet randomForest
#' @examples 
#' #############################
#' ## === Example for q=1 === ##
#' #############################
#' data(data.q1) 
#' # = First unity was treated on t=51 by adding a constant equal 3
#' 
#' data=list(data.q1) # = Even if q=1 the data must be in a list
#' 
#' ## == Fitting the ArCo using linear regression == ##
#' 
#' # = creating fn and p.fn function = #
#' fn=function(X,y){
#'     return(lm(y~X))
#' }
#' p.fn=function(model,newdata){
#'     b=coef(model)
#'     return(cbind(1,newdata)%*%b)
#' }
#' 
#' t0a=estimate_t0(data = data,fn = fn, p.fn = p.fn, treated.unity = 1 )
#' 
#' 
#' #############################
#' ## === Example for q=2 === ##
#' #############################
#' 
#' # = First unity was treated on t=51 by adding constants 15 and -10 
#' # for the first and second variables
#' data(data.q2) # data is already a list
#' 
#' ## == Detecting lambda0 using the package glmnet via LASSO and crossvalidation == ##
#' 
#' t0b=estimate_t0(data = data.q2,fn = fn, p.fn = p.fn, treated.unity = 1, start=0.4)



estimate_t0=function (data, fn, p.fn, start = 0.3, end = 0.95, treated.unity = 1, 
          lag = 0, Xreg = NULL) 
{
  for (i in 1:length(data)) {
    if (is.null(colnames(data[[i]]))) {
      colnames(data[[i]]) = paste("V", i, "-U", 1:ncol(data[[i]]), 
                                  sep = "")
    }
  }
  for (i in 1:length(data)) {
    aux = length(unique(colnames(data[[i]])))
    k = ncol(data[[i]])
    if (aux < k) {
      colnames(data[[i]]) = paste("V", i, "-U", 1:ncol(data[[i]]), 
                                  sep = "")
    }
  }
  T = nrow(data[[1]])
  starting.point = floor(start * T)
  ending.point = min(floor(end * T),nrow(data[[1]])-lag+1)
  save.delta = matrix(0, T, length(data))
  for (i in starting.point:ending.point) {
    m = fitArCo(data = data, fn = fn, p.fn = p.fn, treated.unity = treated.unity, 
                lag = lag, t0 = i, Xreg = Xreg)
    delta = m$delta[,2]
    save.delta[i, ] = delta
  }
  delta.norm = sqrt(rowSums((save.delta)^2))
  t0 = which(delta.norm == max(delta.norm))
  delta = delta.norm[t0]
  return(c(to = t0, delta.norm = delta))
}
