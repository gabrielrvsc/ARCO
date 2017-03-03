#' Estimates the ARCO using the model selected by the user
#' 
#'   This description may be useful to clarify the notation and understand how the arguments must be supplied to the functions. 
#' * units: Each unity is indexed by a number between 1,...,n. They are for exemple: countries, states, municipalities, firms, etc. 
#' * Variables:  For each unity and for every time period t=1,...,T we observe q_i >= 1 variables. They are for example: GDP, inflation, sales, etc.
#' * Intervention:  The intervention took place only in the treated unity at time t0=L0*T, where L0 is in (0,1).
#' @param data A list of matrixes or dataframes of length q. Each matrix is T X n and it contains observations of a single variable for all units and all periods of time. Even in the case of a single variable (q=1), the matrix must be inside a list.
#' @param fn The function used to estimate the first stage model. This function must receive only two arguments in the following order: X (independent variables), y (dependent variable). If the model requires additional arguments they must be supplied inside the function fn.
#' @param p.fn The function used to estimate the predict using the first stage model. This function also must receive only two arguments in the following order: model (model estimated in the first stage), newdata (out of sample data to estimate the second stage). If the prediction requires additional arguments they must be supplied inside the function p.fn.
#' @param treated.unity Single number indicating the unity where the intervention took place.
#' @param t0 Single number indicating the intervention period.
#' @param lag Number of lags in the first stage model. Default is 0, i.e. only contemporaneous variables are used.
#' @param Xreg Exogenous controls.
#' @param display Default is TRUE. Shows the counterfactual plot automaticaly.
#' @param HACweights Vector of weights for the robust covariance matrix of the delta statistics. Default is 1 for the lag 0 and 0 for all other lags.
#' @param alpha Significance level for the delta.
#' @param ... Aditional parameters used only if display = TRUE.
#' @keywords ARCO
#' @export
#' @import Matrix glmnet randomForest
#' @importFrom graphics abline lines par plot
#' @importFrom stats cov embed qnorm
#' @examples 
#' #############################
   ## === Example for q=1 === ##
   #############################
#' data(data.q1) # = First unity was treated on t=51 by adding a constant equal 3
#' data=list(data.q1) # = Even if q=1 the data must be in a list
#' ## == Fitting the ARCO using linear regression == ##
#' # = creating fn and p.fn function = #
#' fn=function(X,y){
#' return(lm(y~X))
#' }
#' p.fn=function(model,newdata){
#' b=coef(model)
#' return(cbind(1,newdata) %*% b)}
#' ARCO=fitARCO(data = data,fn = fn, p.fn = p.fn, treated.unity = 1 , t0 = 51)
#' 
#' #############################
   ## === Example for q=2 === ##
   #############################
#' 
#' # = First unity was treated on t=51 by adding constants 3 and -3 for the first and second variables
#' data(data.q2) # data is already a list
#' ## == Fitting the ARCO using the package randomForest == ##
#' require(randomForest)
#' ## == Bartlett kernel weights for two lags == ##
#' l=2
#' w <- seq(1, 0, by = -(1/(l + 1)))[1:(l+1)]
#' ARCO2=fitARCO(data = data.q2,fn = randomForest, p.fn = predict,
#' treated.unity = 1 , t0 = 51, HACweights = w)
#' ## == Fitting the ARCO using the package glmnet via LASSO and crossvalidation == ##
#' require(glmnet)
#' ## == Bartlett kernel weights for two lags == ##
#' ARCO3=fitARCO(data = data.q2,fn = cv.glmnet, p.fn = predict, 
#' treated.unity = 1 , t0 = 51, HACweights = w)
#'
#' @references Carvalho, C., Masini, R., Medeiros, M. (2016) "ARCO: An Artificial Counterfactual Approach For High-Dimensional Panel Time-Series Data.".





fitArCo=function (data, fn, p.fn, treated.unity, t0, lag = 0, Xreg = NULL, HACweights = 1, alpha = 0.05) 
{
  if (is.null(names(data))) {
    names(data) = paste("Variable", 1:length(data), sep = "")
  }
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
  if (length(data) == 1) {
    Y = matrix(data[[1]][, treated.unity], ncol = 1)
    X = data[[1]][, -treated.unity]
    X = as.matrix(X)
    colnames(X) = paste(names(data), colnames(data[[1]])[-treated.unity], 
                        sep = ".")
  }else {
    Y = Reduce("cbind", lapply(data, function(x) x[, treated.unity]))
    X = Reduce("cbind", lapply(data, function(x) x[, -treated.unity]))
    aux = rep(NA, length(data))
    for (i in 1:length(data)) {
      aux[i] = colnames(data[[i]])[-treated.unity]
    }
    colnames(X) = paste(aux, names(data), sep = ".")
  }
  Y.raw = Y
  if (lag != 0) {
    aux1 = sort(rep(0:lag, ncol(X)))
    aux = paste(rep(colnames(X), ncol(X)), "lag", aux1, sep = ".")
    X = embed(X, lag + 1)
    colnames(X) = aux
    Y = tail(Y, nrow(X))
  }
  if (length(Xreg) != 0) {
    X = cbind(X, tail(Xreg, nrow(X)))
  }
  if (is.vector(Y)) {
    Y = matrix(Y, length(Y), 1)
  }
  y.fit = matrix(Y[1:(t0 - 1 - lag), ], ncol = length(data))
  y.pred = matrix(Y[-c(1:(t0 - 1 - lag)), ], ncol = length(data))
  x.fit = X[1:(t0 - 1 - lag), ]
  x.pred = X[-c(1:(t0 - 1 - lag)), ]
  save.cf = matrix(NA, nrow(y.pred), length(data))
  save.fitted = matrix(NA, nrow(Y), length(data))
  model.list = list()
  for (i in 1:length(data)) {
    model = fn(x.fit, y.fit[, i])
    model.list[[i]] = model
    contra.fact = p.fn(model, x.pred)
    save.cf[, i] = contra.fact
    save.fitted[, i] = p.fn(model, X)
  }
  delta.aux = tail(Y.raw, nrow(save.cf)) - save.cf
  delta = colMeans(delta.aux)
  aux = matrix(0, nrow(X), length(data))
  aux[(t0 - lag):nrow(aux), ] = 1
  vhat = Y - (save.fitted + t(t(aux) * delta))
  v1 = matrix(vhat[1:(t0 - lag - 1), ], ncol = length(data))
  v2 = matrix(vhat[(t0 - lag):nrow(vhat), ], ncol = length(data))
  tau01 = cov(v1) * HACweights[1]
  tau02 = cov(v2) * HACweights[1]
  tauk1bart = 0
  tauk2bart = 0
  M = length(HACweights) - 1
  if (M >= nrow(v2)) {
    stop("HAC lags bigger than treatment.")
  }
  if (M > 0) {
    for (k in 1:M) {
      aux = cov(v1[(1 + k):nrow(v1), ], v1[1:(nrow(v1) - 
                                                k), ])
      tauk1bart = tauk1bart + (aux + t(aux)) * HACweights[k + 
                                                            1]
      aux = cov(v2[(1 + k):nrow(v2), ], v2[1:(nrow(v2) - 
                                                k), ])
      tauk2bart = tauk2bart + (aux + t(aux)) * HACweights[k + 
                                                            1]
    }
  }
  tauT1 = tau01 + tauk1bart
  tauT2 = tau02 + tauk2bart
  sigmahat = (tauT1/(t0 - 1) + tauT2/(nrow(X) - t0 + 1)) * 
    nrow(X)
  w = sqrt(diag(sigmahat))
  W = nrow(X)*t(delta)%*%solve(sigmahat)%*%delta
  p.value=1-pchisq(W,length(delta))
  uI = delta + (w * qnorm(1 - alpha/2))/sqrt(nrow(X))
  lI = delta - (w * qnorm(1 - alpha/2))/sqrt(nrow(X))
  delta.stat = cbind(LB = lI, delta = delta, UB = uI)
  
  names(model.list) = names(data)
  colnames(save.cf) = names(data)
  rownames(save.cf) = tail(rownames(Y.raw), nrow(save.cf))
  colnames(save.fitted) = names(data)
  rownames(save.fitted) = head(rownames(Y), nrow(save.fitted))
  rownames(delta.stat) = names(data)
  save.fitted=head(save.fitted,nrow(save.fitted)-nrow(save.cf))
  
  return(list(cf = save.cf, fitted=save.fitted, model = model.list, delta = delta.stat,p.value=p.value , data=data, t0=t0, treated.unity=treated.unity))
}
  
