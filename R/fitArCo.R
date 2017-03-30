#' Estimates the ArCo using the model selected by the user
#' 
#' Estimates the Artificial Counterfactual unsing any model supplied by the user, calculates the most relevant statistics and allows for the counterfactual confidence intervals to be estimated by block bootstrap. \cr
#' The model must be supplied by the user through the arguments fn and p.fn. The first determines which function will be used to estimate the model and the second determines the forecasting function. For more details see the examples and the description on the arguments.
#' 
#' @details This description may be useful to clarify the notation and understand how the arguments must be supplied to the functions.
#' \itemize{
#' \item{units: }{Each unity is indexed by a number between \eqn{1,\dots,n}. They are for exemple: countries, states, municipalities, firms, etc.}
#' \item{Variables: }{For each unity and for every time period \eqn{t=1,\dots,T} we observe \eqn{q_i \ge 1} variables. They are for example: GDP, inflation, sales, etc.}
#' \item{Intervention: }{The intervention took place only in the treated unity at time \eqn{t_0=\lambda_0*T}, where \eqn{\lambda_0} is in (0,1).}
#' }
#' 
#' @param data A list of matrixes or data frames of length q. Each matrix is T X n and it contains observations of a single variable for all units and all periods of time. Even in the case of a single variable (q=1), the matrix must be inside a list.
#' @param fn The function used to estimate the first stage model. This function must receive only two arguments in the following order: X (independent variables), y (dependent variable). If the model requires additional arguments they must be supplied inside the function fn.
#' @param p.fn The forecasting function used to estimate the counterfactual using the first stage model (normally a predict funtion). This function also must receive only two arguments in the following order: model (model estimated in the first stage), newdata (out of sample data to estimate the second stage). If the prediction requires additional arguments they must be supplied inside the function p.fn.
#' @param treated.unity Single number indicating the unity where the intervention took place.
#' @param t0 Single number indicating the intervention period.
#' @param lag Number of lags in the first stage model. Default is 0, i.e. only contemporaneous variables are used.
#' @param Xreg Exogenous controls.
#' @param alpha Significance level for the delta confidence bands.
#' @param boot.cf Should bootstrap confidence intervals for the counterfactual be calculated (default=FALSE). 
#' @param R Number of bootstrap replications in case boot.cf=TRUE.
#' @param l Block length for the block bootstrap.  
#' @param VCOV.type Type of covariance matrix for the delta. "iid" for standard covariance matrix, "var" or "varhac" to use prewhitened covariance matrix using VAR models, "varhac" selects the order of the VAR automaticaly and "nw" for Newey West. In the last case the user may select the kernel type and combine the kernel with the VAR prewhitening. For more details see Andrews and Monahan (1992).
#' @param VCOV.lag Lag used on the robust covariance matrix if VCOV.type is different from "iid".
#' @param bandwidth.kernel Kernel bandwidth. If NULL the bandwidth is automatically calculated.
#' @param kernel.type Kernel to be used for VCOV.type="nw".
#' @param VHAC.max.lag Maximum lag of the VAR in case VCOV.type="varhac".  
#' @param prewhitening.kernel If TRUE and VCOV.type="nw", the covariance matrix is calculated with prewhitening (default=FALSE).
#' @return An object with S3 class fitArCo.
#' \item{cf}{estimated counterfactual}
#' \item{fitted.values}{In sample fitted values for the pre-treatment period.}
#' \item{model}{A list with q estimated models, one for each variable. Each element in the list is the output of the fn function.}
#' \item{delta}{The delta statistics and its confidence interval.}
#' \item{p.value}{ArCo p-value.}
#' \item{data}{The data used.}
#' \item{t0}{The intervention period used.}
#' \item{treated.unity}{The treated unity used.}
#' \item{boot.cf}{A list with the bootstrap result (boot.cf=TRUE) or logical FALSE (boot.cf=FALSE). In the first case, each element in the list refers to one bootstrap replication of the counterfactual, i. e. the list length is R.}
#' \item{call}{The matched call.}
#' @keywords ArCo
#' @export
#' @import Matrix glmnet
#' @importFrom stats cov embed qnorm
#' @examples 
#' #############################
#' ## === Example for q=1 === ##
#' #############################
#' data(data.q1)
#' # = First unity was treated on t=51 by adding 
#' # a constant equal to one standard deviation
#' 
#' data=list(data.q1) # = Even if q=1 the data must be in a list
#' 
#' ## == Fitting the ArCo using linear regression == ##
#' # = creating fn and p.fn function = #
#' fn=function(X,y){
#' return(lm(y~X))
#' }
#' p.fn=function(model,newdata){
#' b=coef(model)
#' return(cbind(1,newdata) %*% b)
#' }
#' 
#' ArCo=fitArCo(data = data,fn = fn, p.fn = p.fn, treated.unity = 1 , t0 = 51)
#' 
#' #############################
#' ## === Example for q=2 === ##
#' #############################
#' 
#' # = First unity was treated on t=51 by adding constants of one standard deviation
#' # for the first and second variables
#' 
#' data(data.q2) # data is already a list
#' 
#' ## == Fitting the ArCo using the package glmnet == ##
#' ## == Quadratic Spectral kernel weights for two lags == ##
#' 
#' ## == Fitting the ArCo using the package glmnet == ##
#' ## == Bartlett kernel weights for two lags == ##
#' require(glmnet)
#' set.seed(123)
#' ArCo2=fitArCo(data = data.q2,fn = cv.glmnet, p.fn = predict,treated.unity = 1 , t0 = 51, 
#'              VCOV.type = "nw",kernel.type = "QuadraticSpectral",VCOV.lag = 2)
#'
#' @references Carvalho, C., Masini, R., Medeiros, M. (2016) "ArCo: An Artificial Counterfactual Approach For High-Dimensional Panel Time-Series Data.".
#' 
#' Andrews, D. W., & Monahan, J. C. (1992). An improved heteroskedasticity and autocorrelation consistent covariance matrix estimator. Econometrica: Journal of the Econometric Society, 953-966.
#' @seealso \code{\link{plot}}, \code{\link{estimate_t0}}, \code{\link{panel_to_ArCo_list}}

fitArCo=function (data, fn, p.fn, treated.unity, t0, lag = 0, Xreg = NULL, alpha = 0.05, boot.cf = FALSE, R = 100, l = 3,VCOV.type=c("iid","var","nw","varhac"),VCOV.lag=1,bandwidth.kernel=NULL,kernel.type=c("QuadraticSpectral","Truncated","Bartlett","Parzen","TukeyHanning"),VHAC.max.lag=5,prewhitening.kernel=FALSE) 
{
  VCOV.type=match.arg(VCOV.type)
  kernel.type=match.arg(kernel.type)
  if (boot.cf == TRUE) {
    if (R < 10) {
      stop("Minimum number of bootstrap samples is 10.")
    }
  }
  if (is.null(names(data))) {
    names(data) = paste("Variable", 1:length(data), sep = "")
  }
  for (i in 1:length(data)) {
    if (is.null(colnames(data[[i]]))) {
      colnames(data[[i]]) = paste("Unity", 1:ncol(data[[i]]), 
                                  sep = "")
    }
  }
  for (i in 1:length(data)) {
    aux = length(unique(colnames(data[[i]])))
    k = ncol(data[[i]])
    if (aux < k) {
      colnames(data[[i]]) = paste("Unity", 1:ncol(data[[i]]), 
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
    aux = list()
    for (i in 1:length(data)) {
      aux[[i]] = paste(names(data)[i], colnames(data[[i]])[-treated.unity], 
                       sep = ".")
    }
    colnames(X) = unlist(aux)
  }
  Y.raw = Y
  if (lag != 0) {
    aux1 = sort(rep(0:lag, ncol(X)))
    aux = paste(rep(colnames(X), lag + 1), "lag", aux1, sep = ".")
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
  T=nrow(X)
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
  boot.list = FALSE
  if (boot.cf == TRUE) {
    serie = cbind(y.fit, x.fit)
    q = length(data)
    bootfunc = function(serie) {
      y.fit = serie[, 1:q]
      x.fit = serie[, -c(1:q)]
      if (is.vector(y.fit)) {
        y.fit = matrix(y.fit, ncol = 1)
      }
      save.cf.boot = matrix(NA, nrow(x.pred), q)
      for (i in 1:q) {
        model.boot = fn(x.fit, y.fit[, i])
        contra.fact.boot = p.fn(model.boot, x.pred)
        save.cf.boot[, i] = contra.fact.boot
      }
      return(as.vector(save.cf.boot))
    }
    boot.cf = boot::tsboot(serie, bootfunc, R = R, l = 3, 
                           sim = "fixed")
    boot.stat = boot.cf$t
    boot.list = list()
    for (i in 1:nrow(boot.stat)) {
      boot.list[[i]] = matrix(boot.stat[i, ], ncol = q)
    }
  }
  delta.aux = tail(Y.raw, nrow(save.cf)) - save.cf
  delta = colMeans(delta.aux)
  aux = matrix(0, T, length(data))
  aux[(t0 - lag):nrow(aux), ] = 1
  vhat = Y - (save.fitted + t(t(aux) * delta))
  v1 = matrix(vhat[1:(t0 - lag - 1), ], ncol = length(data))
  v2 = matrix(vhat[(t0 - lag):nrow(vhat), ], ncol = length(data))
  
  t0lag=t0-lag
  sigmahat=T*switch(VCOV.type,
                    iid = cov(v1)/(t0lag-1) + cov(v2)/(T-t0lag),
                    var = VAR(v1,VCOV.lag)$LR/(t0lag-1) + VAR(v2,VCOV.lag)$LR/(T-t0lag),
                    nw  = neweywest(v1,NULL,kernel.type,prewhitening.kernel,VCOV.lag)/(t0lag-1) + neweywest(v2,NULL,kernel.type,prewhitening.kernel,VCOV.lag)/(T-t0lag),
                    varhac = VARHAC(v1,VHAC.max.lag)/(t0lag-1) + VARHAC(v2,VHAC.max.lag)/(T-t0lag)
  )
  
  
  
  w = sqrt(diag(sigmahat))
  W = T * t(delta) %*% solve(sigmahat) %*% delta
  p.value = 1 - stats::pchisq(W, length(delta))
  
  uI = delta + (w * qnorm(1 - alpha/2))/sqrt(T)
  lI = delta - (w * qnorm(1 - alpha/2))/sqrt(T)
  delta.stat = cbind(LB = lI, delta = delta, UB = uI)
  names(model.list) = names(data)
  colnames(save.cf) = names(data)
  rownames(save.cf) = tail(rownames(Y.raw), nrow(save.cf))
  colnames(save.fitted) = names(data)
  rownames(save.fitted) = head(rownames(Y), nrow(save.fitted))
  rownames(delta.stat) = names(data)
  save.fitted = head(save.fitted, nrow(save.fitted) - nrow(save.cf))
  if (typeof(boot.list) == "list") {
    NAboot = Reduce(sum, boot.list)
    if (is.na(NAboot)) {
      warning("Some of the boostrap counterfactuals may have returned NA values. \n \n              A possible cause is the number of observations being close the number of variables if the lm function was used.")
    }
  }
  result = list(cf = save.cf, fitted.values = save.fitted, model = model.list, 
                delta = delta.stat, p.value = p.value, data = data, t0 = t0, 
                treated.unity = treated.unity, boot.cf = boot.list, call = match.call())
  class(result) = "fitArCo"
  return(result)
}