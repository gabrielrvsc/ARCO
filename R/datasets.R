#' A generated dataset used in the examples
#' 
#' This data contains 100 observations of 20 variables generated using the dgp from Carvalho, Masini and Medeiros (2016). Each variable is one ArCo unit. The intervention took place on the first unit at t0=51 by adding a constant equal to 0.628, which is one standard deviation of the treated unit variable before the intervention.
#'
#' @docType data
#' @keywords datasets
#' @name data.q1
#' @usage data(data.q1)
#' @format A matrix with 100 rows and 20 variables.
#' @references Carvalho, C., Masini, R., Medeiros, M. (2016) "ArCo: An Artificial Counterfactual Approach For High-Dimensional Panel Time-Series Data.".
NULL


#' A dataset used in the examples
#'
#' This data is a list with two matrixes, each have 100 observations and 6 variables. It was generatedthe using the dgp from Carvalho, Masini and Medeiros (2016). In the ArCo context, each matrix is a variable and each variable in the matrixes is an unit. The intervention took place on the first unit at t0=51 by adding constants of 0.840 and 0.511 (one standar deviation before the intervention) on variables (matrixes) 1 and 2. 
#'
#' @docType data
#' @keywords datasets
#' @name data.q2
#' @usage data(data.q2)
#' @format A list with 2 matrixes of 100 rows and 6 variables.
#' @references Carvalho, C., Masini, R., Medeiros, M. (2016) "ArCo: An Artificial Counterfactual Approach For High-Dimensional Panel Time-Series Data.".
NULL

#' Dataset used on the empirical example by Carvalho, Masini and Medeiros (2016).
#' 
#' This is the data from the \emph{nota fiscal paulista} (NFP) example from Carvalho, Masini and Medeiros (2016). The variables are the food away from home component of the inflation and the GDP for 9 metropolitan areas in Brazil. Each variable is represented by a matrix inside the list. The treated unit is the Sao Paulo metropolitan area, which is the first column in each matrix. The treatment took place at \eqn{t_0=34}. 
#'
#' @docType data
#' @keywords datasets
#' @name inflationNFP
#' @usage data(inflationNFP)
#' @format A list with two matrixes of 56 rows and 9 variables.
#' @references Carvalho, C., Masini, R., Medeiros, M. (2016) "ArCo: An Artificial Counterfactual Approach For High-Dimensional Panel Time-Series Data.".
NULL