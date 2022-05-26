#' @title Quantile-based clustering algorithm
#'
#' @description This function allows to run the $k$-quantile clustering algorithm, allowing for different constraints: common theta and unscaled variables (CU), common theta and scaled variables (CS), variable-wise theta and unscaled variables (VU) and the variable-wise theta and scaled variables (VS).
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param k The number of clusters. The default is k=2.
#' @param method The chosen constrained method. The options are: CU (Common theta and Unscaled variables), CS (Common theta and Scaled variables), VU (Variable-wise theta and Unscaled variables), VS (Variable-wise theta and Scaled variables).The default is the unconstrained method, VS.
#' @param eps The relative convergence tolerances for objective function. The default is set to 1e-8.
#' @param it.max A number that gives integer limits on the number of the algorithm iterations. By default, it is set to 100.
#' @param B The number of times the initialization step is repeated; the default is 30.
#' @param lambda The initial value for lambda_j, the variable scaling parameters, for models CS and VS. By default, lambdas are set to be equal to 1.
#' @details Algorithm CU: Common theta and Unscaled variables. A common value of theta for all the variables is assumed. Algorithm CS: Common theta and Scaled variables via lambda_j. A common value of theta is taken but variables are scaled through lambda_j. Algorithm VU: Variable-wise theta_j and Unscaled variables. A different theta for every single variable is estimated to better accomodate different degree of skeweness in the data. Algorithm VS: Variable-wise theta_j and Scaled variables via lambda_j. A different theta for every single variable is estimated to better accomodate different degree of skeweness in the data and variables are scaled through lambda_j.
#' @return A list containing the following elements:
#' \item{method}{The chosen parameterization.}
#' \item{k}{The number of clusters.}
#' \item{cl}{A vector whose [i]th entry is classification of observation i in the test data.}
#' \item{qq}{A matrix whose [h,j]th entry is the theta-quantile of variable j in cluster h.}
#' \item{theta}{A vector whose [j]th entry is the percentile theta for variable j.}
#' \item{Vseq}{The values of the objective function V at each step of the algorithm.}
#' \item{V}{The final value of the objective function V.}
#' \item{lambda}{A vector containing the scaling factor for each variable.}
#' @references Hennig, C., Viroli, C., Anderlucci, L. (2019) "Quantile-based clustering" \emph{Electronic Journal of Statistics}, 13 (2) 4849-4883  <doi:10.1214/19-EJS1640>
#' @examples out <- kquantiles(iris[,-5],k=3,method="VS")
#' out$theta
#' out$qq
#'
#' table(out$cl)
#' @export

kquantiles<-function(data,k=2,method="VS",eps=1e-8,it.max=100,B=30,lambda=NULL){
  p=ncol(data)
  if (!(method %in% c("VS","VU","CS","CU"))) {
    stop("Error: method must be either 'CU', 'CS', 'VU', 'VS'.")
  }

  if (method=="VS"){
    if (is.null(lambda)) lambda<-rep(1,p)
    out<-alg.VS(data,k,eps,it.max,B,lambda)
  }
  if (method=="VU"){
    out<-alg.VU(data,k,eps,it.max,B)
  }
  if (method=="CS"){
    if (is.null(lambda)) lambda<-rep(1,p)
    out<-alg.CS(data,k,eps,it.max,B,lambda)
  }

  if (method=="CU"){
    out<-alg.CU(data,k,eps,it.max,B)
  }

  return(out)
  }

