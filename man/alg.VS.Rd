% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alg.VS.R
\name{alg.VS}
\alias{alg.VS}
\title{VS quantile-based clustering algorithm}
\usage{
alg.VS(data, k = 2, eps = 1e-08, it.max = 100, B = 30, lambda = rep(1, p))
}
\arguments{
\item{data}{A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.}

\item{k}{The number of clusters. The default is k=2.}

\item{eps}{The relative convergence tolerances for objective function. The default is set to 1e-8.}

\item{it.max}{A number that gives integer limits on the number of the VS algorithm iterations. By default, it is set to 100.}

\item{B}{The number of times the initialization step is repeated; the default is 30.}

\item{lambda}{The initial value for lambda_j, the variable scaling parameters. By default, lambdas are set to be equal to 1.}
}
\value{
A list containing the following elements:
\item{method}{The chosen parameterization, VS, Variable-wise theta_j and Scaled variables}
\item{k}{The number of clusters.}
\item{cl}{A vector whose [i]th entry is classification of observation i in the test data.}
\item{qq}{A matrix whose [h,j]th entry is the theta-quantile of variable j in cluster h.}
\item{theta}{A vector whose [j]th entry is the percentile theta for variable j.}
\item{Vseq}{The values of the objective function V at each step of the algorithm.}
\item{V}{The final value of the objective function V.}
\item{lambda}{A vector containing the scaling factor for each variable.}
}
\description{
This function allows to run the VS (Variable-wise theta_j and Scaled variables through lambda_j) version of the quantile-based clustering algorithm.
}
\details{
Algorithm VS: Variable-wise theta_j and Scaled variables via lambda_j. A different theta for every single variable is estimated to better accomodate different degree of skeweness in the data and variables are scaled through lambda_j.
}
\examples{
out <- alg.VS(iris[,-5],k=3)
out$theta
out$qq
out$lambda

table(out$cl)
}
\references{
Hennig, C., Viroli, C., Anderlucci, L. (2019) "Quantile-based clustering" \emph{Electronic Journal of Statistics}, 13 (2) 4849-4883  <doi:10.1214/19-EJS1640>
}
