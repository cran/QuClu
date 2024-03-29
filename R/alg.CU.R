#' @title CU quantile-based clustering algorithm
#'
#' @description This function allows to run the CU (Common theta and Unscaled variables) version of the quantile-based clustering algorithm.
#'
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param k The number of clusters. The default is k=2.
#' @param eps The relative convergence tolerances for objective function. The default is set to 1e-8.
#' @param it.max A number that gives integer limits on the number of the CU algorithm iterations. By default, it is set to 100.
#' @param B The number of times the initialization step is repeated; the default is 30.
#' @details Algorithm CU: Common theta and Unscaled variables. A common value of theta for all the variables is assumed.
#' This strategy directly generalizes the conventional k-means to other moments of the distribution to better accommodate skewness in the data.
#' @return A list containing the following elements:
#' \item{method}{The chosen parameterization, CU, Common theta and Unscaled variables}
#' \item{k}{The number of clusters.}
#' \item{cl}{A vector whose [i]th entry is classification of observation i in the test data.}
#' \item{qq}{A matrix whose [h,j]th entry is the theta-quantile of variable j in cluster h.}
#' \item{theta}{A vector whose [j]th entry is the percentile theta for variable j.}
#' \item{Vseq}{The values of the objective function V at each step of the algorithm.}
#' \item{V}{The final value of the objective function V.}
#' \item{lambda}{A vector containing the scaling factor for each variable.}
#' @references Hennig, C., Viroli, C., Anderlucci, L. (2019) "Quantile-based clustering" \emph{Electronic Journal of Statistics}, 13 (2) 4849-4883  <doi:10.1214/19-EJS1640>
#' @examples out <- alg.CU(iris[,-5],k=3)
#' out$theta
#' out$qq
#'
#' table(out$cl)
#' @export

alg.CU=function(data,k=2,eps=1e-8,it.max=100,B=30)
{
  data<-as.matrix(data)
  numobs=nrow(data)
  p=ncol(data)

  qq=matrix(0,k,p)
  QQ=matrix(0,numobs,k)
  QQ0=array(0,c(numobs,p,k))
  VV.temp=NULL
  VV=0

  #######################################################################
  ## initialization
  #######################################################################

  ## 1) equispaced interval
  VV[1]=Inf

  for (hh in 1:B) {theta=stats::runif(1)
  for (i in 1:k) {
    if (k==1) qq[i,]=apply(data,2,stats::quantile,prob=theta/2)
    else qq[i,]=apply(data,2,stats::quantile,prob=(i-1)/(k-1)*0.5+theta/2)
  QQ[,i]=rowSums((theta+((1-2*theta)*(data<t(matrix(qq[i,],p,numobs)))))*abs(data-t(matrix(qq[i,],p,numobs))))-p*log(theta*(1-theta))}
  cl=apply(QQ,1,which.min)
  VV.temp=sum(QQ[cbind(seq_along(cl),cl)])
  if (VV.temp<VV[1]) { VV[1]=VV.temp
  cl.true=cl
  theta.true=theta}
  }

  cl=cl.true
  theta=theta.true

  ### clustering step
  ratio=5
  h=1
  while ((ratio>eps) & (h<it.max)) {h=h+1

  ### a) compute nk
  nk=table(cl)
  if (length(nk)<k) nk=table(factor(cl,levels=1:k))   ## se la classificazione degenera va in errore

  ### b) compute the new barycenters
  for (i in 1:k) {if (nk[i]>0) qq[i,]=apply(data[cl==i,,drop=FALSE],2,stats::quantile,theta)
  QQ[,i]=rowSums((theta+((1-2*theta)*(data<t(matrix(qq[i,],p,numobs)))))*abs(data-t(matrix(qq[i,],p,numobs))))-p*log(theta*(1-theta))}

  ### c) compute theta
  theta=stats::optim(theta,fn.cu,method="L-BFGS-B",lower=0.0001,upper=0.999,data=data,k=k,cl=cl,qq=qq)$par

  ### d) compute z
  cl=apply(QQ,1,which.min)
  VV[h]=sum(QQ[cbind(seq_along(cl),cl)])

  ratio=(VV[h-1]-VV[h])/VV[h-1]
  if (h<5) ratio=2*eps
  }

  ### d) estimate lambda
  for (i in 1:k)  QQ0[,,i]=as.matrix((theta+((1-2*theta)*(data<t(matrix(qq[i,],p,numobs)))))*abs(data-t(matrix(qq[i,],p,numobs))))

  select<-function(QQ0,cl) return(QQ0[cbind(seq_along(cl),cl)])
  den=colMeans(apply(QQ0,2,select,cl))
  den=ifelse(den==0,eps,den)
  lambda=1/den

  theta<-rep(theta,p)
  names(theta)<-names(lambda)<-colnames(qq) <- colnames(data)
  return(list(method="CU",k=k,Vseq=VV,V=VV[h],cl=cl,qq=qq,theta=theta,lambda=lambda))
}



fn.cu=function(theta,data,k,cl,qq)
{VV=0
p=ncol(data)
for (i in 1:k) if (sum(cl==i)>0) {
  nn=sum(cl==i)
  xx=data[cl==i,,drop=FALSE]
  a=rowSums((theta+((1-2*theta)*(xx<t(matrix(qq[i,],p,nn)))))*abs(xx-t(matrix(qq[i,],p,nn))))
  VV=VV +sum(a)
}
numobs=length(cl)
VV=VV-p*numobs*log(theta*(1-theta))
return(VV)
}

