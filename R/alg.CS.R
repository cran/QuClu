#' @title CS quantile-based clustering algorithm
#'
#' @description This function allows to run the CS (Common theta and Scaled variables through lambda_j) version of the quantile-based clustering algorithm.
#'
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param k The number of clusters. The default is k=2.
#' @param eps The relative convergence tolerances for objective function. The default is set to 1e-8.
#' @param it.max A number that gives integer limits on the number of the CS algorithm iterations. By default, it is set to 100.
#' @param B The number of times the initialization step is repeated; the default is 30.
#' @param lambda The initial value for lambda_j, the variable scaling parameters. By default, lambdas are set to be equal to 1.
#' @details Algorithm CS: Common theta and Scaled variables via lambda_j. A common value of theta is taken but variables are scaled through lambda_j.
#' @return A list containing the following elements:
#' \item{cl}{A vector whose [i]th entry is classification of observation i in the test data.}
#' \item{qq}{A matrix whose [h,j]th entry is the theta-quantile of variable j in cluster h.}
#' \item{theta}{The estimated common theta.}
#' \item{Vseq}{The values of the objective function V at each step of the algorithm.}
#' \item{V}{The final value of the objective function V.}
#' \item{lambda}{A vector containing the scaling factor for each variable.}
#' @references C. Hennig, C. Viroli, L. Anderlucci (2018). \emph{Quantile-based clustering}. \url{http://arxiv.org/abs/1806.10403}
#' @examples out <- alg.CS(iris[,-5],k=3)
#' out$theta
#' out$qq
#' out$lambda
#'
#' table(out$cl)
#' @export

alg.CS=function(data,k=2,eps=1e-8,it.max=100,B=30,lambda=rep(1,p))
{
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
  for (i in 1:k) {qq[i,]=apply(data,2,stats::quantile,prob=(i-1)/(k-1)*0.5+theta/2)
  QQ[,i]=rowSums((theta+((1-2*theta)*(data<t(matrix(qq[i,],p,numobs)))))*abs(data-t(matrix(qq[i,],p,numobs))))-p*log(theta*(1-theta))}
  cl=apply(QQ,1,which.min)
  VV.temp=sum(QQ[cbind(seq_along(cl),cl)])
  if (VV.temp<VV[1]) { VV[1]=VV.temp
  cl.true=cl
  theta.true=theta}
  }

  cl=cl.true
  theta=theta.true


  ###########################

  ### clustering step
  ratio=5
  h=1
  while ((ratio>eps) & (h<it.max)) {
    h=h+1

    ### a) compute nk
    nk=table(cl)
    if (length(nk)<k) nk=table(factor(cl,levels=1:k))   ## se la classificazione degenera va in errore

    ### b) compute the new barycenters
    for (i in 1:k) {if (nk[i]>0) qq[i,]=apply(data[cl==i,,drop=FALSE],2,stats::quantile,theta)
    QQ0[,,i]=as.matrix((theta+((1-2*theta)*(data<t(matrix(qq[i,],p,numobs)))))*abs(data-t(matrix(qq[i,],p,numobs))))
    QQ[,i]=rowSums(matrix(lambda,numobs,p,byrow=TRUE)*QQ0[,,i])-sum(log(lambda*theta*(1-theta)))}

    ### c) compute theta
    theta=stats::optim(theta,fn.cs,method="L-BFGS-B",lower=0.0001,upper=0.999,data=data,k=k,cl=cl,qq=qq,lambda=lambda)$par


    ### d) estimate lambda
    select<-function(QQ0,cl) return(QQ0[cbind(seq_along(cl),cl)])
    den=colMeans(apply(QQ0,2,select,cl))
    den=ifelse(den==0,eps,den)
    lambda=1/den

    ### e) compute z
    cl=apply(QQ,1,which.min)
    VV[h]=sum(QQ[cbind(seq_along(cl),cl)])

    ratio=(VV[h-1]-VV[h])/VV[h-1]
    if (h<5) ratio=2*eps
  }

  names(lambda)<-colnames(qq)<-colnames(data)

  return(list(Vseq=VV,V=VV[h],cl=cl,qq=qq,theta=theta,lambda=lambda))
}


fn.cs=function(theta,data,k,cl,qq,lambda)
{VV=0
p=ncol(data)
numobs=length(cl)
for (i in 1:k) if (sum(cl==i)>0) {
  nn=sum(cl==i)
  xx=data[cl==i,,drop=FALSE]
  a=rowSums(matrix(lambda,nn,p,byrow=TRUE)*(theta+((1-2*theta)*(xx<t(matrix(qq[i,],p,nn)))))*abs(xx-t(matrix(qq[i,],p,nn))))
  VV=VV +sum(a)-sum(nn*log(lambda*theta*(1-theta)))
}

return(VV)

}
