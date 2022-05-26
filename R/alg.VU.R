#' @title VU quantile-based clustering algorithm
#'
#' @description This function allows to run the VU (Variable-wise theta_j and Unscaled variables) version of the quantile-based clustering algorithm.
#' @param data A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations and columns correspond to variables.
#' @param k The number of clusters. The default is k=2.
#' @param eps The relative convergence tolerances for objective function. The default is set to 1e-8.
#' @param it.max A number that gives integer limits on the number of the VU algorithm iterations. By default, it is set to 100.
#' @param B The number of times the initialization step is repeated; the default is 30.
#' @details Algorithm VU: Variable-wise theta_j and Unscaled variables. A different theta for every single variable is estimated to better accomodate different degree of skeweness in the data.
#' @return A list containing the following elements:
#' \item{method}{The chosen parameterization, VU, Variable-wise theta_j and Unscaled variables}
#' \item{k}{The number of clusters.}
#' \item{cl}{A vector whose [i]th entry is classification of observation i in the test data.}
#' \item{qq}{A matrix whose [h,j]th entry is the theta-quantile of variable j in cluster h.}
#' \item{theta}{A vector whose [j]th entry is the percentile theta for variable j.}
#' \item{Vseq}{The values of the objective function V at each step of the algorithm.}
#' \item{V}{The final value of the objective function V.}
#' \item{lambda}{A vector containing the scaling factor for each variable.}
#' @references Hennig, C., Viroli, C., Anderlucci, L. (2019) "Quantile-based clustering" \emph{Electronic Journal of Statistics}, 13 (2) 4849-4883  <doi:10.1214/19-EJS1640>
#' @examples out <- alg.VU(iris[,-5],k=3)
#' out$theta
#' out$qq
#'
#' table(out$cl)
#' @export

alg.VU=function(data,k=2,eps=1e-8,it.max=100,B=30)
{
  data<-as.matrix(data)
  numobs=nrow(data)
  p=ncol(data)
  ### init partition e quantiles
  qq=matrix(0,k,p)
  VV=0
  VV.temp=NULL

  QQ=array(0,c(numobs,k,p))
  QQ0=array(0,c(numobs,p,k))

  #######################################################################
  ## initialization
  #######################################################################

  ## 1) equispaced interval
  VV[1]=Inf

  for (hh in 1:B) {theta=stats::runif(p)
  for (j in 1:p) for (i in 1:k) {
    if (k==1) qq[i,j]=stats::quantile(data[,j],theta[j]/2)
    else qq[i,j]=stats::quantile(data[,j],prob=(i-1)/(k-1)*0.5+theta[j]/2)
  QQ[,i,j]=(theta[j]+(1-2*theta[j])*(data[,j] < matrix(qq[i,j],numobs)))*abs(data[,j]-matrix(qq[i,j],numobs))-log(theta[j]*(1-theta[j]))}
  cl=apply(apply(QQ,c(1,2),sum),1,which.min)
  conta=0
  for (j in 1:p) conta=conta+sum(QQ[cbind(seq_along(cl),cl,j)])
  VV.temp=conta
  if (VV.temp<VV[1]) { VV[1]=VV.temp
  cl.true=cl
  theta.true=theta}
  }

  cl=cl.true
  theta=theta.true

  ### clustering step
  ratio=5
  h=1

  while ((ratio>eps) & (h<it.max)) {
    h=h+1

    ### a) compute nk
    nk=table(cl)
    if (length(nk)<k) nk=table(factor(cl,levels=1:k))   ## se la classificazione degenera va in errore

    ### b) compute the new barycenters
    for (j in 1:p) for (i in 1:k) {if (nk[i]>0) qq[i,j]=stats::quantile(data[cl==i,j],theta[j])
    QQ[,i,j]=(theta[j]+(1-2*theta[j])*(data[,j] < matrix(qq[i,j],numobs)))*abs(data[,j]-matrix(qq[i,j],numobs))-log(theta[j]*(1-theta[j]))}

    ### c) compute theta
    for (j in 1:p) theta[j]=stats::optim(theta[j],fn.vu,method="L-BFGS-B",lower=0.0001,upper=0.999,data=data[,j,drop=FALSE],k=k,cl=cl,qq=qq[,j,drop=FALSE])$par

    ### d) compute z
    cl=apply(apply(QQ,c(1,2),sum),1,which.min)

    conta=0
    for (j in 1:p) conta=conta+sum(QQ[cbind(seq_along(cl),cl,j)])
    VV[h]=conta

    ratio=(VV[h-1]-VV[h])/VV[h-1]
    if (h<5) ratio=2*eps
  }

  ### d) estimate lambda
  for (j in 1:p) for (i in 1:k) QQ0[,j,i]=(theta[j]+(1-2*theta[j])*(data[,j] < matrix(qq[i,j],numobs)))*abs(data[,j]-matrix(qq[i,j],numobs))-log(theta[j]*(1-theta[j]))
  select<-function(QQ0,cl) return(QQ0[cbind(seq_along(cl),cl)])
  den=colMeans(apply(QQ0,2,select,cl))
  den=ifelse(den==0,eps,den)
  lambda=1/den

  names(lambda)<-names(theta)<-colnames(qq)<-colnames(data)
  return(list(method="VU",k=k,Vseq=VV,V=VV[h],cl=cl,qq=qq,theta=theta,lambda=lambda))
}

fn.vu=function(theta,data,k,cl,qq)
{VV=0
p=ncol(data)
for (i in 1:k) if (sum(cl==i)>0) {
  nn=sum(cl==i)
  xx=data[cl==i,,drop=FALSE]
  a=rowSums((theta+((1-2*theta)*(xx<t(matrix(qq[i,],p,nn)))))*abs(xx-t(matrix(qq[i,],p,nn))))
  VV=VV +sum(a)
}
numobs=length(cl)
VV=VV-numobs*log(theta*(1-theta))
return(VV)

}

