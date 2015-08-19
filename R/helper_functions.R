
#' Vector to matrix
#'
#' @param vc vector
#' @param dm number of rows
#'
#' @return matrix
#' @export
#'
#' @examples VTM(1:10,3)
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

#' logit function
#' @export
logit = function(xx){log(xx/(1-xx))}

#' g inverse logit function
#' @export
g = function(xx){exp(xx)/(1+exp(xx))}

#' d1g
#' @export
d1g = function(xx){g(xx)*(1-g(xx))}

#' d2g
#' @export
d2g = function(xx){d1g(xx)*(1-2*g(xx))}

#' g.logit same as g
#' @export
g.logit = function(xx){exp(xx)/(exp(xx)+1)}

#' dg.logit
#' @export
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}

#' Intercept of Logistic Regression
#' @export
Intercept.GLM = function(yi,bxi){glm(yi~bxi,family="binomial")}

#' ProbTestPos
#' @export
ProbTestPos <- function(tpr,fpr,prev)
{
  tpr*prev + fpr*(1-prev)
}

#' L2Norm
#' @export
L2Norm <- function(data,coef,intercept=F)
{
  yy = data[,1]; xx = data[,-1]; if(intercept){xx=cbind(1,xx)}
  apply((yy - xx%*%t(coef))^2,2,mean)
}


#' S.FUN
#' @export
S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

#' Sinv.FUN
#' @export
Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth)
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

#' Create autocorrelation matrix
#' @export
autocorr.mat <- function(p = 100, rho = 0.9) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}


