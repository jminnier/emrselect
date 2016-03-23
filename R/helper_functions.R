## ==================================================================================== ##
# emrselect R package for Automated Feature Selection with Electronic Medical Record Data
# Copyright (C) 2015-2016  Jessica Minnier
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# You may contact the author of this code, Jessica Minnier, at <minnier@ohsu.edu>
## ==================================================================================== ##

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

#' Sum.I
#' @export
sum.I <-function(yy,FUN,Yi,Vi=NULL,ties.method = "first") ## ties either 'f' or 'a'
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi),ties.method=ties.method)[1:length(yy)]-rank(yy,ties.method=ties.method)
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

