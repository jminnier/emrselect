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

#' Automated feature selection
#'
#' @param dat.X data.frame or matrix of covariates
#' @param dat.S data.frame or matrix of surrogate markers
#' @param b0 number of perturbation resampling draws
#' @param sub.n
#' @param Scov.diag if TRUE restrict covariance of S mclust to
#' be diagonal, default FALSE
#'
#' @return list
#' @export
#'
emrselect <- function(dat.X,dat.S,b0=100,sub.n=NULL,Scov.diag=FALSE) {
  dat.X = as.matrix(dat.X)
  A0 = ncol(as.matrix(dat.S))
  mclust.modelNames=NULL
  if(Scov.diag) {
    #mclust.modelNames = c("EII","VII","EEI","VEI","EVI","VVI")
    mclust.modelNames="EEI"
    }
  ## ========================================================================== ##
  ## multiple surroage analysis: clustering of the multivariate surrogate first ##
  ##                             then use predicted prob as pseudo outcome      ##
  ## ========================================================================== ##
  tmpout_all <- kern_varselect(dat.S=dat.S,dat.X=dat.X,b0=b0,
                              sub.n=sub.n,
                              mclust.modelNames=mclust.modelNames)
  ## ==================================================================================== ##
  ##                     compare to marginal single surrogate analysis                     ##
  ## **1: clustering of a univariate surrogate; then use predicted prob as pseudo outcome ##
  ## **2: threshold each of the univariate surrogate to as a pseudo outcome
  ## ==================================================================================== ##
  bptb1.list = bptb2.list = as.list(1:A0); phat1.marg = phat2.marg = bhat1.marg = bhat2.marg = NULL
  ## single Sk at a time, either first do clustring, or threshold by 1
  for(kk in 1:A0){ #A0 number of S variables
    tmpdat.S = dat.S[,kk,drop=F]
    tmpout1 <- kern_varselect(dat.S=tmpdat.S,dat.X=dat.X,b0=b0,
                              sub.n=sub.n,
                              mclust.modelNames=NULL)

    tmpdat.S = 1*(tmpdat.S >=1)
    tmpout2 <- kern_varselect(dat.S=tmpdat.S,dat.X=dat.X,b0=b0,
                             sub.n=sub.n,
                             mclust.modelNames=NULL)

    bptb1.list[[kk]] = tmpout1$tmpb
    bhat1.marg = cbind(bhat1.marg, tmpout1$bhat); phat1.marg = cbind(phat1.marg,tmpout1$phat)
    bptb2.list[[kk]] = tmpout2$tmpb
    bhat2.marg = cbind(bhat2.marg, tmpout2$bhat); phat2.marg = cbind(phat2.marg,tmpout2$phat)
  }

  return(list(
    "multiple_surrogate_clustering"=list("phat"=tmpout_all$phat,"bhat"=tmpout_all$bhat,"bptb"=tmpout_all$tmpb),
    "univariate_surrogate_clustering"=list("phat"=apply(phat1.marg,1,mean),"bhat"=apply(bhat1.marg,1,mean),
                                           "phat_all"=phat1.marg,"bhat_all"=bhat1.marg,"bptb_all"=bptb1.list),
    "univariate_surrogate_thresholding"=list("phat"=apply(phat2.marg,1,mean),"bhat"=apply(bhat2.marg,1,mean),
                                             "phat_all"=phat2.marg,"bhat_all"=bhat2.marg,"bptb_all"=bptb2.list)
  )
  )
}

#' Kernel function to perform surrogate clustering,
#' regularized regression estimation, and
#' perturbation resampling
#'
#' @param dat.S
#' @param dat.X
#' @param b0
#' @param sub.n numeric, if nonnull bootstrap a subsample of rows
#' @param mclust.modelNames see ?mclustModelNames in the mclust package
#'
#' @return list
#' @export
#'
kern_varselect <- function(dat.S, dat.X, b0, sub.n=NULL,
                           sub.sample.type="perturb",
                           mclust.modelNames=NULL) {
  ## need dat.S; dat.X both in matrix form ##
  dat.S = as.matrix(dat.S);
  tmpind = 1:nrow(dat.S)
  if(!is.null(sub.n)) {
    tmpind = sample(1:nrow(dat.S),size = sub.n)
  }
  dat.X1 = dat.X[tmpind,,drop=FALSE]
  dat.S1 = dat.S[tmpind,,drop=FALSE]

  if(length(unique(dat.S[,1])) > 2){ #if dat.S is not binary, approximate
    fit.type = "approx"
    mclustfit = mclust::Mclust(dat.S1, G=2, modelNames=mclust.modelNames)
    pi.S = ProbD.S(dat.S1,par=mclustfit$par);
  }else{
    fit.type = "exact"
    mclustfit = NULL
    pi.S = dat.S1[,1]
  }

  bhat = try(Est.ALASSO.GLM.new(cbind(pi.S,dat.X1),fit.type=fit.type))
  if(class(bhat)=="try-error") {bhat=rep(NA,ncol(dat.X))}

  tmpb = matrix(NA,nrow=b0,ncol=ncol(dat.X)+1); phat=NULL
  if(b0>0) {
  for(bb in 1:b0){
    tmpdat = cbind(pi.S,dat.X)
    tmpWi = rexp(nrow(dat.X))
    if(!is.null(sub.n)) {
      tmpind = sample(1:nrow(dat.X),size = sub.n)
      dat.X1 = dat.X[tmpind,,drop=FALSE]
      dat.S1 = dat.S[tmpind,,drop=FALSE]
      if(length(unique(dat.S1[,1])) > 2){ #if dat.S is not binary, approximate
        fit.type = "approx"
        mclustfit = mclust::Mclust(dat.S1, G=2, modelNames=mclust.modelNames)
        pi.S = ProbD.S(dat.S1,par=mclustfit$par);
      }else{
        fit.type = "exact"
        mclustfit = NULL
        pi.S = dat.S1[,1]
      }
      tmpdat = cbind(pi.S,dat.X1)
      tmpWi = rep(1,length(tmpind))
      if(sub.sample.type=="perturb") {tmpWi = rexp(length(tmpind))}
    }
    junk = Est.ALASSO.GLM.new(tmpdat, Wi=tmpWi,fit.type=fit.type);
    if(length(junk)==ncol(dat.X)+1){tmpb[bb,] = junk}
  }#end bb loop
  phat = apply(tmpb[,-1]==0,2,mean,na.rm=T) #phat does not include intercept
  }#end b0>0

  return(list("tmpb"=tmpb,"bhat"=bhat,"phat"=phat,"mclustfit.par"=mclustfit$par))
}

logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}


#' Calculate Probability Y=1 based on mclust output
#'
#' @param Si data.frame or matrix of surrogate markers
#' @param par parameter output from mclust function
#'
#' @return numeric value of length nrow(Si) estimating probability Y=1
#' @export
#'
ProbD.S = function(Si,par){
  par.list = list("pro"=par$pro, "mu"=matrix(par$mean,ncol=2),"var"=list(1,1));
  Si = as.matrix(Si); k1 = which.max(apply(par.list$mu,2,mean)); k0 = setdiff(1:2,k1)
  if(ncol(Si)==1){
    sig2 = par$variance$sigmasq; if(length(sig2)==1){sig2=rep(sig2,2)}; par.list$var = as.list(sig2)
  }else{
    for(kk in 1:2){par.list$var[[kk]] = par$variance$sigma[,,kk]}
  }
  tmp1 = dmvnorm(Si,mean=par.list$mu[,k1],sigma=as.matrix(par.list$var[[k1]]))*par$pro[k1]
  tmp0 = dmvnorm(Si,mean=par.list$mu[,k0],sigma=as.matrix(par.list$var[[k0]]))*par$pro[k0]
  tmp1/(tmp1+tmp0)
}




#' Predict D=1 via mclust on S and X
#'
#' @param dat.Xt
#' @param dat.St
#' @param dat.Xv
#' @param dat.Sv
#' @param betahat
#' @param mclustfit.par
#' @param mclust.modelNames
#'
#' @return list of predictive probabilities based on X, S, and S+X for training and validation sets
#' @export
#'
#' @examples
#' #update with an example dataset
#' #kernfit <- kern_varselect(dat.S=dat.S,dat.X=dat.X,b0=0,mclust.modelNames="EEI")
#' #predall <- ProbD.SX(dat.Xt,dat.St,dat.Xv,dat.Sv,betahat=kernfit$bhat[-1])
#'
ProbD.SX = function(dat.Xt,dat.St,dat.Xv,dat.Sv,
                    betahat,
                    sub.n=NULL,
                    mclust.modelNames="EEI") {

  dat.Xt = as.matrix(dat.Xt)
  dat.St = as.matrix(dat.St)

  pi.Xt = dat.Xt%*%betahat #betahat no intercept
  pi.Xv = dat.Xv%*%betahat

  # if sub.n is given, use a subset of data to train clustering parameters
  tmpind = 1:nrow(dat.S)
  if(!is.null(sub.n)) {
    tmpind = sample(1:nrow(dat.S),size = sub.n)
  }
  dat.X1 = dat.Xt[tmpind,,drop=FALSE]
  dat.S1 = dat.St[tmpind,,drop=FALSE]
  pi.X1 = dat.X1%*%betahat


  mclustfit = mclust::Mclust(dat.S1, G=2, modelNames=mclust.modelNames)
  pi.St = ProbD.S(cbind(dat.St),par=mclustfit$par)
  pi.Sv = ProbD.S(cbind(dat.Sv),par=mclustfit$par)

  mclustfit = mclust::Mclust(cbind(dat.S1,g.logit(pi.X1)), G=2, modelNames=mclust.modelNames)
  pi.SXt = ProbD.S(cbind(dat.St,g.logit(pi.Xt)),par=mclustfit$par)
  pi.SXv = ProbD.S(cbind(dat.Sv,g.logit(pi.Xv)),par=mclustfit$par)

  return(list("pi.SXt"=pi.SXt,"pi.SXv"=pi.SXv,"pi.Xt"=pi.Xt,"pi.Xv"=pi.Xv,"pi.St"=pi.St,"pi.Sv"=pi.Sv))
}
