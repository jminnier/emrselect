
#' Automated feature selection
#'
#' @param dat.X data.frame or matrix of covariates
#' @param dat.S data.frame or matrix of surrogate markers
#' @param b0 number of perturbation resampling draws
#'
#' @return list
#' @export
#'
emrselect <- function(dat.X,dat.S,b0=100) {
  dat.X = as.matrix(dat.X)
  A0 = ncol(as.matrix(dat.S))
  ## ========================================================================== ##
  ## multiple surroage analysis: clustering of the multivariate surrogate first ##
  ##                             then use predicted prob as pseudo outcome      ##
  ## ========================================================================== ##
  tmpout_all <- kern_varselect(dat.S=dat.S,dat.X=dat.X,b0=b0)
  ## ==================================================================================== ##
  ##                     compare to marginal single surrogate analysis                     ##
  ## **1: clustering of a univariate surrogate; then use predicted prob as pseudo outcome ##
  ## **2: threshold each of the univariate surrogate to as a pseudo outcome
  ## ==================================================================================== ##
  bptb1.list = bptb2.list = as.list(1:A0); phat1.marg = phat2.marg = bhat1.marg = bhat2.marg = NULL
  ## single Sk at a time, either first do clustring, or threshold by 1
  for(kk in 1:A0){ #A0 number of S variables
    tmpdat.S = dat.S[,kk,drop=F]; tmpout1 <- kern_varselect(dat.S=tmpdat.S,dat.X=dat.X,b0=b0)
    tmpdat.S = 1*(tmpdat.S >=1);      tmpout2 <- kern_varselect(dat.S=tmpdat.S,dat.X=dat.X,b0=b0)
    bptb1.list[[kk]] = tmpout1$tmpb; bhat1.marg = cbind(bhat1.marg, tmpout1$bhat); phat1.marg = cbind(phat1.marg,tmpout1$phat)
    bptb2.list[[kk]] = tmpout2$tmpb; bhat2.marg = cbind(bhat2.marg, tmpout2$bhat); phat2.marg = cbind(phat2.marg,tmpout2$phat)
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
#'
#' @return list
#' @export
#'
kern_varselect <- function(dat.S, dat.X, b0) {
  ## need dat.S; dat.X both in matrix form ##
  dat.S = as.matrix(dat.S);
  if(length(unique(dat.S[,1])) > 2){ #if dat.S is not binary, approximate
    fit.type = "approx"; tmpfit = mclust::Mclust(dat.S,G=2); pi.S = ProbD.S(dat.S,par=tmpfit$par);
  }else{
    fit.type = "exact"; pi.S = dat.S[,1]
  }
  tmpb = matrix(NA,nrow=b0,ncol=ncol(dat.X)+1)
  bhat = Est.ALASSO.GLM.new(cbind(pi.S,dat.X),fit.type=fit.type)
  for(bb in 1:b0){
    junk = Est.ALASSO.GLM.new(cbind(pi.S,dat.X), Wi=rexp(nrow(dat.X)),fit.type=fit.type);
    if(length(junk)==ncol(dat.X)+1){tmpb[bb,] = junk}}
  phat = apply(tmpb[,-1]==0,2,mean,na.rm=T) #phat does not include intercept
  return(list("tmpb"=tmpb,"bhat"=bhat,"phat"=phat))
}

logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}

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