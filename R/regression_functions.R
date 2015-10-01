

Est.ALASSO.GLM.new = function(data,Wi=NULL,fit.type,nopen.ind=NULL,BIC.factor=0.1,offset=NULL,regularize=T){
  if(fit.type=="exact"){
    Est.ALASSO.GLM(data,Wi=Wi,nopen.ind=nopen.ind,BIC.factor=BIC.factor,offset=offset,fam0="binomial",regularize=regularize)
  }else{
    Est.ALASSO.GLM.Approx(data,Wi=Wi,nopen.ind=nopen.ind,BIC.factor=BIC.factor,offset=offset)
  }
}

#' GLM with approximation
#'
#' @param data 1st column y; remaining x
#' @param Wi weights for resampling if interested in obtaining standard errors
#' @param rtn
#' @param nopen.ind indexes which subset of x should not be penalized
#' @param BIC.factor
#' @param offset
#'
#' @return bhat
#' @export
#'
Est.ALASSO.GLM.Approx = function(data,
                                 Wi=NULL,
                                 rtn="EST",
                                 nopen.ind=NULL,
                                 BIC.factor=0.1,
                                 offset=NULL){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x)
  if(is.null(Wi)){Wi=rep(1,nn)}; if(is.null(offset)){offset=rep(0,nn)}

  #initial beta
  btilde = glm(y~x,family=binomial,weight=Wi)

  #Jacobian
  Ahat = solve(summary(btilde)$cov.unscaled); btilde = btilde$coef
  Ahat.half = svd(Ahat); Ahat.half = Ahat.half$u%*%diag(sqrt(Ahat.half$d))%*%t(Ahat.half$u)
  Xtilde = Ahat.half; Ytilde = Ahat.half%*%btilde

  #trick to force lasso to estimate adaptive lasso
  w.b = 1/abs(btilde); Xtilde.t = Xtilde/VTM(w.b,nrow(Xtilde))

  #lasso
  tmpfit = lars::lars(Xtilde.t,Ytilde,type="lasso",normalize=F,intercept=F)
  lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
  b.all = lars::predict.lars(tmpfit, s=lam.all, type="coefficients",mode="lambda")$coef
  b.all = b.all/VTM(w.b,nrow(b.all)); m0 = length(lam.all)
  df.all = apply(b.all[,-1,drop=F]!=0,1,sum)+1;

  ## =============================================================================================================== ##
  ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
  ## =============================================================================================================== ##
  BIC.lam = -2*logitlik.fun(t(b.all),dat=data)+min(nn^BIC.factor,log(nn))*df.all
  m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  bhat
}



#' Logistic regression or logistic adaptive LASSO estimator.
#' Initial estimator obtained via ridge
#'
#' @param data 1st column y; remaining x
#' @param Wi weights for resampling if interested in obtaining standard errors
#' @param rtn
#' @param nopen.ind
#' @param regularize binary; if FALSE estimates GLM, if TRUE estimates adaptive LASSO GLM
#' @param BIC.factor
#' @param offset
#' @param fam0
#'
#' @return beta
#' @export
#'
Est.ALASSO.GLM = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize,
                          BIC.factor=0.1,offset=NULL,fam0="binomial"){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x)
  if(is.null(Wi)){Wi=rep(1,nn)}; if(is.null(offset)){offset=rep(0,nn)}
  if(regularize){
    ## ================================================================================ ##
    ## ridge regression initial estimator
    ## ================================================================================ ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef;
    bini = as.vector(coef(glmnet::glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family=fam0,offset=offset)))

    ## ================================================================================ ##
    ## adaptive weights for aLASSO
    ## ================================================================================ ##
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))

    ## ================================================================================ ##
    ## glmpath provides solution path for a range of penalty parameters ##
    ## ================================================================================ ##
    tmpfit = glmpath::glmpath(x.t,y,nopenalty.subset=nopen.ind,family=fam0,weight=Wi,standardize=F,min.lambda=0,
                     lambda2=lam.ridge,offset=offset)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = glmpath::predict.glmpath(tmpfit, s=lam.all, type="coefficients",mode="lambda",offset=offset)
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all betas (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(glmpath::predict.glmpath(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda",offset=offset),2,sum)+min(sum(y)^BIC.factor,log(sum(y)))*df.all
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=fam0,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  #out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  out = c("b"=bhat,"bini"=bini,"lamhat"=lamhat)
  if(rtn=="EST"){return(bhat)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}


#' Cross-validated LASSO
#'
#' @param data
#'
#' @return beta
#' @export
#'
CV.FUN <- function(data)
{
  y = data[,1]; x = data[,-1]; nn=length(y)
  w = 1/abs(lm(y~x)$coef[-1]); x.w = x/VTM(w,nrow(x))
  s0.all = (1:100)/100
  fit.all = l1ce(y~x.w, standardize=F, bound=s0.all)
  K = 5; L2cv = NULL
  for(k in 1:K)
  {
    indv = 1:floor(nn/K) + (k-1)*floor(nn/K)
    indt = setdiff(1:nn,indv)
    fitk = l1ce(y~x.w,subset=indt, standardize=F, bound=s0.all)
    L2cv = rbind(L2cv, L2Norm(cbind(y,x.w)[indv,],coef(fitk),intercept=T))
  }
  L2cv = apply(L2cv,2,mean)
  s0 = min(s0.all[L2cv==min(L2cv)])
  bhat  = l1ce(y~x.w, standardize=F, bound=s0)
  list("b"=coef(bhat)/w, "s0"=s0, "lam0" =bhat$Lagrangian,"b0"=bhat$bound)
}


