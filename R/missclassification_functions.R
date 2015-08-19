
ROC.FUN.ALASSO.632boot  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                                    FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T){
  ## ========================================================================= ##
  ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
  ## ==== 1st column is binary disease status;                            ==== ##
  ## ==== 2nd column and on are individual markes                         ==== ##
  ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
  ## ========================================================================= ##
  nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)};
  if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}
  n.set <- 1; pp = ncol(data);  yyi.vmat = matrix(NA, nrow=nn, ncol=rep)
  ## === Apparent Accuracy === ##
  betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
  if(length(betahat)>1){
    yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp])
    q.nonzero = sum(abs(betahat[2:(pp)])!=0)
    if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
    if(yes.seed){set.seed(315)}
    if(yes.CV&q.nonzero>0){
      ## === 0.632 bootstrap === ##
      roc.cv = NULL;
      for(i in 1:rep){
        tmpind=sample(1:nn,replace=T); ind.v = setdiff(1:nn,unique(tmpind))
        dat.t = data[tmpind,]; wgti.t = wgti[tmpind]; dat.v = data[ind.v,];  wgti.v = wgti[ind.v]
        beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti.t,rtn="EST",regularize=regularize,nopen.ind=nopen.ind),silent=T)
        if(length(beta.t)>1){
          beta.t = beta.t[1:pp]; yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t); yyi.vmat[ind.v,i]=yyi.v
          roc.k = try(ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti.v),silent=T)
          if(length(roc.k)>1){roc.cv = cbind(roc.cv, roc.k)}  } }
      roc.cv = apply(roc.cv,1,mean)*0.632+rochat*0.368
    }else{roc.cv = NULL}
  }else{roc.cv=beta.hat=NULL}
  if(rtn=="EST"){
    if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
    return(out)
  }else{return(list("rochat"=rochat, "roc.cv"=roc.cv, "beta"=betahat, "Si"=yyi, "yyi.v"=yyi.vmat))}
}


ROC.Est.FUN <- function(Di,yyi,yy0=0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0
  for(k in 1:pp)
  {
    yy = yy0;
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y;
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR);
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
    #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}


AUC.FUN = function(data)
{
  dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd)
  x0 = xx[dd==0]; x1 = xx[dd==1]
  sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
}

#' PPV.FUN
#' @export
PPV.FUN <- function(fpr,SE, mu0){ 1/(1+fpr/SE*(1-mu0)/mu0)}

#' NPV.Project
#' @export
NPV.Project <- function(npv.e0y0,p.e1,p.y0.e0,npv.e1=1)
{
  ## P(D=0 | E=1 or E=0&Y=0) = {P(D=0|E=1)P(E=1)+P(D=0|E=0&Y=0)P(E=0&Y=0)}/{P(E=1)+P(E=0&Y=0)
  npv.all = (npv.e1 * p.e1 + npv.e0y0*p.y0.e0*(1-p.e1))/(p.e1+p.y0.e0*(1-p.e1))
  npv.all
}