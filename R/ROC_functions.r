
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

ROC.FUN.Sup.Par = function(St,Yt,fpr=jump.u){
  mhati = g.logit(St); ss = unique(sort(St)); mu1 = mean(Yt); mu0 = 1-mu1; nt = length(Yt)
  S1.ss = sum.I(ss, "<=", St, mhati)/sum(mhati)
  S0.ss = sum.I(ss, "<=", St, 1-mhati)/sum(1-mhati); auc = sum(S1.ss[-1]*(S0.ss[-nt]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr)$y}); colnames(out)=tmpnm
  list(auc,out)
}


ROC.FUN.SSL.NP = function(Sv,St,Yt,fpr=jump.u,yes.CV=F,Xt=NULL,Xv=NULL,rep=10,regularize=yes.regularize){
  nv = length(Sv); nt = length(Yt); Pt = sum.I(St, ">=", Sv)/nv
  Pv = sum.I(Sv,">=",Sv)/nv; bw = sd(Pt)/nt^0.3; ss = sort(Pv)
  mhat.v = predict(locfit(Yt~lp(Pt, deg=0, h=bw),ev=Pv))
  mu1 = mean(Yt); mu0 = 1-mu1;
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Pv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Pv, 1-mhat.v)/sum(1-mhat.v);
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Pvi = sum.I(cbind(1,Xv)%*%bhat.v,">=",Sv)/nv
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Pvi,1-mvi)/sum(1-mvi)
      }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }
  auc = sum(S1.ss[-1]*(S0.ss[-nv]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=ss,"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out)=tmpnm
  list(auc,out)
}

ROC.FUN.SSL.Par = function(Xv,dat.t,fpr=jump.u,yes.CV=F,rep=10,wgt=NULL,rtn="list",bptb=NULL,bptb2=NULL,regularize=yes.regularize){
  if(is.null(wgt)){wgt=rep(1,nrow(dat.t))}
  if(is.null(bptb)){bptb = Est.ALASSO.GLM(dat.t,rtn="EST",regularize=regularize,Wi=wgt)[1+0:ncol(Xv)]}
  Yt = dat.t[,1]; if(is.null(bptb2)){bptb2 = bptb};
  Sv = c(cbind(1,Xv)%*%bptb); mhat.v = g.logit(c(cbind(1,Xv)%*%bptb2));
  ss = unique(sort(Sv)); mu1 = mean(mhat.v); mu0 = 1-mu1; nv = length(Sv); nt = length(Yt)
  if(!yes.CV){
    S1.ss = sum.I(ss, "<=", Sv, mhat.v)/sum(mhat.v)
    S0.ss = sum.I(ss, "<=", Sv, 1-mhat.v)/sum(1-mhat.v);
  }else{
    K.fold = 2; nt.t = floor(nt/K.fold); S1.ss = S0.ss = matrix(NA,nrow=length(ss),ncol=K.fold*rep)
    for(rr in 1:rep){
      dat.t = dat.t[sample(1:nt),]
      for(kk in 1:K.fold){
        indt.v = 1:nt.t + (kk-1)*nt.t; indt.t = setdiff(1:nt, indt.v)
        bhat.t = Est.ALASSO.GLM(dat.t[indt.t,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        bhat.v = Est.ALASSO.GLM(dat.t[indt.v,],rtn="EST",regularize=regularize,BIC.power=0.1)[1+0:ncol(Xv)];
        mvi = g.logit(cbind(1,Xv)%*%bhat.t); Svi = cbind(1,Xv)%*%bhat.v
        S1.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,mvi)/sum(mvi)
        S0.ss[,kk+K.fold*(rr-1)] = sum.I(ss,"<=",Svi,1-mvi)/sum(1-mvi)
    }}
    S1.ss = apply(S1.ss,1,mean); S0.ss = apply(S0.ss,1,mean)
  }
  auc = sum(S1.ss[-1]*(S0.ss[-nv]-S0.ss[-1]))
  PPV.ss = S1.ss*mu1/(S1.ss*mu1+S0.ss*mu0); NPV.ss = (1-S0.ss)*mu0/((1-S0.ss)*mu0+(1-S1.ss)*mu1)
  out=cbind("cut"=g.logit(ss),"FPR"=S0.ss,"TPR"=S1.ss,"PPV"=PPV.ss,"NPV"=NPV.ss); tmpnm=colnames(out)
  out=sapply(1:ncol(out),function(kk){approx(out[,"FPR"],out[,kk],fpr,rule=2)$y}); colnames(out)=tmpnm
  if(rtn=="vec"){return(c(auc,out))}else{return(list(auc,out))}
}

logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}




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

ROC.FUN.ALASSO.cv  <- function(ind=NULL,data, wgti=NULL, Wi=NULL, yy0 = 0.5, nopen.ind=NULL,
                            FPR0=seq(.01,.99,by=.01),rtn="ALL",yes.CV=F,yes.seed=F,rep=10,regularize=T,lam.ridge=0)
  {
    ## ========================================================================= ##
    ## ==== ind: for bootstrap methods; data[ind,] returns bootstrap samples === ##
    ## ==== 1st column is binary disease status;                            ==== ##
    ## ==== 2nd column and on are individual markes                         ==== ##
    ## ==== yy0: prespecified cut-off;  FPR0: prespecified Specificity level ==== ##
    ## ========================================================================= ##

    nn <- nrow(data); if(is.null(wgti)){wgti=rep(1,nn)};
    if(!is.null(ind)){data <- data[ind,]; wgti = wgti[ind]}
    n.set <- 1; pp = ncol(data)

    ## ========================= ##
    ## === Apparent Accuracy === ##
    ## ========================= ##

  	betahat = try(Est.ALASSO.GLM(data,Wi=wgti, rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
    if(length(betahat)>1)
      {
        yyi = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:pp]);
        q.nonzero = sum(abs(betahat[2:(pp)])!=0)
        if(q.nonzero==0){rochat = NULL}else{rochat = ROC.Est.FUN(data[,1],yyi,yy0,FPR0,wgti=wgti)}
        if(yes.seed){set.seed(315)}
        if(yes.CV & q.nonzero >0)
        {
            ## =============================== ##
            ## === K-fold cross validation === ##
            ## =============================== ##
            K = 2; nv = floor(nn/K); bias.cv = NULL
            for(i in 1:rep)
              {
                tmpind=sample(1:nn); data = data[tmpind,]; wgti = wgti[tmpind]
                for(k in 1:K)
                 {
                    ind.v = 1:nv + (k-1)*nv; ind.t = setdiff(1:nn,ind.v)
                    dat.t = data[ind.t,]; dat.v = data[ind.v,]
                    ## ============================================================================== ##
                    ## ==== Calculating (1) Coef for Combining Markers (beta) with Training Data ==== ##
                    ## ====             (2) Combined Marker Value (yyi.new) with Validation Data ==== ##
                    ## ============================================================================== ##
                    beta.t = try(Est.ALASSO.GLM(dat.t,Wi=wgti[ind.t],rtn="EST",regularize=regularize,nopen.ind=nopen.ind,lam.ridge=lam.ridge),silent=T)
                    if(length(beta.t)>1)
                      {
                        beta.t = beta.t[1:pp]
                        yyi.v = g.logit(cbind(1,as.matrix(dat.v[,-1]))%*%beta.t)
                        yyi.t = g.logit(cbind(1,as.matrix(dat.t[,-1]))%*%beta.t)
                        bias.k = try(ROC.Est.FUN(dat.t[,1],yyi.t,yy0,FPR0,wgti[ind.t])-ROC.Est.FUN(dat.v[,1],yyi.v,yy0,FPR0,wgti[ind.v]),silent=T)
                        if(length(bias.k)>1){bias.cv = cbind(bias.cv, bias.k)}
                      }
                }
              }
            print(ncol(bias.cv)); bias.cv = apply(bias.cv,1,mean,trim=0.05,na.rm=T)/2; roc.cv = rochat - bias.cv
        }else{
            roc.cv = NULL
            }
     }else{
      roc.cv=beta.hat=NULL
     }
   if(rtn=="EST"){
     if(yes.CV){out=c(roc.cv,betahat)}else{out=c(rochat,betahat)}
     return(out)
   }else{return(list(rochat, roc.cv, betahat, yyi))}
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

#' AUC function
#'
#' @param data
#'
#' @return AUC
#' @export
#'
#' @examples AUC.FUN(cbind(rbinom(10,size=2,prob=.5),rnorm(10)))
AUC.FUN = function(data)
  {
	dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd)
	x0 = xx[dd==0]; x1 = xx[dd==1]
  	sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
  }

chat.FUN.632boot <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data); data = as.matrix(data); nn = length(yy)
    set.seed(1202); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:20)
      {
      	ind.t = sample(1:nn,rep=T); ind.v = setdiff(1:nn,ind.t)
        datat = data[ind.t,]; datav = data[ind.v,]
        beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
        if(length(beta)>1){
            phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
            c.fpr   = as.numeric(try(quantile(phatv[datav[,1]==0],1-fpr0),silent=T))
            c.all = cbind(c.all,c.fpr); ##c.t.all = c(c.t.all, c.fpr.t);
         }
      }
    phat0 = c(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat0[data[,1]==0],1-fpr0)
    c.fpr.bc  = 0.368*c.fpr.0 + 0.632*apply(as.matrix(c.all),1,mean,na.rm=T)
    c(g.logit(c.fpr.0),g.logit(c.fpr.bc))
  }

chat.FUN <- function(data,fpr0=0.05,betahat,regularize=T,nopen.ind=NULL)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    phat.0 = (cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    set.seed(1202); nt = floor(nrow(data)/2); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      {
        data = data[sample(1:nrow(data)),]
        for(k in 1:2)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
        	beta = try(Est.ALASSO.GLM(datat,nopen.ind=nopen.ind,regularize=regularize,rtn="EST")[1:p0],silent=T)
            if(length(beta)>1){
                phatv = (cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = (cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0)
                c.all = rbind(c.all,c.fpr); c.t.all = rbind(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    c.fpr.bc  = c.fpr.0 + apply(c.all-c.t.all,2,mean,trim=0.05)
    g.logit(c(c.fpr.0,c.fpr.bc))
  }

Predict.FUN <- function(data,newdata,fpr0=0.05,betahat)
  {
    yy = data[,1]; xx = as.matrix(data[,-1]); p0 = ncol(data)
    set.seed(1202); nt = floor(nrow(data)/3); tmpb.all = c.t.all = c.all = NULL
    for(rep in 1:10)
      {
        data = data[sample(1:nrow(data)),]
        for(k in 1:3)
          {
            tmpind = 1:nt+(k-1)*nt
            datat = data[-tmpind,]; datav = data[tmpind,]
            beta = try(Est.ALASSO.GLM(datat,nopen.ind=c(1,2))[[1]][1:p0])
            if(length(beta)>1){
                phatv = g.logit(cbind(1,as.matrix(datav[,-1]))%*%beta)
                phatt = g.logit(cbind(1,as.matrix(datat[,-1]))%*%beta)
                c.fpr   = quantile(phatv[datav[,1]==0],1-fpr0)
                c.fpr.t = quantile(phatt[datat[,1]==0],1-fpr0)
                c.all = c(c.all,c.fpr); c.t.all = c(c.t.all, c.fpr.t)
                tmpb.all = cbind(tmpb.all,beta)}
          }
      }
    #pnew.all = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%tmpb.all)
    pnew.0 = g.logit(cbind(1,as.matrix(newdata[,-1]))%*%betahat[1:p0])
    phat.0 = g.logit(cbind(1,as.matrix(data[,-1]))%*%betahat[1:p0])
    c.fpr.0   = quantile(phat.0[data[,1]==0],1-fpr0)
    c.fpr.bc  = c.fpr.0 + mean(c.all-c.t.all)
    data.frame("patient_num"=newdata$patient,"Prob.RA"=round(pnew.0,5),
    "Cut-off.0"=round(c.fpr.0,5), "Cut-off.bc"=round(c.fpr.bc,5))
  }


