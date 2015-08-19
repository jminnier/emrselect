
#' Simulate multiG
#'
#' @param nn
#' @param b.G
#' @param b.S
#' @param rho.S
#'
#' @return dat
#' @export
#'
SIM.FUN.MultiG = function(nn,b.G,b.S,rho.S){## p: # of SNPs; K: # of algorithms
  p = length(b.G)-1; G = matrix(c(rnorm(nn*2), rbinom(nn*(p-2),prob=0.3,size=2)),nrow=nn)
  Y = rbinom(nn,size=1,prob=g(b.G[1]+c(G%*%b.G[-1])))
  A = length(b.S); S = mvrnorm(nn,mu=c(-1.5,-1,rep(-0.5,A-2)),Sigma=rho.S+(1-rho.S)*diag(A)) + VTM(b.S,nn)*Y;
  dat=cbind(Y,G,S); colnames(dat)=c("Y",paste0("G",1:p),paste0("S",1:A))
  dat
}

#' Simulate data
#'
#' @param nn
#' @param b.G
#' @param b.S
#' @param rho.S
#' @param rho.G
#' @param ng.cts
#' @param Gtype
#' @param Stype
#'
#' @return dat
#' @export
#'
#' @examples
#'
#' b0.G = rep(c(0.3,0.5,0.3,0.2,0),c(1,2,3,2,50)); b0.G = rbind(-1,b0.G); p = ncol(b0.G)
#' b0.S = c(1.5,1.5,2)
#' mydata0 = SIM.FUN(NN,b.G=b0.G,b.S=b0.S,rho.S=0.3, rho.G=0.2);
SIM.FUN = function(nn,b.G,b.S,rho.S,rho.G,ng.cts,
                   Gtype="binomial", # or lognormal
                   Stype="cts" #or count
){

  ## Y is binary
  Y = rbinom(nn,size=1,prob=0.3)

  ## Mean of G
  tmp.mean = VTM(b.G[1,],nn) + Y*VTM(b.G[2,],nn)
  tmp.prob = g(tmp.mean)

  ## Surrogates S
  A = length(b.S); p = ncol(b.G)
  S = mvrnorm(nn,mu=c(-1.5,-1,rep(-0.5,A-2)),Sigma=rho.S+(1-rho.S)*diag(A)) + VTM(b.S,nn)*Y

  #ggplot(melt(data.frame(S)),aes(x=value,color=variable))+geom_density()
  #S2 = round(exp(S-1))
  #ggplot(melt(data.frame(S2)),aes(x=value,color=variable))+geom_density()
  if(Stype=="count") {S = round(exp(S-1))}#else leave S as ctns normal

  ## Simulate S
  if(Gtype=="binomial") {
    G = cbind(matrix(rnorm(nn*ng.cts,mean=tmp.mean[,(1:ng.cts)]),nrow=nn), #cts G
              matrix(rbinom(nn*(p-ng.cts),size=2,prob=c(tmp.prob[,-(1:ng.cts)])),nrow=nn)) #count G
  }else if(Gtype=="lognormal") {
    sig.G = autocorr.mat(p=p,rho=rho.G)
    GN = t(apply(tmp.mean, 1, function(k) rmvnorm(1,mean=k,sigma=sig.G)))
    G = cbind(GN[,1:ng.cts],round(exp(GN[,-(1:ng.cts)])))
  }else if(Gtype=="normal") {
    GN = t(apply(tmp.mean, 1, function(k) rmvnorm(1,mean=k,sigma=sig.G)))
  }
  dat=cbind(Y,G,S); colnames(dat)=c("Y",paste0("G",1:p),paste0("S",1:A))
  dat
}

SIM.FUN.old2 = function(nn,b.G,b.S,rho.S){## p: # of SNPs; K: # of algorithms
  Y = rbinom(nn,size=1,prob=0.3)
  tmp.prob = g(VTM(b.G[1,],nn) + Y*VTM(b.G[2,],nn))
  A = length(b.S)
  S = MASS::mvrnorm(nn,mu=c(-1.5,-1,rep(-0.5,A-2)),
                    Sigma=rho.S+(1-rho.S)*diag(A)) + VTM(b.S,nn)*Y;
  p = ncol(b.G)
  G = cbind(rnorm(nn,mean=logit(tmp.prob)),matrix(rbinom(nn*(p-1),size=2,prob=c(tmp.prob)),nrow=nn))
  dat=cbind(Y,G,S); colnames(dat)=c("Y",paste0("G",1:p),paste0("S",1:A))
  dat
}


SIM.FUN.old = function(nn,rtn="data.t")
{
  xx = mvrnorm(nn,mu=rep(0,p.x),Sigma=Sig0.X);
  icd.B = rbinom(nn,size=1,prob=g.logit(xx[,1]*3+2))
  xx[,1] = (xx[,1]*(xx[,1]>0)+(rexp(nn,rate=0.1)+5)*rbinom(nn,size=1,prob=0.1))*icd.B
  prob.x = g.logit(-alp0+c(xx%*%beta0)-3*(1-icd.B)+0.2*(xx[,1]>15))
  yy = rbinom(nn,prob=prob.x,size=1); dat = cbind(yy,xx)
  if(rtn=="data.t"){return(dat)}else{
    zz = rbinom(nn, size=2, prob=g.logit(log(maf)+gam.z*yy))
    return(list(dat,cbind("D"=yy,"P.x"=prob.x,"G"=zz)))}
}
