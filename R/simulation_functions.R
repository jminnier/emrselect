
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
#'
#' @return dat
#' @export
#'
#' @examples
#'
#' b0.G = rep(c(0.3,0.5,0.3,0.2,0),c(1,2,3,2,50)); b0.G = rbind(-1,b0.G); p = ncol(b0.G)
#' b0.S = c(1.5,1.5,2); rho0.S = 0.3; A0 = length(b0.S);
#' mydata0 = SIM.FUN(NN,b.G=b0.G,b.S=b0.S,rho.S=rho0.S);
SIM.FUN = function(nn,b.G,b.S,rho.S){## p: # of SNPs; K: # of algorithms
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
