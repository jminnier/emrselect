
logit = function(xx){log(xx/(1-xx))};  g = function(xx){exp(xx)/(1+exp(xx))};
d1g = function(xx){g(xx)*(1-g(xx))}; d2g = function(xx){d1g(xx)*(1-2*g(xx))}

g.logit = function(xx){exp(xx)/(exp(xx)+1)}
logit = function(xx){log(xx/(1-xx))}
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}
Intercept.GLM = function(yi,bxi){glm(yi~bxi,family="binomial")}


ProbTestPos <- function(tpr,fpr,prev)
{
  tpr*prev + fpr*(1-prev)
}


L2Norm <- function(data,coef,intercept=F)
{
  yy = data[,1]; xx = data[,-1]; if(intercept){xx=cbind(1,xx)}
  apply((yy - xx%*%t(coef))^2,2,mean)
}


PPV.FUN <- function(fpr,SE, mu0){ 1/(1+fpr/SE*(1-mu0)/mu0)}
NPV.Project <- function(npv.e0y0,p.e1,p.y0.e0,npv.e1=1)
{
  ## P(D=0 | E=1 or E=0&Y=0) = {P(D=0|E=1)P(E=1)+P(D=0|E=0&Y=0)P(E=0&Y=0)}/{P(E=1)+P(E=0&Y=0)
  npv.all = (npv.e1 * p.e1 + npv.e0y0*p.y0.e0*(1-p.e1))/(p.e1+p.y0.e0*(1-p.e1))
  npv.all
}

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

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth)
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

