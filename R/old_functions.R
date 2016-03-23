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

Par.FUN = function(dat){
  Y = dat[,"Y"]
  S = dat[,substring(colnames(dat),1,1)=="S",drop=F]
  G = dat[,substring(colnames(dat),1,1)=="G"]

  # surrogate
  uni.S = unique(sort(c(S))); K = length(uni.S); nn = nrow(dat); A = ncol(S)
  Dmat = Imat = matrix(0,nrow=nn*A,ncol=K*A);
  for(aa in 1:A){
    Imat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=1*(S[,aa]==VTM(uni.S,nn));Dmat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=1
  }
  bet = glm(Y~G,family=binomial)$coef
  alp1 = matrix(t(rep(Y,A))%*%Imat/t(rep(Y,A))%*%Dmat,nrow=K);
  alp0 = matrix(t(rep(1-Y,A))%*%Imat/t(rep(1-Y,A))%*%Dmat,nrow=K);
  row.names(alp1)=paste0("SE",1:K); row.names(alp0)=paste0("SP",1:K)
  colnames(alp1)=colnames(alp0)=paste0("S",1:A)
  list("bet"=bet,"alp1"=alp1,"alp0"=alp0)
}

Est.FUN2 = function(dat,bet.ini,alp1.ini,alp0.ini,alp.known=F){ ##alp1,alp0 is K x A, remove first category
  p = length(bet.ini)-1; Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  uni.S = unique(sort(c(Si))); K=length(uni.S) ## remove first category
  Z = cbind(1,Gi[rep(1:nn,A),]) ## (nn*A) x (p+1) matrix
  eps = 1; Dmat = Imat = matrix(0,nrow=nn*A,ncol=K*A);
  for(aa in 1:A){Imat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=1*(Si[,aa]==VTM(uni.S,nn));Dmat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=1}
  bet = c(bet.ini); alp1 = c(alp1.ini); alp0 = c(alp0.ini); theta.old = c(bet,alp1,alp0)
  while(eps > 1e-5){
    eps2 = 1
    while(eps2 > 1e-4){
      betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
      P = pi*c(Imat%*%alp1)+(1-pi)*c(Imat%*%alp0)
      d1P = c(Imat%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
      score = c("S.bet"=t(d1P/P)%*%Z); Info = -t(Z)%*%(Z*d1P^2/P^2) ##Info = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2);
      bet.new = bet - c(solve(Info)%*%score); eps2 = sqrt(sum((bet.new-bet)^2)); bet = bet.new}
    Yihat.a = g(betG+log((Imat%*%alp1)/(Imat%*%alp0))) ## n x A matrix
    if(!alp.known){
      alp1.new = c(matrix(c(t(Yihat.a)%*%(Imat))/c(t(Yihat.a)%*%Dmat),ncol=A));
      alp0.new = c(matrix(c(t(1-Yihat.a)%*%Imat)/c(t(1-Yihat.a)%*%Dmat),ncol=A))
    }else{alp1.new = alp1; alp0.new = alp0}
    score = c("S.alp1"=t(pi/P)%*%Imat,"S.alp0"=t((1-pi)/P)%*%Imat)
    eps = sqrt(sum(c(alp1.new-alp1,alp0.new-alp0)^2));
    alp1 = alp1.new; alp0 = alp0.new; print(c(eps,score)); print(cbind(alp1,c(alp1.ini),alp0,c(alp0.ini)))
    #print(c(step,eps)); print(bet); print(cbind(matrix(alp1,nrow=A,byrow=T),matrix(alp0,nrow=A,byrow=T)))
    #if(step>2000){if(eps.vec[pmin(step,5000)]>eps.vec[pmin(step,5000)-200]){break}}
  }
  alp1 = matrix(alp1,ncol=A); alp0 = matrix(alp0,ncol=A);colnames(alp0)=colnames(alp1) = paste0("S",1:A);
  row.names(alp1) = paste0("SE",1:K); row.names(alp0)=paste0("SP",1:K)
  if(eps > 1e-5){bet=bet*NA;alp=alp*NA}
  list("bet"=bet,"alp1"=alp1,"alp0")
}

Est.FUN = function(dat,bet.ini,alp1.ini,alp0.ini,alp.known=F){ ##alp1,alp0 is K-1 x A, remove first category
  p = length(bet.ini)-1; Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  uni.S = unique(sort(c(Si))); s1 = uni.S[1]; uni.S = uni.S[-1]; K1=length(uni.S); K=K1+1## remove first category
  Z = cbind(1,Gi[rep(1:nn,A),]) ## (nn*A) x (p+1) matrix
  Dmat = Imat = matrix(0,nrow=nn*A,ncol=K*A);Itil = matrix(0,nrow=nn*A,ncol=K1*A); I1=rep(0,nn*A)
  for(aa in 1:A){tmpI = 1*(Si[,aa]==VTM(c(s1,uni.S),nn));Imat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=tmpI;
  tmpI = tmpI[,-1,drop=F]-1*(Si[,aa]==s1);Itil[1:nn+(aa-1)*nn,1:K1+(aa-1)*K1]=tmpI
  I1[1:nn+(aa-1)*nn]=1*(Si[,aa]==s1);Dmat[1:nn+(aa-1)*nn,1:K+(aa-1)*K]=1}
  step=eps=1; bet = c(bet.ini); alp1 = c(alp1.ini); alp0 = c(alp0.ini); eps.vec=rep(0,5000)
  while(eps > 1e-5){
    eps2=step1=1
    while(eps2 > 1e-5){
      betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
      P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
      d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
      score = c("S.bet"=t(d1P/P)%*%Z); Info = -t(Z)%*%(Z*d1P^2/P^2) ##Info = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2);
      bet.new = bet - c(solve(Info)%*%score); eps2 = sqrt(sum((bet.new-bet)^2)); bet = bet.new; step1=step1+1
      if(step1>100){break}}
    Yihat.a = g(betG+log((I1+Itil%*%alp1)/(I1+Itil%*%alp0))) ## n x A matrix
    if(!alp.known){
      alp1.new = c(matrix(c(t(Yihat.a)%*%(Imat))/c(t(Yihat.a)%*%Dmat),ncol=A)[-1,]);
      alp0.new = c(matrix(c(t(1-Yihat.a)%*%Imat)/c(t(1-Yihat.a)%*%Dmat),ncol=A)[-1,])
    }else{alp1.new = alp1; alp0.new = alp0}
    score = c("S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
    eps = sqrt(sum(c(alp1.new-alp1,alp0.new-alp0)^2)); eps.vec[step]=eps;
    alp1 = alp1.new; alp0 = alp0.new; ##print(c(eps,step,step1,score)); print(cbind(alp1,c(alp1.ini),alp0,c(alp0.ini)))
    #print(c(step,eps)); print(bet); print(cbind(matrix(alp1,nrow=A,byrow=T),matrix(alp0,nrow=A,byrow=T)))
    if(step>2000){if(eps.vec[pmin(step,5000)]>(min(eps.vec[1:step])*10)&max(abs(score))>1){break}}; step=step+1;
  }
  alp1 = matrix(alp1,ncol=A); alp0 = matrix(alp0,ncol=A);colnames(alp0)=colnames(alp1) = paste0("S",1:A);
  row.names(alp1) = paste0("SE",1+1:K1); row.names(alp0)=paste0("SP",1+1:K1)
  if(eps > 1e-5){bet=bet*NA;alp1=alp1*NA;alp0=alp0*NA}
  list("bet"=bet,"alp1"=alp1,"alp0"=alp0)
}

Est.FUN.uniS = function(dat,bet.ini,alp.ini){
  p = length(bet.ini)-1; Gi = cbind(1,dat[,1:p]); Si = dat[,-(1:p)];  eps = 1; k.step = 0; nn = nrow(dat);
  uni.S = sort(unique(c(Si))); K = length(uni.S); alp1.ini = alp.ini[1:K]; alp0.ini = alp.ini[-(1:K)]
  Ii = 1*(Si== VTM(uni.S, nn)) ## (nn x K)
  bet = bet.ini; alp1 = alp1.ini; alp0 = alp0.ini   ## bet: intercept, bet (2 x p matrix); alp: A x K matrix
  while(eps > 5e-5){
    alp1.ini = alp1; alp0.ini = alp0; bet.ini = bet; k.step = k.step+1
    pi.i = g(c(Gi%*%bet)); pi.bar = mean(pi.i);
    Yihat = g(logit(pi.i) + log(Ii%*%alp1) - log(Ii%*%alp0))
    Yihat[Ii%*%alp1 == 0] = 0; Yihat[Ii%*%alp0 == 0] = 1;
    bet = glm(Yihat ~ Gi[,-1],family=binomial)$coef
    alp1 = c(t(Yihat)%*%Ii)/sum(Yihat); alp0 = c(t(1-Yihat)%*%Ii)/sum(1-Yihat)
    eps = sum(abs(c(bet-bet.ini,alp1-alp1.ini,alp0-alp0.ini)));
    if(is.na(eps)){browser()}
    if(floor(k.step/500)==ceiling(k.step/500)){print(c(k.step,eps,round(bet,3))); print(round(c(alp1,alp0),3))}
  }; print(k.step)
  alp = c(alp1,alp0); names(alp) = c(paste0("SE",1:K),paste0("SP",1:K))
  list("bet"=bet,"alp"=alp)
}



Par.Marg.FUN = function(dat){
  Y = dat[,"Y"]; S = dat[,substring(colnames(dat),1,1)=="S",drop=F]; G = dat[,substring(colnames(dat),1,1)=="G"];
  A = ncol(S); uni.S = unique(sort(c(S))); K = length(uni.S); tmpn = nrow(dat);
  bet = apply(G,2,function(g,y){glm(y~g,family=binomial)$coef},y=Y)
  alp1 = matrix(t(Y)%*%(S[,rep(1:A,K)]==VTM(rep(uni.S,rep(A,K)),tmpn)),nrow=A)/sum(Y)     ## AxK
  alp0 = matrix(t(1-Y)%*%(S[,rep(1:A,K)]==VTM(rep(uni.S,rep(A,K)),tmpn)),nrow=A)/sum(1-Y) ## AxK
  row.names(alp1) = row.names(alp0) = paste0("S",1:A);
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1:K),paste0("SP",1:K))
  list("bet"=bet,"alp"=alp)
}

Est.Marg = function(dat,bet.ini,alp1.ini,alp0.ini,alp.known=F){ ##alp1,alp0 is A x K-1, remove first category
  p = ncol(bet.ini); Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  K1 = ncol(alp1.ini); uni.S = unique(sort(c(Si))); s1 = uni.S[1]; uni.S = uni.S[-1] ## remove first category
  Z = matrix(0,nrow=p*nn,ncol=p*2); for(jj in 1:p){Z[(1:nn)+nn*(jj-1),1:2+2*(jj-1)]=cbind(1,Gi[,jj])}
  Z = Z[rep(1:nrow(Z),A),] ## (nn*p)*A x 2p matrix
  Dmat = Imat = matrix(0,nrow=nn*p*A,ncol=K*A);
  Itil = matrix(0,nrow=nn*p*A,ncol=K1*A); eps = 1; step=0; I1 = rep(0,nn*p*A)
  for(aa in 1:A){tmpI = 1*(Si[,aa]==VTM(c(s1,uni.S),nn)); Imat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=tmpI[rep(1:nn,p),]
  tmpI = tmpI[,-1,drop=F]- 1*(Si[,aa]==s1);      Itil[1:(nn*p)+(aa-1)*(nn*p),1:K1+(aa-1)*K1]=tmpI[rep(1:nn,p),]
  I1[1:(nn*p)+(aa-1)*(nn*p)]=1*(Si[rep(1:nn,p),aa]==s1);Dmat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=1}
  bet = c(bet.ini); alp1 = c(t(alp1.ini)); alp0 = c(t(alp0.ini));
  while(eps > 1e-5){
    eps2 = 1; step1= 1
    while(eps2 > 1e-5){
      betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
      P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
      d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
      score = c("S.bet"=t(d1P/P)%*%Z); Info = -t(Z)%*%(Z*d1P^2/P^2) ##Info = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2);
      bet.new = bet - c(solve(Info)%*%score); eps2 = sqrt(sum((bet.new-bet)^2)); bet = bet.new; step1=step1+1}
    Yihat.ja = g(betG+log((I1+Itil%*%alp1)/(I1+Itil%*%alp0))) ## n x A matrix
    if(!alp.known){
      alp1.new = c(matrix(c(t(Yihat.ja)%*%(Imat))/c(t(Yihat.ja)%*%Dmat),ncol=A)[-1,]);
      alp0.new = c(matrix(c(t(1-Yihat.ja)%*%Imat)/c(t(1-Yihat.ja)%*%Dmat),ncol=A)[-1,])
    }else{alp1.new = alp1; alp0.new = alp0}
    score = c("S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
    eps = sqrt(sum(c(alp1.new-alp1,alp0.new-alp0)^2));
    alp1 = alp1.new; alp0 = alp0.new; ##print(c(eps,score)); print(c(bet,alp1,alp0))
    #print(c(step,eps)); print(bet); print(cbind(matrix(alp1,nrow=A,byrow=T),matrix(alp0,nrow=A,byrow=T)))
    #if(step>2000){if(eps.vec[pmin(step,5000)]>eps.vec[pmin(step,5000)-200]){break}}
  }
  alp1 = t(matrix(alp1,ncol=A)); alp0 = t(matrix(alp0,ncol=A))
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1+1:K1),paste0("SP",1+1:K1))
  bet = matrix(bet,ncol=p); if(eps > 1e-5){bet=bet*NA;alp=alp*NA}
  list("bet"=bet,"alp"=alp)
}

Loglik.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  neglik = -sum(log(P)); print(neglik)
  neglik
}

Score.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
  -c("S.bet"=t(d1P/P)%*%Z,"S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
}

Info.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
  Info = matrix(0,ncol=p*2+A*K1*2,nrow=p*2+A*K1*2)
  tmpout = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2); indr=indc=1:(2*p); Info[indr,indc]=tmpout
  tmpout = -t(Z)%*%(Itil*(pi*d1P-P*d1pi)/P^2);indr=2*p+1:(A*K1); Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
  tmpout = -t(Z)%*%(Itil*((1-pi)*d1P+P*d1pi)/P^2);indr=indr+A*K1;Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
  tmpout = -t(Itil)%*%(Itil*pi^2/P^2); indr=indc=2*p+1:(A*K1); Info[indr,indc]=tmpout
  tmpout = -t(Itil)%*%(Itil*pi*(1-pi)/P^2);indr=indr+A*K1; Info[indr,indc]=Info[indc,indr]=t(tmpout)
  tmpout = -t(Itil)%*%(Itil*(1-pi)^2/P^2); indc=indr; Info[indr,indc]=tmpout
  -Info
}

Est.Marg = function(dat,bet.ini,alp1.ini,alp0.ini,alp.known=F){ ##alp1,alp0 is A x K-1, remove first category
  p = ncol(bet.ini); Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  K1 = ncol(alp1.ini); uni.S = unique(sort(c(Si))); s1 = uni.S[1]; uni.S = uni.S[-1] ## remove first category
  Z = matrix(0,nrow=p*nn,ncol=p*2); for(jj in 1:p){Z[(1:nn)+nn*(jj-1),1:2+2*(jj-1)]=cbind(1,Gi[,jj])}
  Z = Z[rep(1:nrow(Z),A),] ## (nn*p)*A x 2p matrix
  Dmat = Imat = matrix(0,nrow=nn*p*A,ncol=K*A);
  Itil = matrix(0,nrow=nn*p*A,ncol=K1*A); eps = 1; step=0; I1 = rep(0,nn*p*A)
  for(aa in 1:A){tmpI = 1*(Si[,aa]==VTM(c(s1,uni.S),nn)); Imat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=tmpI[rep(1:nn,p),]
  tmpI = tmpI[,-1,drop=F]- 1*(Si[,aa]==s1);      Itil[1:(nn*p)+(aa-1)*(nn*p),1:K1+(aa-1)*K1]=tmpI[rep(1:nn,p),]
  I1[1:(nn*p)+(aa-1)*(nn*p)]=1*(Si[rep(1:nn,p),aa]==s1);Dmat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=1}
  bet = c(bet.ini); alp1 = c(t(alp1.ini)); alp0 = c(t(alp0.ini)); theta.old = c(bet,alp1,alp0)
  while(eps > 1e-5){
    eps2 = 1
    while(eps2 > 1e-5){
      betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
      P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
      d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
      score = c("S.bet"=t(d1P/P)%*%Z); Info = -t(Z)%*%(Z*d1P^2/P^2) ##Info = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2);
      bet.new = bet - c(solve(Info)%*%score); eps2 = sqrt(sum((bet.new-bet)^2)); bet = bet.new}
    Yihat.ja = g(betG+log((I1+Itil%*%alp1)/(I1+Itil%*%alp0))) ## n x A matrix
    if(!alp.known){
      alp1.new = c(matrix(c(t(Yihat.ja)%*%(Imat))/c(t(Yihat.ja)%*%Dmat),ncol=A)[-1,]);
      alp0.new = c(matrix(c(t(1-Yihat.ja)%*%Imat)/c(t(1-Yihat.ja)%*%Dmat),ncol=A)[-1,])
    }else{alp1.new = alp1; alp0.new = alp0}
    score = c("S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
    eps = sqrt(sum(c(alp1.new-alp1,alp0.new-alp0)^2));
    alp1 = alp1.new; alp0 = alp0.new; ##print(c(eps,score)); print(c(bet,alp1,alp0))
    #print(c(step,eps)); print(bet); print(cbind(matrix(alp1,nrow=A,byrow=T),matrix(alp0,nrow=A,byrow=T)))
    #if(step>2000){if(eps.vec[pmin(step,5000)]>eps.vec[pmin(step,5000)-200]){break}}
  }
  alp1 = t(matrix(alp1,ncol=A)); alp0 = t(matrix(alp0,ncol=A))
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1+1:K1),paste0("SP",1+1:K1))
  bet = matrix(bet,ncol=p); if(eps > 1e-5){bet=bet*NA;alp=alp*NA}
  list("bet"=bet,"alp"=alp)
}

Loglik.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  neglik = -sum(log(P)); print(neglik)
  neglik
}

Score.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
  -c("S.bet"=t(d1P/P)%*%Z,"S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
}

Info.Marg = function(theta, Itil,I1,Z,p,A,K1){ ##alp1,alp0 is A x K-1, remove first category
  bet = theta[1:(p*2)]; alp1 = theta[1:(A*K1)+2*p]; alp0 = theta[1:(A*K1)+2*p+A*K1]
  betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
  P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
  d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
  Info = matrix(0,ncol=p*2+A*K1*2,nrow=p*2+A*K1*2)
  tmpout = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2); indr=indc=1:(2*p); Info[indr,indc]=tmpout
  tmpout = -t(Z)%*%(Itil*(pi*d1P-P*d1pi)/P^2);indr=2*p+1:(A*K1); Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
  tmpout = -t(Z)%*%(Itil*((1-pi)*d1P+P*d1pi)/P^2);indr=indr+A*K1;Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
  tmpout = -t(Itil)%*%(Itil*pi^2/P^2); indr=indc=2*p+1:(A*K1); Info[indr,indc]=tmpout
  tmpout = -t(Itil)%*%(Itil*pi*(1-pi)/P^2);indr=indr+A*K1; Info[indr,indc]=Info[indc,indr]=t(tmpout)
  tmpout = -t(Itil)%*%(Itil*(1-pi)^2/P^2); indc=indr; Info[indr,indc]=tmpout
  -Info
}

Est.Marg.NR = function(dat,bet.ini,alp1.ini,alp0.ini){ ##alp1,alp0 is A x K-1, remove first category
  p = ncol(bet.ini); Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  K1 = ncol(alp1.ini); uni.S = unique(sort(c(Si))); s1 = uni.S[1]; uni.S = uni.S[-1] ## remove first category
  Z = matrix(0,nrow=p*nn,ncol=p*2); for(jj in 1:p){Z[(1:nn)+nn*(jj-1),1:2+2*(jj-1)]=cbind(1,Gi[,jj])}
  Z = Z[rep(1:nrow(Z),A),] ## (nn*p)*A x 2p matrix
  Itil = matrix(0,nrow=nn*p*A,ncol=K1*A); eps = 1; step=0; I1 = rep(0,nn*p*A)
  for(aa in 1:A){tmpI = 1*(Si[,aa]==VTM(uni.S,nn))-1*(Si[,aa]==s1);
  Itil[1:(nn*p)+(aa-1)*(nn*p),1:K1+(aa-1)*K1]=tmpI[rep(1:nn,p),]
  I1[1:(nn*p)+(aa-1)*(nn*p)]=1*(Si[rep(1:nn,p),aa]==s1)}
  bet = c(bet.ini); alp1 = c(t(alp1.ini)); alp0 = c(t(alp0.ini)); theta.old = c(bet,alp1,alp0)
  while(eps > 1e-5){
    bet = theta.old[1:(2*p)]; alp1 = theta.old[2*p+1:(A*K1)]; alp0 = theta.old[2*p+A*K1+1:(A*K1)]
    betG = Z%*%bet; pi = c(g(betG)); d1pi = c(d1g(betG)); d2pi = c(d2g(betG))
    P = I1 + pi*c(Itil%*%alp1)+(1-pi)*c(Itil%*%alp0)
    d1P = c(Itil%*%(alp1-alp0)); d2P = d2pi*d1P; d1P = d1pi*d1P
    score = c("S.bet"=t(d1P/P)%*%Z,"S.alp1"=t(pi/P)%*%Itil,"S.alp0"=t((1-pi)/P)%*%Itil)
    Info = matrix(0,ncol=p*2+A*K1*2,nrow=p*2+A*K1*2)
    tmpout = -t(Z)%*%(Z*(d1P^2-P*d2P)/P^2); indr=indc=1:(2*p); Info[indr,indc]=tmpout
    tmpout = -t(Z)%*%(Itil*(pi*d1P-P*d1pi)/P^2);indr=2*p+1:(A*K1); Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
    tmpout = -t(Z)%*%(Itil*((1-pi)*d1P+P*d1pi)/P^2);indr=indr+A*K1;Info[indr,indc]=t(tmpout); Info[indc,indr]=tmpout
    tmpout = -t(Itil)%*%(Itil*pi^2/P^2); indr=indc=2*p+1:(A*K1); Info[indr,indc]=tmpout
    tmpout = -t(Itil)%*%(Itil*pi*(1-pi)/P^2);indr=indr+A*K1; Info[indr,indc]=Info[indc,indr]=t(tmpout)
    tmpout = -t(Itil)%*%(Itil*(1-pi)^2/P^2); indc=indr; Info[indr,indc]=tmpout
    theta.new = theta.old - c(solve(Info)%*%score)
    eps = sqrt(sum((theta.new-theta.old)^2)); theta.old = theta.new; print(eps)
  }
  alp1 = t(matrix(alp1,ncol=A)); alp0 = t(matrix(alp0,ncol=A))
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1+1:K1),paste0("SP",1+1:K1))
  bet = matrix(bet,ncol=p); if(eps > 1e-5){bet=bet*NA;alp=alp*NA}
  list("bet"=bet,"alp"=alp)
}

Est.Marg.EM = function(dat,bet.ini,alp.ini,alp.known=F){
  p = ncol(bet.ini); Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  nn = nrow(dat); A = ncol(Si)
  K = ncol(alp.ini)/2; alp1.ini = alp.ini[,1:K,drop=F]; alp0.ini = alp.ini[,-(1:K),drop=F]
  uni.S = unique(sort(c(Si))); ## (nn*A x K)
  Z = matrix(0,nrow=p*nn,ncol=p*2); for(jj in 1:p){Z[(1:nn)+nn*(jj-1),1:2+2*(jj-1)]=cbind(1,Gi[,jj])}
  Z = Z[rep(1:nrow(Z),A),] ## (nn*p)*A x p matrix
  Dmat = Imat = matrix(0,nrow=nn*p*A,ncol=K*A); eps = 1; step=0
  for(aa in 1:A){tmpI = 1*(Si[,aa]==VTM(uni.S,nn)); Imat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=tmpI[rep(1:nn,p),]
  Dmat[1:(nn*p)+(aa-1)*(nn*p),1:K+(aa-1)*K]=1}
  bet = c(bet.ini); alp1 = c(t(alp1.ini)); alp0 = c(t(alp0.ini)); eps.vec = rep(0,5000)
  ## bet: intercept, bet (2 x p matrix); alp: A x K matrix
  while(eps > 1e-5){
    alp1.ini = alp1; alp0.ini = alp0; bet.ini = bet
    Yihat.ja = g(Z%*%bet+log((Imat%*%alp1)/(Imat%*%alp0))) ## n x A matrix
    bet = glm(Yihat.ja~Z-1,family=binomial)$coef
    Yihat.ja = g(Z%*%bet+log((Imat%*%alp1)/(Imat%*%alp0))) ## n x A matrix
    if(!alp.known){
      alp1 = c(t(Yihat.ja)%*%Imat)/c(t(Yihat.ja)%*%Dmat);
      alp0 = c(t(1-Yihat.ja)%*%Imat)/c(t(1-Yihat.ja)%*%Dmat)}
    eps = mean(abs(c(bet-bet.ini,alp1-alp1.ini,alp0-alp0.ini))); step=step+1; eps.vec[step]=eps
    print(c(step,eps)); print(bet); print(cbind(matrix(alp1,nrow=A,byrow=T),matrix(alp0,nrow=A,byrow=T)))
    if(step>2000){if(eps.vec[pmin(step,5000)]>eps.vec[pmin(step,5000)-200]){break}}
  }
  alp1 = t(matrix(alp1,ncol=A)); alp0 = t(matrix(alp0,ncol=A))
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1:K),paste0("SP",1:K))
  bet = matrix(bet,ncol=p); if(eps > 1e-5){bet=bet*NA;alp=alp*NA}
  list("bet"=bet,"alp"=alp)
}


Est.Marg.FUN.old = function(dat,bet.ini,alp.ini,alp.known=F){
  p = ncol(bet.ini); Gi = dat[,1:p]; Si = dat[,-(1:p),drop=F];  eps = 1; nn = nrow(dat); A = ncol(Si)
  K = ncol(alp.ini)/2; alp1.ini = alp.ini[,1:K,drop=F]; alp0.ini = alp.ini[,-(1:K),drop=F]
  uni.S = unique(sort(c(Si))); IiA = 1*(c(Si) == VTM(uni.S, length(c(Si)))) ## (nn*A x K)
  bet = bet.ini; alp1 = alp1.ini; alp0 = alp0.ini   ## bet: intercept, bet (2 x p matrix); alp: A x K matrix
  while(eps > 1e-5){
    alp1.ini = alp1; alp0.ini = alp0; bet.ini = bet
    alp1m = alp1[rep(1:A,rep(nn,A)),]; alp0m = alp0[rep(1:A,rep(nn,A)),]
    term1 = Gi*VTM(bet[2,],nn)+VTM(bet[1,],nn); ## n x p
    term2 = log(matrix((IiA*alp1m)%*%rep(1,K),nrow=nn)/matrix((IiA*alp0m)%*%rep(1,K),nrow=nn)) ## n x A matrix
    Yihat.j = matrix(apply(matrix(g(term1[,rep(1:p,A)]+term2[,rep(1:A,rep(p,A))]),nrow=nn*p),1,mean),nrow=nn)
    Yihat.a = matrix(apply(matrix(g(term1[,rep(1:p,rep(A,p))]+term2[,rep(1:A,p)]),nrow=nn*A),1,mean),nrow=nn)
    bet = sapply(1:p,function(j,Y,G){glm(Y[,j]~G[,j],family=binomial)$coef},Y=Yihat.j,G=Gi)
    if(!alp.known){
      alp1 = matrix(t(rep(1,nn))%*%(matrix(IiA,nrow=nn)*   Yihat.a[,rep(1:A,K)]), nrow=A)/apply(Yihat.a,  2,sum)
      alp0 = matrix(t(rep(1,nn))%*%(matrix(IiA,nrow=nn)*(1-Yihat.a[,rep(1:A,K)])),nrow=A)/apply(1-Yihat.a,2,sum)}
    eps = mean(abs(c(bet-bet.ini,alp1-alp1.ini,alp0-alp0.ini))); ##print(eps); print(bet); print(cbind(alp1,alp0))
  }
  alp = cbind(alp1,alp0); row.names(alp) = paste0("S",1:A); colnames(alp) = c(paste0("SE",1:K),paste0("SP",1:K))
  list("bet"=bet,"alp"=alp)
}


