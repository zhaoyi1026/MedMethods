##############################################################################
# MedMethods method module: heteromed
# Heterogeneous mediation effects
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from hetero_mediation/V2/HeterMed.R ----
#############################
# method packages

# [MedMethods] removed at assembly: library(genlasso)
# [MedMethods] removed at assembly: library(MASS)
#############################

#############################
# methods functions: estimation
med.inter.ITE<-function(X,Z,alpha0,alpha1,beta0,beta1,gamma0,gamma1)
{
  n<-length(X)
  
  # calculate NIE and NDE
  otmp<-matrix(NA,n,3)
  colnames(otmp)<-c("Treatment","NIE","NDE")
  otmp[,1]<-X
  otmp[,2]<-2*(beta0+beta1*X)*(Z%*%alpha1)
  otmp[,3]<-2*Z%*%gamma1+2*beta1*(Z%*%alpha0+apply(Z,2,function(x){return(x*X)})%*%alpha1)
  
  return(otmp)
}

genlasso.tune<-function(y,out.gen)
{
  n<-length(y)
  
  n.out<-ncol(out.gen$beta)
  
  out.res<-apply(out.gen$fit,2,function(x){return(y-x)})
  rss<-colSums(out.res^2)
  
  if(n-10>out.gen$df[n.out])
  {
    sigma2<-rss[n.out]/(n-out.gen$df[n.out])
  }else
  {
    sigma2<-rss[n.out]/n
  }
  
  # Cp
  Cp.vec<-(rss/sigma2)+2*out.gen$df-n
  
  # BIC
  BIC.vec<-n*log(rss/n)+out.gen$df*log(n)
  
  # GCV
  GCV.vec<-rss/((n-out.gen$df)^2)
  
  re<-cbind(lambda=out.gen$lambda,Cp=Cp.vec,BIC=BIC.vec,GCV=GCV.vec)
  colnames(re)<-c("lambda","Cp","BIC","GCV")
  
  return(re)
}

med.inter<-function(X,M,Y,Z,method=c("OLS","genlasso"),genlasso.tune.method=c("Cp","BIC","GCV"))
{
  # X: n by 1 treatment vector
  # M: n by 1 mediator vector
  # Y: n by 1 outcome vector
  # Z: n by p covariates, first column ones
  
  p<-ncol(Z)
  n<-length(X)
  
  if(is.null(colnames(Z))==TRUE)
  {
    colnames(Z)<-c("Intercept",paste0("Z",1:(p-1)))
  }
  
  
  if(method[1]=="OLS")
  {
    idx1<-which(X==1)
    idx2<-which(X==-1)
    
    # M model
    ytmp<-c(M[idx1],M[idx2])
    xtmp<-rbind(cbind(Z[idx1,],matrix(0,length(idx1),p)),cbind(matrix(0,length(idx2),p),Z[idx2,]))
    m.theta<-solve(t(xtmp)%*%xtmp)%*%(t(xtmp)%*%ytmp)
    m.theta1<-m.theta[1:p]
    m.theta2<-m.theta[(p+1):(2*p)]
    alpha0.est<-(m.theta1+m.theta2)/2
    alpha1.est<-(m.theta1-m.theta2)/2
    
    # Y model
    W<-cbind(Z,M)
    colnames(W)<-c(colnames(Z),"M")
    ytmp<-c(Y[idx1],Y[idx2])
    xtmp<-rbind(cbind(W[idx1,],matrix(0,length(idx1),p+1)),cbind(matrix(0,length(idx2),p+1),W[idx2,]))
    y.theta<-solve(t(xtmp)%*%xtmp)%*%(t(xtmp)%*%ytmp)
    y.theta1<-y.theta[1:(p+1)]
    y.theta2<-y.theta[(p+2):(2*(p+1))]
    y.phi0<-(y.theta1+y.theta2)/2
    y.phi1<-(y.theta1-y.theta2)/2
    gamma0.est<-y.phi0[1:p]
    beta0.est<-y.phi0[p+1]
    gamma1.est<-y.phi1[1:p]
    beta1.est<-y.phi1[p+1]
  }
  if(method[1]=="genlasso")
  {
    idx1<-which(X==1)
    idx2<-which(X==-1)
    
    # M model
    ytmp<-c(M[idx1],M[idx2])
    xtmp<-rbind(cbind(Z[idx1,],matrix(0,length(idx1),p)),cbind(matrix(0,length(idx2),p),Z[idx2,]))
    m.D.base<-cbind(rep(0,p-1),diag(rep(1,p-1)))
    m.D<-rbind(cbind(m.D.base,m.D.base),cbind(m.D.base,m.D.base*-1))
    out.gen<-genlasso(y=ytmp,X=xtmp,D=m.D)
    tune.tmp<-genlasso.tune(y=ytmp,out.gen=out.gen)
    lambda.idx<-which.min(tune.tmp[,genlasso.tune.method[1]])
    m.theta<-out.gen$beta[,lambda.idx]
    m.theta1<-m.theta[1:p]
    m.theta2<-m.theta[(p+1):(2*p)]
    alpha0.est<-(m.theta1+m.theta2)/2
    alpha1.est<-(m.theta1-m.theta2)/2
    
    # Y model
    W<-cbind(Z,M)
    colnames(W)<-c(colnames(Z),"M")
    ytmp<-c(Y[idx1],Y[idx2])
    xtmp<-rbind(cbind(W[idx1,],matrix(0,length(idx1),p+1)),cbind(matrix(0,length(idx2),p+1),W[idx2,]))
    y.D.base<-cbind(rep(0,p-1),diag(rep(1,p-1)),rep(0,p-1))
    y.D<-rbind(cbind(y.D.base,y.D.base),cbind(y.D.base,y.D.base*-1))
    out.gen<-genlasso(y=ytmp,X=xtmp,D=y.D)
    tune.tmp<-genlasso.tune(y=ytmp,out.gen=out.gen)
    lambda.idx<-which.min(tune.tmp[,genlasso.tune.method[1]])
    y.theta<-out.gen$beta[,lambda.idx]
    y.theta1<-y.theta[1:(p+1)]
    y.theta2<-y.theta[(p+2):(2*(p+1))]
    y.phi0<-(y.theta1+y.theta2)/2
    y.phi1<-(y.theta1-y.theta2)/2
    gamma0.est<-y.phi0[1:p]
    beta0.est<-y.phi0[p+1]
    gamma1.est<-y.phi1[1:p]
    beta1.est<-y.phi1[p+1]
  }
  
  # calculate NIE and NDE
  otmp<-med.inter.ITE(X=X,Z=Z,alpha0=alpha0.est,alpha1=alpha1.est,beta0=beta0.est,beta1=beta1.est,gamma0=gamma0.est,gamma1=gamma1.est)
  
  re<-list(ITE=otmp,alpha0=alpha0.est,alpha1=alpha1.est,beta0=beta0.est,beta1=beta1.est,gamma0=gamma0.est,gamma1=gamma1.est)
  
  return(re)
}
#############################

#############################
# inference

fit.inf.OLS<-function(X,M,Y,Z,out.med.inter,conf.level=0.95)
{
  # X: n by 1 treatment vector
  # M: n by 1 mediator vector
  # Y: n by 1 outcome vector
  # Z: n by p covariates, first column ones
  # out.med.inter: output from med.inter() function with method="OLS"
  
  p<-ncol(Z)
  n<-length(X)
  
  M<-matrix(M,ncol=1)
  colnames(M)<-"M"
  Y<-matrix(Y,ncol=1)
  colnames(Y)<-"Y"
  
  if(is.null(colnames(Z))==TRUE)
  {
    colnames(Z)<-c("Intercept",paste0("Z",1:(p-1)))
  }
  
  U<-diag(X)%*%Z
  V<-diag(X)%*%M
  
  # estimate model error variance
  s2.m<-mean((M-Z%*%out.med.inter$alpha0-U%*%out.med.inter$alpha1)^2)
  s2.y<-mean((Y-Z%*%out.med.inter$gamma0-U%*%out.med.inter$gamma1-M%*%out.med.inter$beta0-V%*%out.med.inter$beta1)^2)
  
  # calculate design matrix quantities
  Q<-(t(Z)%*%Z)/n
  X.mm<-c(mean(X),mean(X^2),mean(X^3),mean(X^4))
  
  Qzm<-Q%*%out.med.inter$alpha0+X.mm[1]*(Q%*%out.med.inter$alpha1)
  Qzv<-X.mm[1]*(Q%*%out.med.inter$alpha0)+X.mm[2]*(Q%*%out.med.inter$alpha1)
  Qum<-X.mm[1]*(Q%*%out.med.inter$alpha0)+X.mm[2]*(Q%*%out.med.inter$alpha1)
  Quv<-X.mm[2]*(Q%*%out.med.inter$alpha0)+X.mm[3]*(Q%*%out.med.inter$alpha1)
  Qm<-(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha0)[1,1]+2*X.mm[1]*(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha1)[1,1]+
    X.mm[2]*(t(out.med.inter$alpha1)%*%Q%*%out.med.inter$alpha1)[1,1]+s2.m
  Qmv<-X.mm[1]*(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha0)[1,1]+2*X.mm[2]*(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha1)[1,1]+
    X.mm[3]*(t(out.med.inter$alpha1)%*%Q%*%out.med.inter$alpha1)[1,1]+X.mm[1]*s2.m
  Qv<-X.mm[2]*(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha0)[1,1]+2*X.mm[3]*(t(out.med.inter$alpha0)%*%Q%*%out.med.inter$alpha1)[1,1]+
    X.mm[4]*(t(out.med.inter$alpha1)%*%Q%*%out.med.inter$alpha1)[1,1]+X.mm[2]*s2.m
  
  Qx<-rbind(cbind(Q,X.mm[1]*Q,Qzm,Qzv),
            cbind(X.mm[1]*Q,X.mm[2]*Q,Qum,Quv),
            cbind(t(Qzm),t(Qum),Qm,Qmv),
            cbind(t(Qzv),t(Quv),t(Qmv),Qv))
  Qx.inv<-ginv(Qx)
  
  Theta<-cbind(c(out.med.inter$alpha0,out.med.inter$alpha1,0,0),c(out.med.inter$gamma0,out.med.inter$gamma1,out.med.inter$beta0,out.med.inter$beta1))
  Sigma<-diag(c(s2.m,s2.y))
  
  # covariance of vec(Theta)
  # Xi<-kronecker(Sigma,Qx)
  # Px<-kronecker(diag(rep(1,2)),Qx.inv)
  # cov.vecTheta<-(Px%*%Xi%*%Px)/n
  cov.vecTheta<-kronecker(Sigma,Qx.inv)/n
  # SE of each parameter
  vecTheta.se<-sqrt(diag(cov.vecTheta))
  Theta.se<-matrix(vecTheta.se,ncol=2)
  
  # organize output
  zv<-qnorm(1-(1-conf.level)/2)
  
  alpha0.out<-data.frame(Estimate=out.med.inter$alpha0,SE=Theta.se[1:p,1])
  alpha0.out$zvalue<-alpha0.out$Estimate/alpha0.out$SE
  alpha0.out$pvalue<-(1-pnorm(abs(alpha0.out$zvalue)))*2
  alpha0.out$LB<-alpha0.out$Estimate-zv*alpha0.out$SE
  alpha0.out$UB<-alpha0.out$Estimate+zv*alpha0.out$SE
  rownames(alpha0.out)<-colnames(Z)
  
  alpha1.out<-data.frame(Estimate=out.med.inter$alpha1,SE=Theta.se[(p+1):(2*p),1])
  alpha1.out$zvalue<-alpha1.out$Estimate/alpha1.out$SE
  alpha1.out$pvalue<-(1-pnorm(abs(alpha1.out$zvalue)))*2
  alpha1.out$LB<-alpha1.out$Estimate-zv*alpha1.out$SE
  alpha1.out$UB<-alpha1.out$Estimate+zv*alpha1.out$SE
  rownames(alpha1.out)<-colnames(Z)
  
  gamma0.out<-data.frame(Estimate=out.med.inter$gamma0,SE=Theta.se[1:p,2])
  gamma0.out$zvalue<-gamma0.out$Estimate/gamma0.out$SE
  gamma0.out$pvalue<-(1-pnorm(abs(gamma0.out$zvalue)))*2
  gamma0.out$LB<-gamma0.out$Estimate-zv*gamma0.out$SE
  gamma0.out$UB<-gamma0.out$Estimate+zv*gamma0.out$SE
  rownames(gamma0.out)<-colnames(Z)
  
  gamma1.out<-data.frame(Estimate=out.med.inter$gamma1,SE=Theta.se[(p+1):(2*p),2])
  gamma1.out$zvalue<-gamma1.out$Estimate/gamma1.out$SE
  gamma1.out$pvalue<-(1-pnorm(abs(gamma1.out$zvalue)))*2
  gamma1.out$LB<-gamma1.out$Estimate-zv*gamma1.out$SE
  gamma1.out$UB<-gamma1.out$Estimate+zv*gamma1.out$SE
  rownames(gamma1.out)<-colnames(Z)
  
  beta0.out<-data.frame(Estimate=out.med.inter$beta0,SE=Theta.se[2*p+1,2])
  beta0.out$zvalue<-beta0.out$Estimate/beta0.out$SE
  beta0.out$pvalue<-(1-pnorm(abs(beta0.out$zvalue)))*2
  beta0.out$LB<-beta0.out$Estimate-zv*beta0.out$SE
  beta0.out$UB<-beta0.out$Estimate+zv*beta0.out$SE
  rownames(beta0.out)<-"M"
  
  beta1.out<-data.frame(Estimate=out.med.inter$beta1,SE=Theta.se[2*p+2,2])
  beta1.out$zvalue<-beta1.out$Estimate/beta1.out$SE
  beta1.out$pvalue<-(1-pnorm(abs(beta1.out$zvalue)))*2
  beta1.out$LB<-beta1.out$Estimate-zv*beta1.out$SE
  beta1.out$UB<-beta1.out$Estimate+zv*beta1.out$SE
  rownames(beta1.out)<-"M"
  
  h1.func<-function(xtmp)
  {
    h1<-c(rep(0,p),2*(out.med.inter$beta0+out.med.inter$beta1*xtmp[1])*xtmp[-1],rep(0,2+2*p),2*t(out.med.inter$alpha1)%*%xtmp[-1],2*xtmp[1]*t(out.med.inter$alpha1)%*%xtmp[-1])
    return((t(h1)%*%cov.vecTheta%*%h1)[1,1])
  }
  h2.func<-function(xtmp)
  {
    h2<-c(2*out.med.inter$beta1*xtmp[-1],2*out.med.inter$beta1*xtmp[1]*xtmp[-1],rep(0,2+p),2*xtmp[-1],0,2*(t(out.med.inter$alpha0)%*%xtmp[-1]+xtmp[1]*(t(out.med.inter$alpha1)%*%xtmp[-1])))
    return((t(h2)%*%cov.vecTheta%*%h2)[1,1])
  }
  IE.var<-apply(cbind(X,Z),1,h1.func)
  DE.var<-apply(cbind(X,Z),1,h2.func)
  
  IE.out<-data.frame(Treatment=X,Estimate=out.med.inter$ITE[,"NIE"],SE=sqrt(IE.var))
  IE.out$LB<-IE.out$Estimate-zv*IE.out$SE
  IE.out$UB<-IE.out$Estimate+zv*IE.out$SE
  
  DE.out<-data.frame(Treatment=X,Estimate=out.med.inter$ITE[,"NDE"],SE=sqrt(DE.var))
  DE.out$LB<-DE.out$Estimate-zv*DE.out$SE
  DE.out$UB<-DE.out$Estimate+zv*DE.out$SE
  
  re<-list(NIE=IE.out,NDE=DE.out,alpha0=alpha0.out,alpha1=alpha1.out,beta0=beta0.out,beta1=beta1.out,gamma0=gamma0.out,gamma1=gamma1.out,Theta=Theta,vecTheta.cov=cov.vecTheta)
  
  return(re)
}

fit.inf.genlasso<-function(X,M,Y,Z,out.med.inter,conf.level=0.95,genlasso.tune.method=c("Cp","BIC","GCV"),B=500,zero.thred=1e-4,lambda.ridge=0.001,verbose=TRUE)
{
  # X: n by 1 treatment vector
  # M: n by 1 mediator vector
  # Y: n by 1 outcome vector
  # Z: n by p covariates, first column ones
  # # out.med.inter: output from med.inter() function with method="genlasso"
  
  p<-ncol(Z)
  n<-length(X)
  
  M<-matrix(M,ncol=1)
  colnames(M)<-"M"
  Y<-matrix(Y,ncol=1)
  colnames(Y)<-"Y"
  
  if(is.null(colnames(Z))==TRUE)
  {
    colnames(Z)<-c("Intercept",paste0("Z",1:(p-1)))
  }
  
  n.trt<-length(which(X==1))
  n.ctl<-length(which(X==-1))
  
  alpha0.est<-matrix(0,B,p)
  colnames(alpha0.est)<-colnames(Z)
  alpha1.est=gamma0.est=gamma1.est<-alpha0.est
  beta0.est=beta1.est<-rep(0,B)
  IE.est=DE.est<-matrix(NA,n,B)
  for(b in 1:B)
  {
    # training index (sample by treatment group)
    idx.tr.trt<-sample(which(X==1),round(n.trt/2),replace=FALSE)
    idx.tr.ctl<-sample(which(X==-1),round(n.ctl/2),replace=FALSE)
    idx.tr<-sort(c(idx.tr.trt,idx.tr.ctl))
    # testing index
    idx.ts<-(1:n)[-idx.tr]
    
    # training data
    X.tr<-X[idx.tr]
    M.tr<-matrix(M[idx.tr,],ncol=ncol(M))
    colnames(M.tr)<-colnames(M)
    Y.tr<-matrix(Y[idx.tr,],ncol=ncol(Y))
    colnames(Y.tr)<-colnames(Y)
    Z.tr<-matrix(Z[idx.tr,],ncol=ncol(Z))
    colnames(Z.tr)<-colnames(Z)
    
    # testing data
    X.ts<-X[idx.ts]
    M.ts<-matrix(M[idx.ts,],ncol=ncol(M))
    colnames(M.ts)<-colnames(M)
    Y.ts<-matrix(Y[idx.ts,],ncol=ncol(Y))
    colnames(Y.ts)<-colnames(Y)
    Z.ts<-matrix(Z[idx.ts,],ncol=ncol(Z))
    colnames(Z.ts)<-colnames(Z)
    
    # training data output
    out.tr<-med.inter(X.tr,M.tr,Y.tr,Z.tr,method="genlasso",genlasso.tune.method=genlasso.tune.method)
    idx.alpha0<-which(abs(out.tr$alpha0)>zero.thred)
    idx.alpha1<-which(abs(out.tr$alpha1)>zero.thred)
    idx.gamma0<-which(abs(out.tr$gamma0)>zero.thred)
    idx.gamma1<-which(abs(out.tr$gamma1)>zero.thred)
    
    # refit using testing
    U.ts<-diag(X.ts)%*%Z.ts
    V.ts<-diag(X.ts)%*%M.ts
    Mmat.ts<-cbind(Z.ts[,idx.alpha0],U.ts[,idx.alpha1])
    otmp.m<-c(ginv(t(Mmat.ts)%*%Mmat.ts+lambda.ridge*diag(rep(1,ncol(Mmat.ts))))%*%(t(Mmat.ts)%*%M.ts))
    alpha0.est[b,idx.alpha0]<-otmp.m[1:length(idx.alpha0)]
    alpha1.est[b,idx.alpha1]<-otmp.m[(length(idx.alpha0)+1):length(otmp.m)]
    Ymat.ts<-cbind(Z.ts[,idx.gamma0],U.ts[,idx.gamma1],M.ts,V.ts)
    otmp.y<-c(ginv(t(Ymat.ts)%*%Ymat.ts+lambda.ridge*diag(rep(1,ncol(Ymat.ts))))%*%(t(Ymat.ts)%*%Y.ts))
    gamma0.est[b,idx.gamma0]<-otmp.y[1:length(idx.gamma0)]
    gamma1.est[b,idx.gamma1]<-otmp.y[(length(idx.gamma0)+1):(length(idx.gamma0)+length(idx.gamma1))]
    beta0.est[b]<-otmp.y[length(idx.gamma0)+length(idx.gamma1)+1]
    beta1.est[b]<-otmp.y[length(idx.gamma0)+length(idx.gamma1)+2]
    
    # estimate IE and DE for each individual
    otmp.med<-med.inter.ITE(X=X,Z=Z,alpha0=alpha0.est[b,],alpha1=alpha1.est[b,],beta0=beta0.est[b],beta1=beta1.est[b],gamma0=gamma0.est[b,],gamma1=gamma1.est[b,])
    IE.est[,b]<-otmp.med[,"NIE"]
    DE.est[,b]<-otmp.med[,"NDE"]
    
    if(verbose)
    {
      print(paste0("Splitting sample ",b," done!"))
    }
  }
  
  # organize output
  zv<-qnorm(1-(1-conf.level)/2)
  
  alpha0.out<-data.frame(Estimate=apply(alpha0.est,2,mean,na.rm=TRUE),SE=apply(alpha0.est,2,sd,na.rm=TRUE))
  alpha0.out$zvalue<-alpha0.out$Estimate/alpha0.out$SE
  alpha0.out$pvalue<-(1-pnorm(abs(alpha0.out$zvalue)))*2
  alpha0.out$LB<-alpha0.out$Estimate-zv*alpha0.out$SE
  alpha0.out$UB<-alpha0.out$Estimate+zv*alpha0.out$SE
  rownames(alpha0.out)<-colnames(Z)
  
  alpha1.out<-data.frame(Estimate=apply(alpha1.est,2,mean,na.rm=TRUE),SE=apply(alpha1.est,2,sd,na.rm=TRUE))
  alpha1.out$zvalue<-alpha1.out$Estimate/alpha1.out$SE
  alpha1.out$pvalue<-(1-pnorm(abs(alpha1.out$zvalue)))*2
  alpha1.out$LB<-alpha1.out$Estimate-zv*alpha1.out$SE
  alpha1.out$UB<-alpha1.out$Estimate+zv*alpha1.out$SE
  rownames(alpha1.out)<-colnames(Z)
  
  gamma0.out<-data.frame(Estimate=apply(gamma0.est,2,mean,na.rm=TRUE),SE=apply(gamma0.est,2,sd,na.rm=TRUE))
  gamma0.out$zvalue<-gamma0.out$Estimate/gamma0.out$SE
  gamma0.out$pvalue<-(1-pnorm(abs(gamma0.out$zvalue)))*2
  gamma0.out$LB<-gamma0.out$Estimate-zv*gamma0.out$SE
  gamma0.out$UB<-gamma0.out$Estimate+zv*gamma0.out$SE
  rownames(gamma0.out)<-colnames(Z)
  
  gamma1.out<-data.frame(Estimate=apply(gamma1.est,2,mean,na.rm=TRUE),SE=apply(gamma1.est,2,sd,na.rm=TRUE))
  gamma1.out$zvalue<-gamma1.out$Estimate/gamma1.out$SE
  gamma1.out$pvalue<-(1-pnorm(abs(gamma1.out$zvalue)))*2
  gamma1.out$LB<-gamma1.out$Estimate-zv*gamma1.out$SE
  gamma1.out$UB<-gamma1.out$Estimate+zv*gamma1.out$SE
  rownames(gamma1.out)<-colnames(Z)
  
  beta0.out<-data.frame(Estimate=mean(beta0.est,na.rm=TRUE),SE=sd(beta0.est,na.rm=TRUE))
  beta0.out$zvalue<-beta0.out$Estimate/beta0.out$SE
  beta0.out$pvalue<-(1-pnorm(abs(beta0.out$zvalue)))*2
  beta0.out$LB<-beta0.out$Estimate-zv*beta0.out$SE
  beta0.out$UB<-beta0.out$Estimate+zv*beta0.out$SE
  rownames(beta0.out)<-"M"
  
  beta1.out<-data.frame(Estimate=mean(beta1.est,na.rm=TRUE),SE=sd(beta1.est,na.rm=TRUE))
  beta1.out$zvalue<-beta1.out$Estimate/beta1.out$SE
  beta1.out$pvalue<-(1-pnorm(abs(beta1.out$zvalue)))*2
  beta1.out$LB<-beta1.out$Estimate-zv*beta1.out$SE
  beta1.out$UB<-beta1.out$Estimate+zv*beta1.out$SE
  rownames(beta1.out)<-"M"
  
  IE.out<-data.frame(Treatment=X,Estimate=apply(IE.est,1,mean,na.rm=TRUE),SE=apply(IE.est,1,sd,na.rm=TRUE))
  IE.out$LB<-IE.out$Estimate-zv*IE.out$SE
  IE.out$UB<-IE.out$Estimate+zv*IE.out$SE
  
  DE.out<-data.frame(Treatment=X,Estimate=apply(DE.est,1,mean,na.rm=TRUE),SE=apply(DE.est,1,sd,na.rm=TRUE))
  DE.out$LB<-DE.out$Estimate-zv*DE.out$SE
  DE.out$UB<-DE.out$Estimate+zv*DE.out$SE
  
  re<-list(NIE=IE.out,NDE=DE.out,alpha0=alpha0.out,alpha1=alpha1.out,beta0=beta0.out,beta1=beta1.out,gamma0=gamma0.out,gamma1=gamma1.out)
  
  return(re)
}

med.inter.inf<-function(X,M,Y,Z,out.med.inter,conf.level=0.95,method=c("OLS","genlasso"),genlasso.tune.method=c("Cp","BIC","GCV"),B=500,zero.thred=1e-4,lambda.ridge=0.001,verbose=TRUE)
{
  if(method=="OLS")
  {
    out<-fit.inf.OLS(X,M,Y,Z,out.med.inter=out.med.inter,conf.level=conf.level)
  }
  if(method=="genlasso")
  {
    out<-fit.inf.genlasso(X,M,Y,Z,out.med.inter=out.med.inter,conf.level=conf.level,
                          genlasso.tune.method=genlasso.tune.method,B=B,zero.thred=zero.thred,lambda.ridge=lambda.ridge,verbose=verbose)
  }
  
  return(out)
}
#############################

