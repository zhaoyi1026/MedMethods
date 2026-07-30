##############################################################################
# MedMethods method module: hdexposure
# Mediation with high-dimensional exposures and mediators
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from HDExposureMediator/HDExposureMediator-main/HD-CauseMediation.R ----
#################################
# Mediation
# High-dimensional exposure
# High-dimensional mediator
#################################

#################################
# soft-thresholding function
soft.thred<-function(mu,lambda)
{
  if(length(mu)==1)
  {
    return(sign(mu)*max(c(abs(mu)-lambda,0)))
  }else
  {
    re<-rep(NA,length(mu))
    for(j in 1:length(mu))
    {
      re[j]<-sign(mu[j])*max(c(abs(mu[j])-lambda,0))
    }
    return(re)
  }
}

# sparse group lasso solution function
# special case for our solution
sg.lasso.iden<-function(z,lambda1,lambda2)
{
  s<-soft.thred(z,lambda2)
  s.norm<-sqrt(sum(s^2))
  
  if(s.norm==0)
  {
    re<-rep(0,length(z))
  }else
  {
    re<-max(c(s.norm-lambda1,0))*s/s.norm 
  }
  
  return(re)
}
#################################

#################################
# penalty functions

# R1: pathway lasso penalty
pen.R1<-function(alpha,beta,gamma,phi=2,delta=0.1)
{
  q<-nrow(alpha)
  p<-ncol(alpha)
  
  mu<-alpha%*%diag(beta)
  
  s1<-sum(apply(mu,1,function(x){return(sum(abs(x)))}))
  s2<-phi*sum(apply(alpha,1,function(x){return(sum(x^2))}))
  s3<-phi*q*sum(beta^2)
  s4<-delta*sum(apply(alpha,1,function(x){return(sum(abs(x)))}))
  s5<-delta*sum(abs(beta))
  s6<-sum(abs(gamma))
  
  R1.sum<-s1+s2+s3+s4+s5+s6
  
  return(R1.sum)
}

pen.R2<-function(alpha,beta,gamma)
{
  q<-nrow(alpha)
  p<-ncol(alpha)
  
  mu<-alpha%*%diag(beta)
  
  s1<-sqrt(p)*sum(apply(mu,1,function(x){return(sqrt(sum(x^2)))}))
  s2<-sum(abs(gamma))
  
  R2.sum<-s1+s2
  
  return(R2.sum)
}

obj.func<-function(X,M,Y,alpha,beta,gamma,lambda=0,pi=1,phi=2,delta=0.1)
{
  ll<-(sum(diag(t(M-X%*%alpha)%*%(M-X%*%alpha)))+sum((Y-X%*%gamma-M%*%beta)^2))/2
  pen1<-pi*lambda*pen.R1(alpha,beta,gamma,phi,delta)
  pen2<-(1-pi)*lambda*pen.R2(alpha,beta,gamma)
  
  re<-data.frame(logLik=ll,R1=pen1,R2=pen2,obj=ll+pen1+pen2)
  return(re)
}
#################################

#################################
# HD-Cause and HD-Mediation function

# On standardized data
HDCauseMediation.std<-function(X,M,Y,lambda=0,pi=1,phi=2,delta=0.1,rho=1,max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,alpha0=NULL,beta0=NULL,gamma0=NULL)
{
  # input: 
  # data: X, M, Y
  # tuning parameters: lambda, pi, phi, delta
  # Laplacian parameter: rho
  
  n<-nrow(X)
  q<-ncol(X)
  p<-ncol(M)
  
  #---------------------------------------------
  # initial setups
  if(is.null(colnames(X)))
  {
    colnames(X)<-paste0("X",1:q)
  }
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL
    gamma.trace<-NULL
    mu.trace<-NULL
    tau.trace<-NULL
    
    objval.trace<-NULL
  }
  
  if(rho.increase)
  {
    rho0<-rho
  }else
  {
    rho0<-0
  }
  
  # set initial values
  if(is.null(alpha0))
  {
    alpha0<-matrix(0,q,p)
  }
  if(is.null(beta0))
  {
    beta0<-rep(0,p)
  }
  if(is.null(gamma0))
  {
    gamma0<-rep(0,q)
  }
  
  tau<-matrix(0,q,p)
  #---------------------------------------------
  
  #---------------------------------------------
  # optimization
  s<-0
  diff<-100
  
  time<-system.time(
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # update mu
      z<-alpha0%*%diag(beta0)-tau/rho
      lambda1<-(1-pi)*lambda*sqrt(p)/rho
      lambda2<-pi*lambda/rho
      mu.new<-matrix(NA,q,p)
      for(j in 1:q)
      {
        mu.new[j,]<-sg.lasso.iden(z[j,],lambda1,lambda2)
      }
      
      # update alpha
      alpha.new<-matrix(NA,q,p)
      for(j in 1:q)
      {
        V<-rho*diag(beta0)%*%diag(beta0)+(2*pi*lambda*phi+(t(X[,j])%*%X[,j])[1,1])*diag(rep(1,p))
        w<-t(M-X[,-j]%*%alpha0[-j,])%*%X[,j]+diag(beta0)%*%tau[j,]+rho*diag(beta0)%*%mu.new[j,]
        alpha.new[j,]<-solve(V)%*%soft.thred(w,pi*lambda*delta)
      }
      
      # update beta
      Vbeta<-t(M)%*%M+2*pi*lambda*phi*q*diag(rep(1,p))
      wbeta<-t(M)%*%(Y-X%*%gamma0)
      for(j in 1:q)
      {
        Vbeta<-Vbeta+rho*(diag(alpha.new[j,])%*%diag(alpha.new[j,]))
        wbeta<-wbeta+diag(alpha.new[j,])%*%tau[j,]+rho*diag(alpha.new[j,])%*%mu.new[j,]
      }
      beta.new<-c(solve(Vbeta)%*%soft.thred(wbeta,pi*lambda*delta))
      
      # update gamma
      Vgamma<-t(X)%*%X
      wgamma<-t(X)%*%(Y-M%*%beta.new)
      gamma.new<-c(solve(Vgamma)%*%soft.thred(wgamma,lambda))
      
      # update tau
      tau<-tau+rho*(mu.new-alpha.new%*%diag(beta.new))
      
      # objective function
      objval<-obj.func(X,M,Y,alpha.new,beta.new,gamma.new,lambda,pi,phi,delta)
      
      if(trace)
      {
        alpha.trace<-cbind(alpha.trace,alpha.new)
        beta.trace<-cbind(beta.trace,beta.new)
        gamma.trace<-cbind(gamma.trace,gamma.new)
        mu.trace<-cbind(mu.trace,mu.new)
        tau.trace<-cbind(tau.trace,tau)
        
        objval.trace<-c(objval.trace,objval$obj)
      }
      
      diff.alpha<-max(abs(alpha.new-alpha0))
      diff.beta<-max(abs(beta.new-beta0))
      diff.gamma<-max(abs(gamma.new-gamma0))
      
      diff<-max(c(diff.alpha,diff.beta,diff.gamma))
      
      # print(data.frame(diff=diff,alpha=diff.alpha,beta=diff.beta,gamma=diff.gamma,obj=objval$obj))
      
      # update parameters
      alpha0<-alpha.new
      beta0<-beta.new
      gamma0<-gamma.new
      
      rho<-rho+rho0
    })
  
  if(s>max.itr)
  {
    warning("Method does not converge!")
  }
  
  alpha.est<-alpha.new
  rownames(alpha.est)<-colnames(X)
  colnames(alpha.est)<-colnames(M)
  beta.est<-beta.new
  names(beta.est)<-colnames(M)
  gamma.est<-gamma.new
  names(gamma.est)<-colnames(X)
  IE.est<-mu.new
  rownames(IE.est)<-colnames(X)
  colnames(IE.est)<-colnames(M)
  #---------------------------------------------
  
  #---------------------------------------------
  if(trace)
  {
    re<-list(lambda=lambda,IE=IE.est,alpha=alpha.est,beta=beta.est,gamma=gamma.est,logLik=objval,converge=(s<=max.itr),time=time,
             alpha.trace=alpha.trace,beta.trace=beta.trace,gamma.trace=gamma.trace,mu.trace=mu.trace,obj.trace=objval.trace)
  }else
  {
    re<-list(lambda=lambda,IE=IE.est,alpha=alpha.est,beta=beta.est,gamma=gamma.est,logLik=objval,converge=(s<=max.itr),time=time)
  }
  return(re)
  #---------------------------------------------
}

# On original scale data
HDCauseMediation<-function(X,M,Y,lambda=0,pi=1,phi=2,delta=0.1,rho=1,standardize=TRUE,max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,alpha0=NULL,beta0=NULL,gamma0=NULL)
{
  # input: 
  # data: X, M, Y
  # tuning parameters: lambda, pi, phi, delta
  # Laplacian parameter: rho
  
  n<-nrow(X)
  q<-ncol(X)
  p<-ncol(M)
  
  if(standardize)
  {
    # standardize data
    X.sd<-apply(X,2,sd)
    M.sd<-apply(M,2,sd)
    Y.sd<-sd(Y)
    
    X.std<-scale(X,center=TRUE,scale=TRUE)
    M.std<-scale(M,center=TRUE,scale=TRUE)
    Y.std<-scale(Y,center=TRUE,scale=TRUE)
    
    re.std<-HDCauseMediation.std(X.std,M.std,Y.std,lambda=lambda,pi=pi,phi=phi,delta=delta,rho=rho,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,alpha0=alpha0,beta0=beta0,gamma0=gamma0)
    
    alpha.est<-matrix(NA,q,p)
    rownames(alpha.est)<-rownames(re.std$alpha)
    colnames(alpha.est)<-colnames(re.std$alpha)
    IE.est<-matrix(NA,q,p)
    rownames(IE.est)<-rownames(re.std$IE)
    colnames(IE.est)<-colnames(re.std$IE)
    for(j in 1:q)
    {
      for(k in 1:p)
      {
        alpha.est[j,k]<-re.std$alpha[j,k]*(M.sd[k]/X.sd[j])
        IE.est[j,k]<-re.std$IE[j,k]*(Y.sd/X.sd[j])
      }
    }
    beta.est<-re.std$beta*(Y.sd/M.sd)
    gamma.est<-re.std$gamma*(Y.sd/X.sd)
    
    objval<-obj.func(X,M,Y,alpha.est,beta.est,gamma.est,lambda,pi,phi,delta)
    
    re<-list(lambda=lambda,IE=IE.est,alpha=alpha.est,beta=beta.est,gamma=gamma.est,logLik=objval,converge=re.std$converge,time=re.std$time,out.scaled=re.std)
    
  }else
  {
    re<-HDCauseMediation.std(X,M,Y,lambda=lambda,pi=pi,phi=phi,delta=delta,rho=rho,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,alpha0=alpha0,beta0=beta0,gamma0=gamma0)
  }
  
  return(re)
}
#################################

#################################
# PCA on X first
HDCauseMediationPCA<-function(X,M,Y,adaptive=FALSE,var.prop=0.9,n.pc=NULL,lambda=0,pi=1,phi=2,delta=0.1,rho=1,standardize=TRUE,
                              max.itr=10000,tol=1e-6,rho.increase=FALSE,trace=FALSE,alpha0=NULL,beta0=NULL,gamma0=NULL)
{
  # input: 
  # data: X, M, Y
  # tuning parameters: lambda, pi, phi, delta
  # Laplacian parameter: rho
  
  n<-nrow(X)
  r<-ncol(X)
  p<-ncol(M)
  
  if(is.null(colnames(X)))
  {
    colnames(X)<-paste0("X",1:r)
  }
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  # step 1: PCA on X
  Sigma.X<-cov(X)
  svd.X<-svd(Sigma.X)
  if(adaptive)
  {
    n.pc<-max(sum(cumsum(svd.X$d)/sum(svd.X$d)<var.prop)+1,4)
  }else
  {
    if(is.null(n.pc))
    {
      if(r>20)
      {
        n.pc<-max(sum(cumsum((svd.X$d))/sum((svd.X$d))<var.prop)+1,4)
      }else
      {
        n.pc<-r
      }
    }
  }
  Gamma.est<-svd.X$u[,1:n.pc]
  colnames(Gamma.est)<-paste0("PC",1:n.pc)
  rownames(Gamma.est)<-colnames(X)
  for(k in 1:n.pc)
  {
    if(Gamma.est[which.max(abs(Gamma.est[,k])),k]<0)
    {
      Gamma.est[,k]<--Gamma.est[,k]
    }
  }
  X.t<-X%*%Gamma.est
  colnames(X.t)<-paste0("XPC",1:n.pc)
  
  # mediation
  re.med<-HDCauseMediation(X.t,M,Y,lambda=lambda,pi=pi,phi=phi,delta=delta,rho=rho,standardize=standardize,max.itr=max.itr,tol=tol,rho.increase=rho.increase,trace=trace,alpha0=alpha0,beta0=beta0,gamma0=gamma0)
  
  # X outcome
  X.varper<-rbind((svd.X$d[1:n.pc])^2/sum((svd.X$d)^2),cumsum((svd.X$d[1:n.pc])^2/sum((svd.X$d)^2)))
  colnames(X.varper)<-paste0("PC",1:n.pc)
  rownames(X.varper)<-c("Per","CumPer")
  re.X<-list(X.nPC=n.pc,X.varper=X.varper,X.Gamma=Gamma.est)
  
  re<-c(re.X,re.med)
  
  return(re)
}
#################################

