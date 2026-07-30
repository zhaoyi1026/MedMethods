##############################################################################
# MedMethods method module: macc
# Multilevel mediation analysis with structured unmeasured confounding
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from macc/macc/R/cma.delta.lm.HL.R ----
cma.delta.lm.HL <-
function(dat,delta,max.itr=500,tol=1e-4,error.indep=FALSE,error.var.equal=FALSE,Sigma.update=FALSE,
                          var.constraint=FALSE)
{
  re<-cma.delta.lm(dat,delta,max.itr,tol,error.indep,error.var.equal,Sigma.update,var.constraint)
  return(re$HL)
}

### ---- from macc/macc/R/cma.delta.lm.R ----
cma.delta.lm <-
function(dat,delta=0,max.itr=500,tol=1e-4,error.indep=FALSE,error.var.equal=FALSE,Sigma.update=FALSE,
                       var.constraint=FALSE)
{
  N<-length(dat)
  K<-1
  
  ######################################
  # Estimate A, B, C and C' for each subject
  At=Bt=Ct=C2t<-matrix(NA,N,K)
  sigma1=sigma2<-matrix(NA,N,K)
  for(i in 1:N)
  {
    dd<-dat[[i]]
    re<-cma.uni.delta(dd,delta=delta)
    
    At[i,1]<-re$D[1,1]
    Ct[i,1]<-re$D[1,2]
    Bt[i,1]<-re$D[2,2]
    C2t[i,1]<-re$Coefficients[4,1]
    
    sigma1[i]<-re$Sigma[1,1]
    sigma2[i]<-re$Sigma[2,2]
  }
  ######################################
  
  ######################################
  # Confidence interval of variance estimate
  fit.A<-gls(At~1)
  fit.B<-gls(Bt~1)
  fit.C<-gls(Ct~1)
  
  if(is.matrix(var.constraint))
  {
    if(nrow(var.constraint)==3)
    {
      Lambda.confint<-var.constraint[1:3,]
      colnames(Lambda.confint)<-c("LB","UB")
      rownames(Lambda.confint)<-c("A","B","C")
    }else
    {
      warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
      Lambda.confint<-matrix(NA,3,2)
      colnames(Lambda.confint)<-c("LB","UB")
      rownames(Lambda.confint)<-c("A","B","C")
      intA<-intervals(fit.A)
      intB<-intervals(fit.B)
      intC<-intervals(fit.C)
      Lambda.confint[1,]<-as.matrix(intA[[2]][c(1,3)])^2
      Lambda.confint[2,]<-as.matrix(intB[[2]][c(1,3)])^2
      Lambda.confint[3,]<-as.matrix(intC[[2]][c(1,3)])^2 
    }
  }else
    if(var.constraint)
    {
      Lambda.confint<-matrix(NA,3,2)
      colnames(Lambda.confint)<-c("LB","UB")
      rownames(Lambda.confint)<-c("A","B","C")
      intA<-intervals(fit.A)
      intB<-intervals(fit.B)
      intC<-intervals(fit.C)
      Lambda.confint[1,]<-as.matrix(intA[[2]][c(1,3)])^2
      Lambda.confint[2,]<-as.matrix(intB[[2]][c(1,3)])^2
      Lambda.confint[3,]<-as.matrix(intC[[2]][c(1,3)])^2 
    }else
    {
      Lambda.confint<-NULL
    }
  ######################################
  
  ######################################
  # Covariance estimate
  Lambda.hat<-matrix(0,3,3)
  if(!error.var.equal)
  {
    if(error.indep)
    {
      fit.A<-lm(At~1)
      fit.B<-lm(Bt~1)
      fit.C<-lm(Ct~1)
      
      Afix<-coef(fit.A)
      Bfix<-coef(fit.B)
      Cfix<-coef(fit.C)
      b.hat<-c(Afix,Bfix,Cfix)
      
      Lambda.hat[1,1]<-(summary(fit.A)$sigma)^2
      Lambda.hat[2,2]<-(summary(fit.B)$sigma)^2
      Lambda.hat[3,3]<-(summary(fit.C)$sigma)^2
    }else
    {
      fit.b<-lm(cbind(At,Bt,Ct)~1)
      aov.b<-Anova(fit.b)
      b.hat<-as.vector(coef(fit.b))
      Lambda.hat<-aov.b$SSPE/aov.b$error.df
    }
  }else
  {
    if(error.indep)
    {
      bt<-rbind(At,Bt,Ct)
      group<-c(rep("A",length(At)),rep("B",length(Bt)),rep("C",length(Ct)))
      fit<-lm(bt~group)
      b.hat<-as.vector(by(bt[,1],group,mean,na.rm=TRUE))
      diag(Lambda.hat)<-rep((summary(fit)$sigma)^2,3)
    }else
    {
      warning("This variance structure is not valid! The errors are assumed to be independent.")
      bt<-rbind(At,Bt,Ct)
      group<-c(rep("A",length(At)),rep("B",length(Bt)),rep("C",length(Ct)))
      fit<-lm(bt~group)
      b.hat<-as.vector(by(bt[,1],group,mean,na.rm=TRUE))
      diag(Lambda.hat)<-rep((summary(fit)$sigma)^2,3)
    }
  }
  ######################################
  
  ######################################
  # max.itr=0: the Two-stage approach
  # max.itr>0: h-likelihood approach
  diff<-100
  s<-0
  while(s<max.itr&diff>=tol)
  {
    bi.new<-matrix(NA,N,3)
    for(i in 1:N)
    {
      dd<-dat[[i]]
      Y<-cbind(dd$M,dd$R)
      X<-cbind(dd$Z,dd$M)
      
      if(Sigma.update)
      {
        Theta<-matrix(c(At[i],0,Ct[i],Bt[i]),2,2)
        S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
        
        sigma1[i]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/(nrow(dd)*(1-delta^2))
        sigma2[i]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/(nrow(dd)*(1-delta^2))
      }
      
      P<-matrix(c(-delta*sqrt(sigma2[i]/sigma1[i]),0,0,1,1,0),2,3)
      Q<-c(0,delta*sqrt(sigma2[i]/sigma1[i]))
      V<-matrix(c(1,0,0),nrow=1)
      # Update parameters
      bi.new[i,]<-solve(t(X%*%P)%*%(X%*%P)/(sigma2[i]*(1-delta^2))+t(dd$Z%*%V)%*%(dd$Z%*%V)/sigma1[i]+solve(Lambda.hat))%*%
        (t(X%*%P)%*%(dd$R-X%*%Q)/(sigma2[i]*(1-delta^2))+t(dd$Z%*%V)%*%dd$M/sigma1[i]+solve(Lambda.hat)%*%b.hat)
    }
    b.new<-apply(bi.new,2,mean)
    Lambda.new<-matrix(0,3,3)
    Lambda.tmp<-t(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))%*%(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))/N
    if(!error.var.equal)
    {
      if(error.indep)
      {
        diag(Lambda.new)<-diag(Lambda.tmp)
      }else
      {
        Lambda.new<-Lambda.tmp
      }
      
      # variance constraint
      if(sum(var.constraint==FALSE)==0)
      {
        # A
        if(Lambda.new[1,1]<Lambda.confint[1,1])
        {
          Lambda.new[1,1]<-Lambda.confint[1,1]
        }else
          if(Lambda.new[1,1]>Lambda.confint[1,2])
          {
            Lambda.new[1,1]<-Lambda.confint[1,2]
          }
        # B
        if(Lambda.new[2,2]<Lambda.confint[2,1])
        {
          Lambda.new[2,2]<-Lambda.confint[2,1]
        }else
          if(Lambda.new[2,2]>Lambda.confint[2,2])
          {
            Lambda.new[2,2]<-Lambda.confint[2,2]
          }
        # C
        if(Lambda.new[3,3]<Lambda.confint[3,1])
        {
          Lambda.new[3,3]<-Lambda.confint[3,1]
        }else
          if(Lambda.new[3,3]>Lambda.confint[3,2])
          {
            Lambda.new[3,3]<-Lambda.confint[3,2]
          }
      }
    }else
    {
      if(error.indep)
      {
        lambda2<-t(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))%*%c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))/(3*N)
        diag(Lambda.new)<-rep(lambda2,3)
      }else
      {
        warning("This variance structure is not valid! The errors are assumed to be independent.")
        lambda2<-t(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))%*%c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))/(3*N)
        diag(Lambda.new)<-rep(lambda2,3)
      }
    }
    
    
    diff<-max(abs(b.hat-b.new))
    
    b.hat<-b.new
    Lambda.hat<-Lambda.new
    At<-bi.new[,1]
    Bt<-bi.new[,2]
    Ct<-bi.new[,3]
    s<-s+1
    
    # print(c(s,diff))
  }
  ######################################
  
  ######################################
  # summary results
  HL<-cma.lm.h(dat,delta=delta,A.i=At,B.i=Bt,C.i=Ct,b=b.hat,Lambda=Lambda.hat,Sigma.update=Sigma.update)
  if(max.itr==0)
  {
    re.HL<-HL$h2
  }else
  {
    re.HL<-HL$h
  }
  
  AB.p<-b.hat[1]*b.hat[2]
  AB.d<-mean(C2t,na.rm=TRUE)-b.hat[3]
  coe.re<-matrix(NA,6,1)
  colnames(coe.re)<-c("Estimate")
  rownames(coe.re)<-c("A","C","B","C2","AB.prod","AB.diff")
  coe.re[,1]<-c(b.hat[1],b.hat[3],b.hat[2],mean(C2t,na.rm=TRUE),AB.p,AB.d)
  sigma.hat<-cbind(sigma1,sigma2)
  colnames(sigma.hat)<-c("E1","E2")
  ######################################
  
  re<-list(delta=delta,Coefficients=coe.re,Lambda=Lambda.hat,Sigma=sigma.hat,HL=re.HL,convergence=(s<max.itr|max.itr==0),
           Var.constraint=Lambda.confint)
  
  return(re)
}

### ---- from macc/macc/R/cma.h.R ----
cma.h <-
function(dat,delta=0,A.ik,B.ik,C.ik,b,u,Phi,Lambda,random.indep=TRUE,u.int=FALSE,Sigma.update=FALSE)
{
  N<-length(dat)
  K<-length(dat[[1]])
  
  sigma1=sigma2=n<-matrix(NA,N,K)
  h11=h12=h13=h14=h2=h3<-0
  for(i in 1:N)
  {
    for(j in 1:K)
    {
      dd<-dat[[i]][[j]]
      n[i,j]<-nrow(dd)
      if(Sigma.update)
      {
        Y<-cbind(dd$M,dd$R)
        X<-cbind(dd$Z,dd$M)
        Theta<-matrix(c(A.ik[i,j],0,C.ik[i,j],B.ik[i,j]),2,2)
        S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
        
        sigma1[i,j]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/(nrow(dd)*(1-delta^2))
        sigma2[i,j]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/(nrow(dd)*(1-delta^2))
      }else
      {
        re<-cma.uni.delta(dd,delta)
        sigma1[i,j]<-re$Sigma[1,1]
        sigma2[i,j]<-re$Sigma[2,2]
      }
      
      P<-matrix(c(-delta*sqrt(sigma2[i,j]/sigma1[i,j]),0,0,1,1,0),nrow=2,ncol=3)
      Q<-c(0,delta*sqrt(sigma2[i,j]/sigma1[i,j]))
      X<-cbind(dd$Z,dd$M)
      b.ik<-c(A.ik[i,j],B.ik[i,j],C.ik[i,j])
      V<-matrix(c(1,0,0),nrow=1)
      h11<-h11-(log(sigma2[i,j]*(1-delta^2))*nrow(dd)/2)
      h12<-h12-(t(dd$R-X%*%P%*%b.ik-X%*%Q)%*%(dd$R-X%*%P%*%b.ik-X%*%Q)/(2*sigma2[i,j]*(1-delta^2)))[1,1]
      h13<-h13-log(sigma1[i,j])*nrow(dd)/2
      h14<-h14-(t(dd$M-dd$Z%*%V%*%b.ik)%*%(dd$M-dd$Z%*%V%*%b.ik)/(2*sigma1[i,j]))[1,1]
      h2<-h2-log(det(Lambda))/2-(t(b.ik-b-u[i,])%*%solve(Lambda)%*%(b.ik-b-u[i,])/2)[1,1]
    }
    h3<-h3-log(det(Phi))/2-(t(u[i,])%*%solve(Phi)%*%u[i,]/2)[1,1]
  }
  const1<--log(2*pi)*sum(n)
  const2<--log(2*pi)*N*K*3/2
  const3<--log(2*pi)*N*3/2
  
  if(random.indep&u.int)
  {
    l.A=l.B=l.C<-0
    V.A<-diag(rep(Lambda[1,1],K,K))+Phi[1,1]*matrix(1,K,K)
    V.B<-diag(rep(Lambda[2,2],K,K))+Phi[2,2]*matrix(1,K,K)
    V.C<-diag(rep(Lambda[3,3],K,K))+Phi[3,3]*matrix(1,K,K)
    for(i in 1:N)
    {
      l.A<-l.A+(-log(det(V.A))/2-t(A.ik[i,]-b[1])%*%solve(V.A)%*%(A.ik[i,]-b[1]))[1,1]
      l.B<-l.B+(-log(det(V.B))/2-t(B.ik[i,]-b[2])%*%solve(V.B)%*%(B.ik[i,]-b[2]))[1,1]
      l.C<-l.C+(-log(det(V.C))/2-t(C.ik[i,]-b[3])%*%solve(V.C)%*%(C.ik[i,]-b[3]))[1,1]
    }
    h2<-l.A+l.B+l.C
    h<-h11+h12+h13+h14+h2
    re<-data.frame(h1=(const1+h11+h12+h13+h14),h2=const2+h2,h=const1+const2+h)
  }else
  {
    if(!random.indep&u.int)
    {
      warning("Independence in random effects assumption is enforced when calculating the likelihood.")
    }
    h<-h11+h12+h13+h14+h2+h3
    re<-data.frame(h1=(const1+h11+h12+h13+h14),h2=const2+h2,h3=const3+h3,h=const1+const2+const3+h)
  }
  return(re)
}

### ---- from macc/macc/R/cma.lm.h.R ----
cma.lm.h <-
function(dat,delta=0,A.i,B.i,C.i,b,Lambda,Sigma.update=FALSE)
{
  N<-length(dat)
  K<-1
  
  sigma1=sigma2=n<-matrix(NA,N,K)
  h11=h12=h2<-0
  for(i in 1:N)
  {
    dd<-dat[[i]]
    n[i,1]<-nrow(dd)
    
    Y<-cbind(dd$M,dd$R)
    X<-cbind(dd$Z,dd$M)
    if(Sigma.update)
    {
      Theta<-matrix(c(A.i[i],0,C.i[i],B.i[i]),2,2)
      S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
      
      sigma1[i]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/(n[i,1]*(1-delta^2))
      sigma2[i]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/(n[i,1]*(1-delta^2))
    }else
    {
      re<-cma.uni.delta(dd,delta)
      sigma1[i]<-re$Sigma[1,1]
      sigma2[i]<-re$Sigma[2,2]
    }
    
    P<-matrix(c(-delta*sqrt(sigma2[i]/sigma1[i]),0,0,1,1,0),2,3)
    Q<-c(0,delta*sqrt(sigma2[i]/sigma1[i]))
    b.i<-c(A.i[i],B.i[i],C.i[i])
    V<-matrix(c(1,0,0),nrow=1)
    
    h11<-h11-n[i,1]*log(sigma2[i]*(1-delta^2))/2-t(dd$R-X%*%(P%*%b.i+Q))%*%(dd$R-X%*%(P%*%b.i+Q))/(2*sigma2[i]*(1-delta^2))
    h12<-h12-n[i,1]*log(sigma1[i])/2-t(dd$M-dd$Z%*%V%*%b.i)%*%(dd$M-dd$Z%*%V%*%b.i)/(2*sigma1[i])
    h2<-h2-log(det(Lambda))/2-t(b.i-b)%*%solve(Lambda)%*%(b.i-b)/2
  }
  const1<--log(2*pi)*sum(n)
  const2<--log(2*pi)*N*3/2
  
  h<-(h11+h12+h2)[1,1]
  
  re<-data.frame(h1=(const1+h11+h12)[1,1],h2=(const2+h2)[1,1],h=(const1+const2+h))
  
  return(re)
}

### ---- from macc/macc/R/cma.lm.R ----
cma.lm <-
function(dat,model.type=c("single","multilevel"),method=c("HL","TS","HL-TS"),delta=NULL,sens.plot=FALSE,
                 interval=c(-0.9,0.9),tol=1e-4,max.itr=500,conf.level=0.95,
                 error.indep=FALSE,error.var.equal=FALSE,Sigma.update=FALSE,var.constraint=FALSE,
                 legend.pos="topright",xlab=expression(delta),ylab=expression(hat(AB)),
                 cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
  
{
  if(model.type[1]=="single")
  {
    if(is.null(delta))
    {
      delta<-0
    }
    re<-cma.uni.delta(dat,delta=delta,conf.level=conf.level)
    if(sens.plot)
    {
      re.sens<-cma.uni.sens(dat,conf.level=conf.level)
      cma.uni.plot(re.sens,re,delta=NULL,legend.pos=legend.pos,xlab=xlab,ylab=ylab,
                   cex.lab=cex.lab,cex.axis=cex.axis,lgd.cex=lgd.cex,lgd.pt.cex=1,plot.delta0=plot.delta0,...)
    }
  }else
    if(model.type[1]=="multilevel")
    {
      if(is.null(delta))
      {
        if(method[1]=="TS")
        {
          system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=0,tol=tol,
                                    error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                    var.constraint=var.constraint,maximum=TRUE))
          re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,error.var.equal=error.var.equal,
                           Sigma.update=Sigma.update,var.constraint=var.constraint)
        }else
        {
          system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=max.itr,tol=tol,
                                    error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                    var.constraint=var.constraint,maximum=TRUE))
          if(method[1]=="HL")
          {
            re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=max.itr,tol=tol,error.indep=error.indep,error.var.equal=error.var.equal,
                             Sigma.update=Sigma.update,var.constraint=var.constraint)
          }
          if(method[1]=="HL-TS")
          {
            re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,error.var.equal=error.var.equal,
                             Sigma.update=Sigma.update,var.constraint=var.constraint)
          }
        }
      }else
      {
        if(method[1]=="TS")
        {
          re<-cma.delta.lm(dat,delta=delta,max.itr=0,tol=tol,error.indep=error.indep,error.var.equal=error.var.equal,
                           Sigma.update=Sigma.update,var.constraint=var.constraint)
        }
        if(method[1]=="HL")
        {
          re<-cma.delta.lm(dat,delta=delta,max.itr=max.itr,tol=tol,error.indep=error.indep,error.var.equal=error.var.equal,
                           Sigma.update=Sigma.update,var.constraint=var.constraint)
        }
      }
    }
  return(re)
}

### ---- from macc/macc/R/cma.lm.ts.R ----
cma.lm.ts <-
function(dat,Sigma.update=FALSE)
{
  N<-length(dat)
  K<-1
  
  A.hat<-rep(NA,N)
  C2.hat=C.tilde=B.tilde=beta=gamma<-A.hat
  
  i<-1
  for(i in 1:N)
  {
    dd<-dat[[i]]
    A.hat[i]<-(t(dd$Z)%*%dd$M)/(t(dd$Z)%*%dd$Z)
    C2.hat[i]<-(t(dd$Z)%*%dd$R)/(t(dd$Z)%*%dd$Z)
    Sigma.B<-t(cbind(dd$M-dd$Z*A.hat[i],dd$R-dd$Z*C2.hat[i]))%*%cbind(dd$M-dd$Z*A.hat[i],dd$R-dd$Z*C2.hat[i])/nrow(dd)
    
    C.tilde[i]<-(t(dd$M)%*%dd$M%*%t(dd$Z)%*%dd$R-t(dd$Z)%*%dd$M%*%t(dd$M)%*%dd$R)/
      (t(dd$Z)%*%dd$Z%*%t(dd$M)%*%dd$M-t(dd$M)%*%dd$Z%*%t(dd$Z)%*%dd$M)
    gamma[i]<-(sqrt(det(Sigma.B))/Sigma.B[1,1])*((t(dd$Z)%*%dd$M)/(t(dd$Z)%*%dd$Z))
    
    B.tilde[i]<-(t(dd$Z)%*%dd$Z%*%t(dd$M)%*%dd$R-t(dd$M)%*%dd$Z%*%t(dd$Z)%*%dd$R)/
      (t(dd$Z)%*%dd$Z%*%t(dd$M)%*%dd$M-t(dd$M)%*%dd$Z%*%t(dd$Z)%*%dd$M)
    beta[i]<-sqrt(det(Sigma.B))/Sigma.B[1,1]
  }
  
  tau<-(sum((beta-mean(beta))*(B.tilde-mean(B.tilde)))-sum((gamma-mean(gamma))*(C.tilde-mean(C.tilde))))/
    ((sum(beta^2)-N*(mean(beta))^2)+(sum(gamma^2)-N*(mean(gamma))^2))
  
  delta.est<-tau/(sqrt(1+tau^2))
  
  Bt<-B.tilde-tau*beta
  Ct<-C.tilde-tau*gamma
  A.est<-mean(A.hat)
  B.est<-mean(B.tilde)-tau*mean(beta)
  C.est<-mean(C.tilde)+tau*mean(gamma)
  b.hat<-c(A.est,B.est,C.est)
  
  Lambda.est<-matrix(0,3,3)
  diag(Lambda.est)<-(t(A.hat-A.est)%*%(A.hat-A.est)+t(Bt-B.est)%*%(Bt-B.est)+t(Ct-C.est)%*%(Ct-C.est))/(3*N)
  
  sigma1=sigma2<-rep(NA,N)
  for(i in 1:N)
  {
    dd<-dat[[i]]
    
    Y<-cbind(dd$M,dd$R)
    X<-cbind(dd$Z,dd$M)
    if(Sigma.update)
    {
      Theta<-matrix(c(A.hat[i],0,Ct[i],Bt[i]),2,2)
      S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
      
      sigma1[i]<-(S[1,1]-delta.est*S[1,2]*sqrt(S[1,1]/S[2,2]))/(nrow(dd)*(1-delta.est^2))
      sigma2[i]<-(S[2,2]-delta.est*S[1,2]*sqrt(S[2,2]/S[1,1]))/(nrow(dd)*(1-delta.est^2))
    }else
    {
      re<-cma.uni.delta(dd,delta.est)
      sigma1[i]<-re$Sigma[1,1]
      sigma2[i]<-re$Sigma[2,2]
    }
  }
  ########################################################################
  HL<-cma.lm.h(dat,delta=delta.est,A.i=A.hat,B.i=Bt,C.i=Ct,b=b.hat,Lambda=Lambda.est,Sigma.update=Sigma.update)
  
  AB.p<-b.hat[1]*b.hat[2]
  AB.d<-mean(C2.hat,na.rm=TRUE)-b.hat[3]
  coe.re<-matrix(NA,6,1)
  colnames(coe.re)<-c("Estimate")
  rownames(coe.re)<-c("A","C","B","C2","AB.prod","AB.diff")
  coe.re[,1]<-c(b.hat[1],b.hat[3],b.hat[2],mean(C2.hat,na.rm=TRUE),AB.p,AB.d)
  sigma.hat<-cbind(sigma1,sigma2)
  colnames(sigma.hat)<-c("E1","E2")
  
  re<-list(delta=delta.est,Coefficients=coe.re,Lambda=Lambda.est,Sigma=sigma.hat,HL=HL$h)
  
  return(re)
}

### ---- from macc/macc/R/cma.mix.R ----
cma.mix <-
function(dat,model.type=c("single","multilevel"),method=c("HL","TS","HL-TS"),delta=NULL,sens.plot=FALSE,
                  interval=c(-0.90,0.90),tol=10e-4,max.itr=50,conf.level=0.95,optimizer=c("bobyqa","Nelder_Mead","optimx"),
                  mix.pkg=c("lme4","nlme"),random.indep=TRUE,random.var.equal=TRUE,u.int=FALSE,Sigma.update=FALSE,
                  var.constraint=FALSE,random.var.update=TRUE,logLik.type=c("logLik","HL"),
                  legend.pos="topright",xlab=expression(delta),ylab=expression(hat(AB)),
                  cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
{
  if(model.type[1]=="single")
  {
    if(is.null(delta))
    {
      delta<-0
    }
    re<-cma.uni.delta(dat,delta=delta,conf.level=conf.level)
    if(sens.plot)
    {
      re.sens<-cma.uni.sens(dat,conf.level=conf.level)
      cma.uni.plot(re.sens,re,delta=NULL,legend.pos=legend.pos,xlab=xlab,ylab=ylab,
                   cex.lab=cex.lab,cex.axis=cex.axis,lgd.cex=lgd.cex,lgd.pt.cex=1,plot.delta0=plot.delta0,...)
    }
  }else
    if(model.type[1]=="multilevel")
    {
      if(is.null(delta))
      {
        if(method[1]=="TS")
        {
          system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=0,optimizer=optimizer,
                                    mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                    u.int=u.int,Sigma.update=Sigma.update,logLik.type=logLik.type,maximum=TRUE))
          re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                          random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int)
        }else
        {
          system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=max.itr,optimizer=optimizer,
                                    mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                    u.int=u.int,Sigma.update=Sigma.update,var.constraint=var.constraint,
                                    random.var.update=random.var.update,logLik.type=logLik.type,maximum=TRUE))
          if(method[1]=="HL-TS")
          {
            re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                            random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int)
          }
          if(method[1]=="HL")
          {
            re<-cma.uni.mix.hl(dat,delta=re1$maximum,tol=tol,max.itr=max.itr,alpha=1-conf.level,random.indep=random.indep,
                               optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                               Sigma.update=Sigma.update,var.constraint=var.constraint,random.var.update=random.var.update)
          }
        }
      }else
      {
        if(method[1]=="TS")
        {
          re<-cma.uni.mix(dat,delta=delta,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                          random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int)
        }
        if(method[1]=="HL")
        {
          re<-cma.uni.mix.hl(dat,delta=delta,tol=tol,max.itr=max.itr,alpha=1-conf.level,random.indep=random.indep,
                             optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                             Sigma.update=Sigma.update,var.constraint=var.constraint,random.var.update=random.var.update)
        }
      }
    }
  return(re)
  if(as.numeric(re$Var.comp[1])<1e-5)
  {
    warning("The variance of A's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[2])<1e-5)
  {
    warning("The variance of C's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[3])<1e-5)
  {
    warning("The variance of B's random effect is less than 1e-5.")
  }
}

### ---- from macc/macc/R/cma.uni.delta.R ----
cma.uni.delta <-
function(dat,delta=0,conf.level=0.95)
{
  Z<-matrix(dat$Z,ncol=1)
  M<-matrix(dat$M,ncol=1)
  R<-matrix(dat$R,ncol=1)
  
  n<-nrow(Z)
  
  z.alpha<-qnorm(1-(1-conf.level)/2)
  
  ########################################
  # total effect model
  X<-cbind(Z,M)
  fit.M<-lm(M~0+Z)
  fit.R2<-lm(R~0+Z)
  
  beta.M<-coef(fit.M)
  beta.R2<-coef(fit.R2)
  
  Sigma.B.hat<-(t(cbind(M-Z%*%beta.M,R-Z%*%beta.R2))%*%cbind(M-Z%*%beta.M,R-Z%*%beta.R2))/n
  sigma12.hat<-Sigma.B.hat[1,1]
  sigma22.hat<-det(Sigma.B.hat)/(Sigma.B.hat[1,1]*(1-delta^2))
  Sigma.hat<-matrix(c(sigma12.hat,delta*sqrt(sigma12.hat*sigma22.hat),delta*sqrt(sigma12.hat*sigma22.hat),sigma22.hat),2,2)
  
  # variance of beta.M and beta.R2
  beta.M.var<-Sigma.B.hat[1,1]*solve(t(Z)%*%Z)
  beta.R2.var<-Sigma.B.hat[2,2]*solve(t(X)%*%X)
  #########################################
  
  #########################################
  # coefficient estimate given delta
  
  # C'
  C2.hat<-coef(fit.R2)[1]
  C2.hat.se<-sqrt(beta.R2.var[1,1])
  
  # A, B, and C
  zz<-(t(Z)%*%Z)[1,1]
  zm<-(t(Z)%*%M)[1,1]
  mz<-(t(M)%*%Z)[1,1]
  zr<-(t(Z)%*%R)[1,1]
  mm<-(t(M)%*%M)[1,1]
  mr<-(t(M)%*%R)[1,1]
  
  A.hat<-zm/zz
  C.hat<-(mm*zr-zm*mr)/(zz*mm-mz*zm)+(delta*sqrt(sigma22.hat/sigma12.hat))*(zm/zz)
  B.hat<-(zz*mr-mz*zr)/(zz*mm-mz*zm)-(delta*sqrt(sigma22.hat/sigma12.hat))
  # inverse of Fisher infromation
  coef.cov<-matrix(0,3,3)
  coef.cov[1,1]<-n*sigma12.hat/zz
  coef.cov[1,2]=coef.cov[2,1]<-n*delta*sqrt(sigma12.hat*sigma22.hat)/zz
  coef.cov[2,2]<-sigma22.hat*(A.hat^2*zz+n*sigma12.hat-delta^2*A.hat^2*zz)/(sigma12.hat*zz)
  coef.cov[2,3]=coef.cov[3,2]<--sigma22.hat*(1-delta^2)*A.hat/sigma12.hat
  coef.cov[3,3]<-sigma22.hat*(1-delta^2)/sigma12.hat
  
  A.hat.se<-sqrt(coef.cov[1,1]/n)
  C.hat.se<-sqrt(coef.cov[2,2]/n)
  B.hat.se<-sqrt(coef.cov[3,3]/n)
  
  ABp.hat<-A.hat*B.hat
  ABp.hat.se<-sqrt(A.hat^2*B.hat.se^2+B.hat^2*A.hat.se^2)
  ABd.hat<-C2.hat-C.hat
  ABd.hat.se<-sqrt(C2.hat.se^2+C.hat.se^2)
  
  cma.re<-matrix(NA,6,4)
  rownames(cma.re)<-c("A","C","B","C2","ABp","ABd")
  colnames(cma.re)<-c("Estimate","SE","LB","UB")
  cma.re[,1]<-c(A.hat,C.hat,B.hat,C2.hat,ABp.hat,ABd.hat)
  cma.re[,2]<-c(A.hat.se,C.hat.se,B.hat.se,C2.hat.se,ABp.hat.se,ABd.hat.se)
  cma.re[,3]<-cma.re[,1]-z.alpha*cma.re[,2]
  cma.re[,4]<-cma.re[,1]+z.alpha*cma.re[,2]
  
  D.hat<-matrix(c(A.hat,0,C.hat,B.hat),2,2)
  
  re<-list(Coefficients=cma.re,D=D.hat,Sigma=Sigma.hat,delta=delta)
  
  return(re)
}

### ---- from macc/macc/R/cma.uni.mix.dhl.R ----
cma.uni.mix.dhl <-
function(dat,delta,tol=10e-4,max.itr=50,alpha=0.05,optimizer=c("bobyqa","Nelder_Mead","optimx"),
                          mix.pkg=c("lme4","nlme"),random.indep=TRUE,random.var.equal=TRUE,u.int=FALSE,
                          Sigma.update=FALSE,var.constraint=FALSE,random.var.update=TRUE,logLik.type=c("logLik","HL"))
{
  if(max.itr==0)
  {
    if(logLik.type[1]=="HL")
    {
      return(cma.uni.mix(dat,delta,conf.level=1-alpha,optimizer=optimizer,mix.pkg=mix.pkg,random.indep=random.indep,
                         random.var.equal=random.var.equal,u.int=u.int)$HL)
    }else
    {
      return(cma.uni.mix(dat,delta,conf.level=1-alpha,optimizer=optimizer,mix.pkg=mix.pkg,random.indep=random.indep,
                         random.var.equal=random.var.equal,u.int=u.int)$logLik)
    }
  }else
  {
    return(cma.uni.mix.hl(dat,delta=delta,tol=tol,max.itr=max.itr,alpha=alpha,random.indep=random.indep,
                          optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                          Sigma.update=Sigma.update,var.constraint=var.constraint,random.var.update=random.var.update)$HL)
  }
}

### ---- from macc/macc/R/cma.uni.mix.hl.R ----
cma.uni.mix.hl <-
function(dat,delta,tol=1e-4,max.itr=50,alpha=0.05,random.indep=TRUE,optimizer=c("bobyqa","Nelder_Mead","optimx"),
                         mix.pkg=c("lme4","nlme"),random.var.equal=TRUE,u.int=FALSE,Sigma.update=FALSE,var.constraint=FALSE,
                         random.var.update=TRUE)
{
  N<-length(dat)
  K<-length(dat[[1]])
  n<-matrix(NA,N,K)
  
  sigma1=sigma2<-matrix(NA,N,K)
  At=Ct=Bt=C2t<-matrix(NA,N,K)
  for(i in 1:N)
  {
    for(k in 1:K)
    {
      dd<-dat[[i]][[k]]
      n[i,k]<-nrow(dd)
      re2<-cma.uni.delta(dd,delta=delta)
      
      sigma1[i,k]<-re2$Sigma[1,1]
      sigma2[i,k]<-re2$Sigma[2,2]
      
      At[i,k]<-re2$Coefficients[1,1]
      Bt[i,k]<-re2$Coefficients[3,1]
      Ct[i,k]<-re2$Coefficients[2,1]
      C2t[i,k]<-re2$Coefficients[4,1]
    }
  }
  coe<-rep(c("A","C","B"),each=N*K)
  Sub<-factor(rep(1:N,K))
  Sub.all<-rep(Sub,3)
  coet<-c(c(At),c(Ct),c(Bt))
  
  if(random.var.equal==FALSE)
  {
    if(random.indep==TRUE)
    {
      vec.At<-c(At)
      vec.Bt<-c(Bt)
      vec.Ct<-c(Ct)
      if(mix.pkg[1]=="nlme")
      {
        if(optimizer[1]=="optimx")
        {
          fit.A<-lme(vec.At~1,random=~1|Sub,control=lmeControl(opt="optim"))
          fit.B<-lme(vec.Bt~1,random=~1|Sub,control=lmeControl(opt="optim"))
          fit.C<-lme(vec.Ct~1,random=~1|Sub,control=lmeControl(opt="optim"))
        }else
        {
          fit.A<-lme(vec.At~1,random=~1|Sub)
          fit.B<-lme(vec.Bt~1,random=~1|Sub)
          fit.C<-lme(vec.Ct~1,random=~1|Sub)
        }
        Afix<-as.numeric(summary(fit.A)$coefficients$fixed[1])
        Bfix<-as.numeric(summary(fit.B)$coefficients$fixed[1])
        Cfix<-as.numeric(summary(fit.C)$coefficients$fixed[1])
        Phi<-diag(c(as.numeric(VarCorr(fit.A)[1,1]),as.numeric(VarCorr(fit.B)[1,1]),as.numeric(VarCorr(fit.C)[1,1])))
        Lambda<-diag(c(as.numeric(VarCorr(fit.A)[2,1]),as.numeric(VarCorr(fit.B)[2,1]),as.numeric(VarCorr(fit.C)[2,1])))
        u<-cbind(as.matrix(ranef(fit.A)),as.matrix(ranef(fit.B)),as.matrix(ranef(fit.C)))
        colnames(u)<-c("A","B","C")
        b0<-c(Afix,Bfix,Cfix)
        cor.AB=cor.AC=cor.BC<-0
        LL<-as.numeric(logLik(fit.A,REML=FALSE)+logLik(fit.B,REML=FALSE)+logLik(fit.C,REML=FALSE))
        if(is.matrix(var.constraint))
        {
          if(nrow(var.constraint)==6)
          {
            Phi.confint<-var.constraint[1:3,]
            Lambda.confint<-var.constraint[4:6,]
            colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
            rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
          }else
          {
            warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            intA<-intervals(fit.A)
            intB<-intervals(fit.B)
            intC<-intervals(fit.C)
            Phi.confint[1,]<-as.matrix(intA[[2]]$Sub[1,c(1,3)])^2
            Phi.confint[2,]<-as.matrix(intB[[2]]$Sub[1,c(1,3)])^2
            Phi.confint[3,]<-as.matrix(intC[[2]]$Sub[1,c(1,3)])^2
            Lambda.confint[1,]<-as.matrix(intA[[3]][c(1,3)])^2
            Lambda.confint[2,]<-as.matrix(intB[[3]][c(1,3)])^2
            Lambda.confint[3,]<-as.matrix(intC[[3]][c(1,3)])^2 
          }
        }else
          if(var.constraint)
          {
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            intA<-intervals(fit.A)
            intB<-intervals(fit.B)
            intC<-intervals(fit.C)
            Phi.confint[1,]<-as.matrix(intA[[2]]$Sub[1,c(1,3)])^2
            Phi.confint[2,]<-as.matrix(intB[[2]]$Sub[1,c(1,3)])^2
            Phi.confint[3,]<-as.matrix(intC[[2]]$Sub[1,c(1,3)])^2
            Lambda.confint[1,]<-as.matrix(intA[[3]][c(1,3)])^2
            Lambda.confint[2,]<-as.matrix(intB[[3]][c(1,3)])^2
            Lambda.confint[3,]<-as.matrix(intC[[3]][c(1,3)])^2
          }
      }else
      {
        if(optimizer[1]=="optimx")
        {
          fit.A<-lmer(vec.At~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          fit.B<-lmer(vec.Bt~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          fit.C<-lmer(vec.Ct~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        }else
        {
          fit.A<-lmer(vec.At~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
          fit.B<-lmer(vec.Bt~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
          fit.C<-lmer(vec.Ct~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
        }
        Afix<-as.numeric(summary(fit.A)$coefficients[1])
        Bfix<-as.numeric(summary(fit.B)$coefficients[1])
        Cfix<-as.numeric(summary(fit.C)$coefficients[1])
        Phi<-diag(c(as.numeric((attr(VarCorr(fit.A)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit.B)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit.C)[[1]],"stddev")^2)[1])))
        Lambda<-diag(c(as.numeric((attr(VarCorr(fit.A),"sc")^2)[1]),
                       as.numeric((attr(VarCorr(fit.B),"sc")^2)[1]),
                       as.numeric((attr(VarCorr(fit.C),"sc")^2)[1])))
        u<-cbind(as.matrix(ranef(fit.A)$Sub),as.matrix(ranef(fit.B)$Sub),as.matrix(ranef(fit.C)$Sub))
        colnames(u)<-c("A","B","C")
        b0<-c(Afix,Bfix,Cfix)
        cor.AB=cor.AC=cor.BC<-0
        LL<-as.numeric(logLik(fit.A,REML=FALSE)+logLik(fit.B,REML=FALSE)+logLik(fit.C,REML=FALSE))
        if(is.matrix(var.constraint))
        {
          if(nrow(var.constraint)==6)
          {
            Phi.confint<-var.constraint[1:3,]
            Lambda.confint<-var.constraint[4:6,]
            colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
            rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
          }else
          {
            warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            intA<-confint(fit.A)
            intB<-confint(fit.B)
            intC<-confint(fit.C)
            Phi.confint[1,]<-as.matrix(intA[1,])^2
            Phi.confint[2,]<-as.matrix(intB[1,])^2
            Phi.confint[3,]<-as.matrix(intC[1,])^2
            Lambda.confint[1,]<-as.matrix(intA[2,])^2
            Lambda.confint[2,]<-as.matrix(intB[2,])^2
            Lambda.confint[3,]<-as.matrix(intC[2,])^2 
          }
        }else
          if(var.constraint)
          {
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            intA<-confint(fit.A)
            intB<-confint(fit.B)
            intC<-confint(fit.C)
            Phi.confint[1,]<-as.matrix(intA[1,])^2
            Phi.confint[2,]<-as.matrix(intB[1,])^2
            Phi.confint[3,]<-as.matrix(intC[1,])^2
            Lambda.confint[1,]<-as.matrix(intA[2,])^2
            Lambda.confint[2,]<-as.matrix(intB[2,])^2
            Lambda.confint[3,]<-as.matrix(intC[2,])^2 
          }
      }
    }else
    {
      if(mix.pkg[1]=="nlme")
      {
        if(is.matrix(var.constraint))
        {
          if(nrow(var.constraint)==6)
          {
            Phi.confint<-var.constraint[1:3,]
            Lambda.confint<-var.constraint[4:6,]
            colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
            rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
          }else
          {
            warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
          }
        }else
          if(var.constraint)
          {
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint 
          }
        if(optimizer[1]=="optimx")
        {
          fit0<-lme(c(coet)~0+coe,random=~0+coe|Sub.all,control=lmeControl(opt="optim"))
        }else
        {
          fit0<-lme(c(coet)~0+coe,random=~0+coe|Sub.all)
        }
        fit<-NULL
        try(fit<-update(fit0,weights=varIdent(form=~1|factor(coe))))
        if(is.null(fit)==TRUE)
        {
          warning("Equal-varaince assumption in random error is applied instead.")
          fit<-fit0
          Lambda<-diag(rep(summary(fit)$sigma^2,3))
          if(!is.matrix(var.constraint))
          {
            if(var.constraint)
            {
              Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-(intervals(fit)[[3]][c(1,3)])^2 
            }
          }
        }else
        {
          lambda.tmp<-summary(fit)$sigma*coef(fit$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE)
          Lambda<-diag(c((lambda.tmp[1])^2,(lambda.tmp[3])^2,(lambda.tmp[2])^2))
          if(is.matrix(var.constraint)==FALSE)
          {
            if(var.constraint==TRUE)
            {
              int.fit<-intervals(fit)
              Lambda.confint[1,]<-(int.fit[[4]][c(1,3)])^2
              Lambda.confint[2,]<-(int.fit[[3]][2,c(1,3)])^2*Lambda[1,1]
              Lambda.confint[3,]<-(int.fit[[3]][1,c(1,3)])^2*Lambda[1,1]
            }
          }
        }
        Afix<-as.numeric(summary(fit)$coefficients$fixed[1])
        Bfix<-as.numeric(summary(fit)$coefficients$fixed[2])
        Cfix<-as.numeric(summary(fit)$coefficients$fixed[3])
        b0<-c(Afix,Bfix,Cfix)
        u<-as.matrix(ranef(fit))
        Phi<-diag(c(as.numeric(VarCorr(fit)[1,1]),as.numeric(VarCorr(fit)[2,1]),as.numeric(VarCorr(fit)[3,1])))
        cor.AB<-as.numeric(VarCorr(fit)[2,3])
        cor.AC<-as.numeric(VarCorr(fit)[3,3])
        cor.BC<-as.numeric(VarCorr(fit)[3,4])
        Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
        Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
        Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
        LL<-as.numeric(logLik(fit,REML=FALSE))
        if(!is.matrix(var.constraint))
        {
          if(var.constraint)
          {
            int.fit<-intervals(fit)
            Phi.confint[1,]<-as.matrix(int.fit[[2]]$Sub.all[1,c(1,3)])^2
            Phi.confint[2,]<-as.matrix(int.fit[[2]]$Sub.all[2,c(1,3)])^2
            Phi.confint[3,]<-as.matrix(int.fit[[2]]$Sub.all[3,c(1,3)])^2 
          } 
        }
      }else
      {
        warning("Equal-varaince assumption in random error is applied instead.")
        if(optimizer[1]=="optimx")
        {
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        }else
        {
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer=optimizer[1]))
        }
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
        b0<-c(summary(fit)$coefficients[1,1],summary(fit)$coefficients[2,1],summary(fit)$coefficients[3,1])
        Phi<-diag(c(as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[2]),
                    as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[3])))
        cor.AB<-attr(VarCorr(fit)[[1]],"correlation")[1,2]
        cor.AC<-attr(VarCorr(fit)[[1]],"correlation")[1,3]
        cor.BC<-attr(VarCorr(fit)[[1]],"correlation")[2,3]
        Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
        Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
        Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
        Lambda<-diag(rep(attr(VarCorr(fit),"sc")^2,3))
        LL<-as.numeric(logLik(fit,REML=FALSE))
        if(is.matrix(var.constraint))
        {
          if(nrow(var.constraint)==6)
          {
            Phi.confint<-var.constraint[1:3,]
            Lambda.confint<-var.constraint[4:6,]
            colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
            rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
          }else
          {
            warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            re.confint<-confint(fit)
            Phi.confint[1,]<-as.matrix(re.confint[1,])^2
            Phi.confint[2,]<-as.matrix(re.confint[4,])^2
            Phi.confint[3,]<-as.matrix(re.confint[6,])^2
            Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-as.matrix(re.confint[7,])^2 
          }
        }else
          if(var.constraint)
          {
            Lambda.confint<-matrix(NA,3,2)
            colnames(Lambda.confint)<-c("LB","UB")
            rownames(Lambda.confint)<-c("A","B","C")
            Phi.confint<-Lambda.confint
            re.confint<-confint(fit)
            Phi.confint[1,]<-as.matrix(re.confint[1,])^2
            Phi.confint[2,]<-as.matrix(re.confint[4,])^2
            Phi.confint[3,]<-as.matrix(re.confint[6,])^2
            Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-as.matrix(re.confint[7,])^2 
          }
      }
    }
  }else
  {
    if(mix.pkg[1]=="nlme")
    {
      if(optimizer[1]=="optimx")
      {
        fit<-lme(c(coet)~0+coe,random=~0+coe|Sub.all,control=lmeControl(opt="optim"))        
      }else
      {
        fit<-lme(c(coet)~0+coe,random=~0+coe|Sub.all)
      }
      Afix<-as.numeric(summary(fit)$coefficients$fixed[1])
      Bfix<-as.numeric(summary(fit)$coefficients$fixed[2])
      Cfix<-as.numeric(summary(fit)$coefficients$fixed[3])
      b0<-c(Afix,Bfix,Cfix)
      u<-as.matrix(ranef(fit))
      Lambda<-diag(rep(summary(fit)$sigma^2,3))
      Phi<-diag(c(as.numeric(VarCorr(fit)[1,1]),as.numeric(VarCorr(fit)[2,1]),as.numeric(VarCorr(fit)[3,1])))
      if(random.indep)
      {
        cor.AB=cor.AC=cor.BC<-0
      }else
      {
        cor.AB<-as.numeric(VarCorr(fit)[2,3])
        cor.AC<-as.numeric(VarCorr(fit)[3,3])
        cor.BC<-as.numeric(VarCorr(fit)[3,4])
      }
      Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
      Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
      Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
      LL<-as.numeric(logLik(fit,REML=FALSE))
      if(is.matrix(var.constraint))
      {
        if(nrow(var.constraint)>=4)
        {
          Phi.confint<-var.constraint[1:3,]
          Lambda.confint<-matrix(rep(var.constraint[4,],3),nrow=3,ncol=2,byrow=TRUE)
          colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
          rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
        }else
        {
          warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
          Lambda.confint<-matrix(NA,3,2)
          colnames(Lambda.confint)<-c("LB","UB")
          rownames(Lambda.confint)<-c("A","B","C")
          Phi.confint<-Lambda.confint
          int.fit<-intervals(fit)
          Phi.confint[1,]<-as.matrix(int.fit[[2]]$Sub.all[1,c(1,3)])^2
          Phi.confint[2,]<-as.matrix(int.fit[[2]]$Sub.all[2,c(1,3)])^2
          Phi.confint[3,]<-as.matrix(int.fit[[2]]$Sub.all[3,c(1,3)])^2
          Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-as.matrix(int.fit[[3]][c(1,3)])^2
        }
      }else
        if(var.constraint)
        {
          Lambda.confint<-matrix(NA,3,2)
          colnames(Lambda.confint)<-c("LB","UB")
          rownames(Lambda.confint)<-c("A","B","C")
          Phi.confint<-Lambda.confint
          int.fit<-intervals(fit)
          Phi.confint[1,]<-as.matrix(int.fit[[2]]$Sub.all[1,c(1,3)])^2
          Phi.confint[2,]<-as.matrix(int.fit[[2]]$Sub.all[2,c(1,3)])^2
          Phi.confint[3,]<-as.matrix(int.fit[[2]]$Sub.all[3,c(1,3)])^2
          Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-as.matrix(intervals(fit)[[3]][c(1,3)])^2
        }
    }else
    {
      if(optimizer[1]=="optimx")
      {
        fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
      }else
      {
        fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control=lmerControl(optimizer=optimizer[1]))
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
      }
      b0<-c(summary(fit)$coefficients[1,1],summary(fit)$coefficients[2,1],summary(fit)$coefficients[3,1])
      Phi<-diag(c(as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[1]),
                  as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[2]),
                  as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[3])))
      if(random.indep==TRUE)
      {
        cor.AB=cor.AC=cor.BC<-0
      }else
      {
        cor.AB<-attr(VarCorr(fit)[[1]],"correlation")[1,2]
        cor.AC<-attr(VarCorr(fit)[[1]],"correlation")[1,3]
        cor.BC<-attr(VarCorr(fit)[[1]],"correlation")[2,3]
      }
      Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
      Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
      Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
      Lambda<-diag(rep(attr(VarCorr(fit),"sc")^2,3))
      LL<-as.numeric(logLik(fit,REML=FALSE))
      if(is.matrix(var.constraint))
      {
        if(nrow(var.constraint)>=4)
        {
          Phi.confint<-var.constraint[1:3,]
          Lambda.confint<-matrix(rep(var.constraint[4,],3),nrow=3,ncol=2,byrow=TRUE)
          colnames(Phi.confint)=colnames(Lambda.confint)<-c("LB","UB")
          rownames(Phi.confint)=rownames(Lambda.confint)<-c("A","B","C")
        }else
        {
          warning("The number of intervals is not correct. The constraint intervals will be estimated instead.")
          Lambda.confint<-matrix(NA,3,2)
          colnames(Lambda.confint)<-c("LB","UB")
          rownames(Lambda.confint)<-c("A","B","C")
          Phi.confint<-Lambda.confint
          re.confint<-confint(fit)
          Phi.confint[1,]<-(re.confint[1,])^2
          Phi.confint[2,]<-(re.confint[4,])^2
          Phi.confint[3,]<-(re.confint[6,])^2
          Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-(re.confint[7,])^2
        }
      }else
        if(var.constraint==TRUE)
        {
          Lambda.confint<-matrix(NA,3,2)
          colnames(Lambda.confint)<-c("LB","UB")
          rownames(Lambda.confint)<-c("A","B","C")
          Phi.confint<-Lambda.confint
          re.confint<-confint(fit)
          Phi.confint[1,]<-(re.confint[1,])^2
          Phi.confint[2,]<-(re.confint[4,])^2
          Phi.confint[3,]<-(re.confint[6,])^2
          Lambda.confint[1,]=Lambda.confint[2,]=Lambda.confint[3,]<-(re.confint[7,])^2
        } 
    }
  }
  
  diff<-100
  s<-0
  while(s<max.itr&diff>=tol)
  {
    if(random.indep&!random.var.equal&u.int)
    {
      V.alpha<-matrix(Phi[1,1],K,K)+diag(rep(Lambda[1,1],K))
      V.beta<-matrix(Phi[2,2],K,K)+diag(rep(Lambda[2,2],K))
      V.gamma<-matrix(Phi[3,3],K,K)+diag(rep(Lambda[3,3],K))
      for(i in 1:N)
      {
        # A
        z1z=z1m=z1z.d=z1m.d<-NULL  
        for(k in 1:K)
        {
          dd<-dat[[i]][[k]]
          inv.O1<-diag(rep(1/sigma1[i,k],nrow(dd)))
          inv.O2<-diag(rep(1/sigma2[i,k],nrow(dd)))
          O12<-diag(rep(sqrt(sigma2[i,k]/sigma1[i,k]),nrow(dd)))
          
          z1z<-c(z1z,(t(dd$Z)%*%inv.O1%*%(dd$Z))[1,1])
          z1m<-c(z1m,(t(dd$Z)%*%inv.O1%*%(dd$M))[1,1])
          z1z.d<-c(z1z.d,(delta^2*t(dd$Z)%*%O12%*%inv.O2%*%O12%*%(dd$Z)/(1-delta^2))[1,1])
          z1m.d<-c(z1m.d,(delta*t(dd$Z)%*%O12%*%inv.O2%*%(dd$R-dd$Z*Ct[i,k]-dd$M*Bt[i,k]-delta*O12%*%dd$M)/(1-delta^2))[1,1])
        }
        At[i,]<-solve(diag(z1z+z1z.d)+solve(V.alpha))%*%(z1m-z1m.d+solve(V.alpha)%*%rep(b0[1],K))
        
        # B
        m2m=m2r.d<-NULL
        for(k in 1:K)
        {
          dd<-dat[[i]][[k]]
          inv.O1<-diag(rep(1/sigma1[i,k],nrow(dd)))
          inv.O2<-diag(rep(1/sigma2[i,k],nrow(dd)))
          O12<-diag(rep(sqrt(sigma2[i,k]/sigma1[i,k]),nrow(dd)))
          
          m2m<-c(m2m,(t(dd$M)%*%inv.O2%*%(dd$M)/(1-delta^2))[1,1])
          m2r.d<-c(m2r.d,(t(dd$M)%*%inv.O2%*%(dd$R-dd$Z*Ct[i,k]-delta*O12%*%(dd$M-dd$Z*At[i,k]))/(1-delta^2))[1,1])
        }
        Bt[i,]<-solve(diag(m2m)+solve(V.beta))%*%(m2r.d+solve(V.beta)%*%rep(b0[2],K))
        
        # C
        z2z=z2r.d<-NULL
        for(k in 1:K)
        {
          dd<-dat[[i]][[k]]
          inv.O1<-diag(rep(1/sigma1[i,k],nrow(dd)))
          inv.O2<-diag(rep(1/sigma2[i,k],nrow(dd)))
          O12<-diag(rep(sqrt(sigma2[i,k]/sigma1[i,k]),nrow(dd)))
          
          z2z<-c(z2z,(t(dd$Z)%*%inv.O2%*%(dd$Z)/(1-delta^2))[1,1])
          z2r.d<-c(z2r.d,(t(dd$Z)%*%inv.O2%*%(dd$R-dd$M*Bt[i,k]-delta*O12%*%(dd$M-dd$Z*At[i,k]))/(1-delta^2))[1,1])
        }
        Ct[i,]<-solve(diag(z2z)+solve(V.gamma))%*%(z2r.d+solve(V.gamma)%*%rep(b0[3],K))
      }
      # update sigma_1ik and sigma_2ik
      if(Sigma.update)
      {
        for(i in 1:N)
        {
          for(k in 1:K)
          {
            dd<-dat[[i]][[k]]
            Y<-cbind(dd$M,dd$R)
            X<-cbind(dd$Z,dd$M)
            Theta<-matrix(c(At[i,k],0,Ct[i,k],Bt[i,k]),2,2)
            S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
            
            sigma1[i,k]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/(nrow(dd)*(1-delta^2))
            sigma2[i,k]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/(nrow(dd)*(1-delta^2))
          }
        }
      }      
      vec.At<-c(At)
      vec.Bt<-c(Bt)
      vec.Ct<-c(Ct)
      fit.A<-lme(vec.At~1,random=~1|Sub,control=lmeControl(opt="optim"))
      fit.B<-lme(vec.Bt~1,random=~1|Sub,control=lmeControl(opt="optim"))
      fit.C<-lme(vec.Ct~1,random=~1|Sub,control=lmeControl(opt="optim"))      
      Afix<-as.numeric(summary(fit.A)$coefficients$fixed[1])
      Bfix<-as.numeric(summary(fit.B)$coefficients$fixed[1])
      Cfix<-as.numeric(summary(fit.C)$coefficients$fixed[1])
      if(random.var.update)
      {
        Phi<-diag(c(as.numeric(VarCorr(fit.A)[1,1]),as.numeric(VarCorr(fit.B)[1,1]),as.numeric(VarCorr(fit.C)[1,1])))
        Lambda<-diag(c(as.numeric(VarCorr(fit.A)[2,1]),as.numeric(VarCorr(fit.B)[2,1]),as.numeric(VarCorr(fit.C)[2,1])))
        if(sum(var.constraint==FALSE)==0)
        {
          # Phi
          if(Phi[1,1]<Phi.confint[1,1])
          {
            Phi[1,1]<-Phi.confint[1,1]
          }else
            if(Phi[1,1]>Phi.confint[1,2])
            {
              Phi[1,1]<-Phi.confint[1,2]
            }
          if(Phi[2,2]<Phi.confint[2,1])
          {
            Phi[2,2]<-Phi.confint[2,1]
          }else
            if(Phi[2,2]>Phi.confint[2,2])
            {
              Phi[2,2]<-Phi.confint[2,2]
            }
          if(Phi[3,3]<Phi.confint[3,1])
          {
            Phi[3,3]<-Phi.confint[3,1]
          }else
            if(Phi[3,3]>Phi.confint[3,2])
            {
              Phi[3,3]<-Phi.confint[3,2]
            }
          
          # Lambda
          if(Lambda[1,1]<Lambda.confint[1,1])
          {
            Lambda[1,1]<-Lambda.confint[1,1]
          }else
            if(Lambda[1,1]>Lambda.confint[1,2])
            {
              Lambda[1,1]<-Lambda.confint[1,2]
            }
          if(Lambda[2,2]<Lambda.confint[2,1])
          {
            Lambda[2,2]<-Lambda.confint[2,1]
          }else
            if(Lambda[2,2]>Lambda.confint[2,2])
            {
              Lambda[2,2]<-Lambda.confint[2,2]
            }
          if(Lambda[3,3]<Lambda.confint[3,1])
          {
            Lambda[3,3]<-Lambda.confint[3,1]
          }else
            if(Lambda[3,3]>Lambda.confint[3,2])
            {
              Lambda[3,3]<-Lambda.confint[3,2]
            }
        }
      }
      u<-cbind(as.matrix(ranef(fit.A)),as.matrix(ranef(fit.B)),as.matrix(ranef(fit.C)))
      colnames(u)<-c("A","B","C")
      
      bnew<-c(Afix,Bfix,Cfix)
      diff<-max(abs(bnew-b0))
      b0<-bnew
      s<-s+1      
    }else
    {
      if((random.var.equal&u.int)|(!random.indep&u.int))
      {
        warning("The u.int=TRUE argument will be ignored and the full h-likelihood method will be applied.")
      }
      bik.mat<-array(NA,c(N,K,3))
      for(i in 1:N)
      {
        for(k in 1:K)
        {
          dd<-dat[[i]][[k]]
          colnames(dd)<-c("Z","M","R")
          n[i,k]<-nrow(dd)
          Z<-matrix(dd[,1])
          M<-matrix(dd[,2])
          R<-matrix(dd[,3])
          
          r<-delta*sqrt(sigma2[i,k]/sigma1[i,k])
          X<-cbind(Z,M)
          Pik<-matrix(c(-r,0,0,1,1,0),2,3)
          Qik<-c(0,r)
          Vik<-matrix(c(1,0,0),1,3)
          bik<-solve(sigma1[i,k]*t(Pik)%*%t(X)%*%X%*%Pik+sigma2[i,k]*(1-delta^2)*t(Vik)%*%t(Z)%*%Z%*%Vik+
                       sigma1[i,k]*sigma2[i,k]*(1-delta^2)*solve(Lambda))%*%
            (sigma1[i,k]*t(Pik)%*%t(X)%*%(R-X%*%Qik)+sigma2[i,k]*(1-delta^2)*t(Vik)%*%t(Z)%*%M+
               sigma1[i,k]*sigma2[i,k]*(1-delta^2)*solve(Lambda)%*%(b0+u[i,]))
          At[i,k]<-bik[1,1]
          Bt[i,k]<-bik[2,1]
          Ct[i,k]<-bik[3,1]
          bik.mat[i,k,]<-bik
        }
      }
      e<-c(sum(At-b0[1])/N,sum(Bt-b0[2])/N,sum(Ct-b0[3])/N)
      for(i in 1:N)
      {
        u[i,]<-solve(solve(Phi)+K*solve(Lambda))%*%solve(Lambda)%*%(apply(cbind(At[i,],Bt[i,],Ct[i,]),2,sum)-K*b0-e)
      }
      bnew<-c(mean(At),mean(Bt),mean(Ct))-apply(u,2,mean)
      if(random.var.update==TRUE)
      {
        Phi<-t(u)%*%u/N
        if(det(Phi)<10e-5|random.indep==TRUE)
        {
          Phi<-diag(as.vector(diag(Phi)))
        }
        if(random.var.equal==FALSE)
        {
          Lambda.new<-matrix(0,3,3)
          for(i in 1:N)
          {
            for(k in 1:K)
            {
              Lambda.new<-Lambda.new+(bik.mat[i,k,]-bnew-u[i,])%*%t(bik.mat[i,k,]-bnew-u[i,])/(N*K)
            }
          }
          Lambda<-diag(as.vector(diag(Lambda.new)))
        }else
        {
          theta2<-0
          for(i in 1:N)
          {
            for(k in 1:K)
            {
              theta2<-theta2+(t(bik.mat[i,k,]-bnew-u[i,])%*%(bik.mat[i,k,]-bnew-u[i,])/(3*N*K))[1,1]
            }
          }
          Lambda<-diag(rep(theta2,3))
        }
        Phi<-t(u)%*%u/N
        if(random.indep)
        {
          Phi<-diag(as.vector(diag(Phi)))
        }
        # add variance constraint
        if(sum(var.constraint==FALSE)==0)
        {
          # Phi
          if(Phi[1,1]<Phi.confint[1,1])
          {
            Phi[1,1]<-Phi.confint[1,1]
          }else
            if(Phi[1,1]>Phi.confint[1,2])
            {
              Phi[1,1]<-Phi.confint[1,2]
            }
          if(Phi[2,2]<Phi.confint[2,1])
          {
            Phi[2,2]<-Phi.confint[2,1]
          }else
            if(Phi[2,2]>Phi.confint[2,2])
            {
              Phi[2,2]<-Phi.confint[2,2]
            }
          if(Phi[3,3]<Phi.confint[3,1])
          {
            Phi[3,3]<-Phi.confint[3,1]
          }else
            if(Phi[3,3]>Phi.confint[3,2])
            {
              Phi[3,3]<-Phi.confint[3,2]
            }
          Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
          Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
          Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
          
          # Lambda
          if(Lambda[1,1]<Lambda.confint[1,1])
          {
            Lambda[1,1]<-Lambda.confint[1,1]
          }else
            if(Lambda[1,1]>Lambda.confint[1,2])
            {
              Lambda[1,1]<-Lambda.confint[1,2]
            }
          if(Lambda[2,2]<Lambda.confint[2,1])
          {
            Lambda[2,2]<-Lambda.confint[2,1]
          }else
            if(Lambda[2,2]>Lambda.confint[2,2])
            {
              Lambda[2,2]<-Lambda.confint[2,2]
            }
          if(Lambda[3,3]<Lambda.confint[3,1])
          {
            Lambda[3,3]<-Lambda.confint[3,1]
          }else
            if(Lambda[3,3]>Lambda.confint[3,2])
            {
              Lambda[3,3]<-Lambda.confint[3,2]
            }
        } 
      }
      # update sigma_1ik and sigma_2ik
      if(Sigma.update)
      {
        for(i in 1:N)
        {
          for(k in 1:K)
          {
            dd<-dat[[i]][[k]]
            Y<-cbind(dd$M,dd$R)
            X<-cbind(dd$Z,dd$M)
            Theta<-matrix(c(At[i,k],0,Ct[i,k],Bt[i,k]),2,2)
            S<-t(Y-X%*%Theta)%*%(Y-X%*%Theta)
            
            sigma1[i,k]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/(nrow(dd)*(1-delta^2))
            sigma2[i,k]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/(nrow(dd)*(1-delta^2))
          }
        }
      }
      
      diff<-max(abs(b0-bnew))
      b0<-bnew
      s<-s+1
    }
  }
  
  HL<-cma.h(dat,delta=delta,A.ik=At,B.ik=Bt,C.ik=Ct,b=b0,u=u,Phi=Phi,Lambda=Lambda,random.indep=random.indep,
            u.int=u.int,Sigma.update=Sigma.update)$h
  
  zc<-qnorm(1-alpha/2)
  coe.re<-matrix(NA,6,4)
  colnames(coe.re)<-c("Estimate","LB","UB","SE")
  rownames(coe.re)<-c("A","C","B","C2","AB.prod","AB.diff")
  # A, C, B, AB.prod, AB.diff
  coe.re[1,4]<-sqrt(Phi[1,1]/N+Lambda[1,1]/(N*K))
  coe.re[2,4]<-sqrt(Phi[3,3]/N+Lambda[3,3]/(N*K))
  coe.re[3,4]<-sqrt(Phi[2,2]/N+Lambda[2,2]/(N*K))
  coe.re[1,1:3]<-c(b0[1],b0[1]-zc*coe.re[1,4],b0[1]+zc*coe.re[1,4])
  coe.re[2,1:3]<-c(b0[3],b0[3]-zc*coe.re[2,4],b0[3]+zc*coe.re[2,4])
  coe.re[3,1:3]<-c(b0[2],b0[2]-zc*coe.re[3,4],b0[2]+zc*coe.re[3,4])
  coe.re[5,1]<-b0[1]*b0[2]
  coe.re[5,4]<-sqrt((b0[1]*coe.re[3,4])^2+(b0[2]*coe.re[1,4])^2)
  coe.re[5,2:3]<-c(coe.re[5,1]-zc*coe.re[5,4],coe.re[5,1]+zc*coe.re[5,4])
  # C2, AB.diff
  if(optimizer[1]=="optimx")
  {
    fit.C2<-lmer(c(C2t)~(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
  }else
  {
    fit.C2<-lmer(c(C2t)~(1|Sub),control=lmerControl(optimizer=optimizer[1]))
  }
  C2fix<-summary(fit.C2)$coefficients[1]
  s2.C2<-c(as.numeric(VarCorr(fit.C2)),attr(VarCorr(fit.C2),"sc")^2)
  e.C2<-C2fix+mean(ranef(fit.C2)$Sub[,"(Intercept)"])
  coe.re[4,4]<-sqrt(s2.C2[1]/N+s2.C2[2]/(N*K))
  coe.re[4,1:3]<-c(C2fix,C2fix-zc*coe.re[4,4],C2fix+zc*coe.re[4,4])
  coe.re[6,1]<-coe.re[4,1]-coe.re[2,1]
  coe.re[6,4]<-sqrt(coe.re[4,4]^2+coe.re[2,4]^2)
  coe.re[6,2:3]<-c(coe.re[6,1]-zc*coe.re[6,4],coe.re[6,1]+zc*coe.re[6,4])  
  s2.C2<-data.frame(C2=s2.C2[1],Error=s2.C2[2],Total=sum(s2.C2))
  colnames(s2.C2)<-c("Random(C')","Var(Error)","Var(C')")
  
  var.comp<-data.frame(delta=delta,A=Phi[1,1],C=Phi[3,3],B=Phi[2,2],Lambda[1,1],Lambda[3,3],Lambda[2,2])
  colnames(var.comp)<-c("delta","Random(A)","Random(C)","Random(B)","Var(Error A)","Var(Error C)","Var(Error B)")
  
  if(random.indep)
  {
    Phi<-diag(as.vector(diag(Phi)))
  }
  cor.comp.temp<-sqrt(diag(1/diag(Phi)))%*%Phi%*%sqrt(diag(1/diag(Phi)))
  cor.comp<-cor.comp.temp
  cor.comp[1,2]=cor.comp[2,1]<-cor.comp.temp[1,3]
  cor.comp[1,3]=cor.comp[3,1]<-cor.comp.temp[1,2]
  colnames(cor.comp)=rownames(cor.comp)<-c("A","C","B")
  
  re<-list(delta=delta,Coefficients=coe.re,Cor.comp=cor.comp,Var.comp=var.comp,Var.C2=s2.C2,logLik=LL,HL=HL)
  if(sum(var.constraint==FALSE)==0)
  {
    re.var.constraint<-rbind(Phi.confint,Lambda.confint)
    rownames(re.var.constraint)<-c("Random(A)","Random(B)","Random(C)","Var(Error A)","Var(Error B)","Var(Error C)")
    re<-list(delta=delta,Coefficients=coe.re,Cor.comp=cor.comp,Var.comp=var.comp,Var.constraint=re.var.constraint,
             Var.C2=s2.C2,logLik=LL,HL=HL,convergence=(s<max.itr|max.itr==0))
  }
  return(re)
}

### ---- from macc/macc/R/cma.uni.mix.R ----
cma.uni.mix <-
function(dat,delta=0,conf.level=0.95,optimizer=c("bobyqa","Nelder_Mead","optimx"),
                      mix.pkg=c("lme4","nlme"),random.indep=FALSE,random.var.equal=TRUE,u.int=FALSE)
{
  N<-length(dat)
  K<-length(dat[[1]])
  
  coe<-rep(c("A","C","B"),each=N*K)
  Sub<-factor(rep(1:N,K))
  Sub.all<-rep(Sub,3)
  
  zc<-qnorm(1-(1-conf.level)/2)
  
  coe.re<-matrix(NA,6,4)
  colnames(coe.re)<-c("Estimate","LB","UB","SE")
  rownames(coe.re)<-c("A","C","B","C2","AB.prod","AB.diff")
  
  # Estimate A and C'
  At<-matrix(NA,length(dat),length(dat[[1]]))
  colnames(At)<-paste("Session",1:K)
  C2t<-At
  for(i in 1:length(dat))
  {
    for(k in 1:length(dat[[i]]))
    {
      dd<-dat[[i]][[k]]
      re<-cma.uni.delta(dd,delta=0)
      
      At[i,k]<-re$D[1,1]
      C2t[i,k]<-re$Coefficients[4,1]
    }
  }
  if(optimizer[1]=="optimx")
  {
    fit.C2<-lmer(c(C2t)~(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
  }else
  {
    fit.C2<-lmer(c(C2t)~(1|Sub),control=lmerControl(optimizer=optimizer[1]))
  }
  C2fix<-summary(fit.C2)$coefficients[1]
  s2.C2<-c(as.numeric(VarCorr(fit.C2)),attr(VarCorr(fit.C2),"sc")^2)
  e.C2<-C2fix+mean(ranef(fit.C2)$Sub[,"(Intercept)"])
  coe.re[4,4]<-sqrt(s2.C2[1]/N+s2.C2[2]/(N*K))
  coe.re[4,1:3]<-c(C2fix,C2fix-zc*coe.re[4,4],C2fix+zc*coe.re[4,4])
  
  # Given value of delta, estimate B and C  
  Ct<-matrix(NA,length(dat),length(dat[[1]]))
  colnames(Ct)<-paste("Session",1:K)
  Bt<-Ct
  for(i in 1:length(dat))
  {
    for(k in 1:length(dat[[i]]))
    {
      dd<-dat[[i]][[k]]
      re<-cma.uni.delta(dd,delta=delta)
      
      Bt[i,k]<-re$D[2,2]
      Ct[i,k]<-re$D[1,2]
    }
  }
  coet<-c(c(At),c(Ct),c(Bt))  
  
  if(!random.var.equal)
  {
    if(random.indep)
    {
      vec.At<-c(At)
      vec.Bt<-c(Bt)
      vec.Ct<-c(Ct)
      if(mix.pkg[1]=="nlme")
      {
        if(optimizer[1]=="optimx")
        {
          fit.A<-lme(vec.At~1,random=~1|Sub,control=lmeControl(opt="optim"))
          fit.B<-lme(vec.Bt~1,random=~1|Sub,control=lmeControl(opt="optim"))
          fit.C<-lme(vec.Ct~1,random=~1|Sub,control=lmeControl(opt="optim"))
        }else
        {
          fit.A<-lme(vec.At~1,random=~1|Sub)
          fit.B<-lme(vec.Bt~1,random=~1|Sub)
          fit.C<-lme(vec.Ct~1,random=~1|Sub)
        }
        Afix<-as.numeric(summary(fit.A)$coefficients$fixed[1])
        Bfix<-as.numeric(summary(fit.B)$coefficients$fixed[1])
        Cfix<-as.numeric(summary(fit.C)$coefficients$fixed[1])
        Phi<-diag(c(as.numeric(VarCorr(fit.A)[1,1]),as.numeric(VarCorr(fit.B)[1,1]),as.numeric(VarCorr(fit.C)[1,1])))
        Lambda<-diag(c(as.numeric(VarCorr(fit.A)[2,1]),as.numeric(VarCorr(fit.B)[2,1]),as.numeric(VarCorr(fit.C)[2,1])))
        u<-cbind(as.matrix(ranef(fit.A)),as.matrix(ranef(fit.B)),as.matrix(ranef(fit.C)))
        colnames(u)<-c("A","B","C")
        b0<-c(Afix,Bfix,Cfix)
        cor.AB=cor.AC=cor.BC<-0
        LL<-as.numeric(logLik(fit.A,REML=FALSE)+logLik(fit.B,REML=FALSE)+logLik(fit.C,REML=FALSE))
      }else
      {
        if(optimizer[1]=="optimx")
        {
          fit.A<-lmer(vec.At~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          fit.B<-lmer(vec.Bt~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          fit.C<-lmer(vec.Ct~1+(1|Sub),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        }else
        {
          fit.A<-lmer(vec.At~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
          fit.B<-lmer(vec.Bt~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
          fit.C<-lmer(vec.Ct~1+(1|Sub),control= lmerControl(optimizer=optimizer[1]))
        }
        Afix<-as.numeric(summary(fit.A)$coefficients[1])
        Bfix<-as.numeric(summary(fit.B)$coefficients[1])
        Cfix<-as.numeric(summary(fit.C)$coefficients[1])
        Phi<-diag(c(as.numeric((attr(VarCorr(fit.A)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit.B)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit.C)[[1]],"stddev")^2)[1])))
        Lambda<-diag(c(as.numeric((attr(VarCorr(fit.A),"sc")^2)[1]),
                       as.numeric((attr(VarCorr(fit.B),"sc")^2)[1]),
                       as.numeric((attr(VarCorr(fit.C),"sc")^2)[1])))
        u<-cbind(as.matrix(ranef(fit.A)$Sub),as.matrix(ranef(fit.B)$Sub),as.matrix(ranef(fit.C)$Sub))
        colnames(u)<-c("A","B","C")
        b0<-c(Afix,Bfix,Cfix)
        cor.AB=cor.AC=cor.BC<-0
        LL<-as.numeric(logLik(fit.A,REML=FALSE)+logLik(fit.B,REML=FALSE)+logLik(fit.C,REML=FALSE))
      }
    }else
    {
      if(mix.pkg[1]=="nlme")
      {
        if(optimizer[1]=="optimx")
        {
          fit0<-lme(c(coet)~0+coe,random=~0+coe|Sub.all,control=lmeControl(opt="optim"))
        }else
        {
          fit0<-lme(c(coet)~0+coe,random=~0+coe|Sub.all)
        }
        fit<-NULL
        try(fit<-update(fit0,weights=varIdent(form=~1|factor(coe))))
        if(is.null(fit))
        {
          warning("Equal-varaince assumption in random error is applied instead.")
          fit<-fit0
          Lambda<-diag(rep(summary(fit)$sigma^2,3))
        }else
        {
          Lambda<-diag(c((summary(fit)$sigma*coef(fit$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE)[1])^2,
                         (summary(fit)$sigma*coef(fit$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE)[3])^2,
                         (summary(fit)$sigma*coef(fit$modelStruct$varStruct, unconstrained=FALSE, allCoef=TRUE)[2])^2))
        }
        Afix<-as.numeric(summary(fit)$coefficients$fixed[1])
        Bfix<-as.numeric(summary(fit)$coefficients$fixed[2])
        Cfix<-as.numeric(summary(fit)$coefficients$fixed[3])
        b0<-c(Afix,Bfix,Cfix)
        u<-as.matrix(ranef(fit))
        Phi<-diag(c(as.numeric(VarCorr(fit)[1,1]),as.numeric(VarCorr(fit)[2,1]),as.numeric(VarCorr(fit)[3,1])))
        cor.AB<-as.numeric(VarCorr(fit)[2,3])
        cor.AC<-as.numeric(VarCorr(fit)[3,3])
        cor.BC<-as.numeric(VarCorr(fit)[3,4])
        Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
        Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
        Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
        LL<-as.numeric(logLik(fit,REML=FALSE))
      }else
      {
        warning("Equal-varaince assumption in random error is applied instead.")
        if(optimizer[1]=="optimx")
        {
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        }else
        {
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer=optimizer[1]))
        }
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
        b0<-c(summary(fit)$coefficients[1,1],summary(fit)$coefficients[2,1],summary(fit)$coefficients[3,1])
        Phi<-diag(c(as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[1]),
                    as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[2]),
                    as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[3])))
        cor.AB<-attr(VarCorr(fit)[[1]],"correlation")[1,2]
        cor.AC<-attr(VarCorr(fit)[[1]],"correlation")[1,3]
        cor.BC<-attr(VarCorr(fit)[[1]],"correlation")[2,3]
        Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
        Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
        Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
        Lambda<-diag(rep(attr(VarCorr(fit),"sc")^2,3))
        LL<-as.numeric(logLik(fit,REML=FALSE))
      }
    }
  }else
  {
    if(mix.pkg[1]=="nlme")
    {
      if(optimizer[1]=="optimx")
      {
        fit<-lme(c(coet)~0+coe,random=~0+coe|Sub.all,control=lmeControl(opt="optim"))        
      }else
      {
        fit<-lme(c(coet)~0+coe,random=~0+coe|Sub.all)
      }
      Afix<-as.numeric(summary(fit)$coefficients$fixed[1])
      Bfix<-as.numeric(summary(fit)$coefficients$fixed[2])
      Cfix<-as.numeric(summary(fit)$coefficients$fixed[3])
      b0<-c(Afix,Bfix,Cfix)
      u<-as.matrix(ranef(fit))
      Lambda<-diag(rep(summary(fit)$sigma^2,3))
      Phi<-diag(c(as.numeric(VarCorr(fit)[1,1]),as.numeric(VarCorr(fit)[2,1]),as.numeric(VarCorr(fit)[3,1])))
      if(random.indep)
      {
        cor.AB=cor.AC=cor.BC<-0
      }else
      {
        cor.AB<-as.numeric(VarCorr(fit)[2,3])
        cor.AC<-as.numeric(VarCorr(fit)[3,3])
        cor.BC<-as.numeric(VarCorr(fit)[3,4])
      }
      Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
      Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
      Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
      LL<-as.numeric(logLik(fit,REML=FALSE))
    }else
    {
      if(optimizer[1]=="optimx")
      {
        fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
      }else
      {
        fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer=optimizer[1]))
        u<-as.matrix(ranef(fit)$Sub.all)     
        cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        s<-0
        while(length(which(abs(cor.t)<1e-06))>0&s<20)
        {
          s<-s+1
          fit<-lmer(c(coet)~0+coe+(0+coe|Sub.all),control= lmerControl(optimizer="optimx", optCtrl=list(method="L-BFGS-B")))
          cor.t<-1-c(attr(VarCorr(fit)[[1]],"correlation")-diag(diag(attr(VarCorr(fit)[[1]],"correlation"))))
        }
      }
      b0<-c(summary(fit)$coefficients[1,1],summary(fit)$coefficients[2,1],summary(fit)$coefficients[3,1])
      Phi<-diag(c(as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[1]),
                  as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[2]),
                  as.numeric((attr(VarCorr(fit)[[1]],"stddev")^2)[3])))
      if(random.indep)
      {
        cor.AB=cor.AC=cor.BC<-0
      }else
      {
        cor.AB<-attr(VarCorr(fit)[[1]],"correlation")[1,2]
        cor.AC<-attr(VarCorr(fit)[[1]],"correlation")[1,3]
        cor.BC<-attr(VarCorr(fit)[[1]],"correlation")[2,3]
      }
      Phi[1,2]=Phi[2,1]<-cor.AB*sqrt(Phi[1,1]*Phi[2,2])
      Phi[1,3]=Phi[3,1]<-cor.AC*sqrt(Phi[1,1]*Phi[3,3])
      Phi[2,3]=Phi[3,2]<-cor.BC*sqrt(Phi[2,2]*Phi[3,3])
      Lambda<-diag(rep(attr(VarCorr(fit),"sc")^2,3))
      LL<-as.numeric(logLik(fit,REML=FALSE))
    }
  }
  
  cor.comp<-matrix(NA,3,3)
  colnames(cor.comp)<-c("A","C","B")
  rownames(cor.comp)<-c("A","C","B")
  if(random.indep==TRUE)
  {
    cor.AC=cor.AB=cor.BC<-0
  }
  cor.comp[1,2]=cor.comp[2,1]<-cor.AC
  cor.comp[1,3]=cor.comp[3,1]<-cor.AB
  cor.comp[2,3]=cor.comp[3,2]<-cor.BC
  diag(cor.comp)<-rep(1,3)
  
  delta.est<-delta
  coe.re[c(1,3,2),4]<-sqrt(diag(Phi)/N+diag(Lambda)/(N*K))
  coe.re[c(1,3,2),1]<-b0
  coe.re[c(1,3,2),2]<-b0-zc*coe.re[c(1,3,2),4]
  coe.re[c(1,3,2),3]<-b0+zc*coe.re[c(1,3,2),4]  
  coe.re[5,1]<-b0[1]*b0[2]
  coe.re[5,4]<-sqrt((b0[1]*coe.re[3,4])^2+(b0[2]*coe.re[1,4])^2)
  coe.re[5,2:3]<-c(coe.re[5,1]-zc*coe.re[5,4],coe.re[5,1]+zc*coe.re[5,4])
  coe.re[6,1]<-C2fix-b0[3]
  coe.re[6,4]<-sqrt(coe.re[4,4]^2+coe.re[2,4]^2)
  coe.re[6,2:3]<-c(coe.re[6,1]-zc*coe.re[6,4],coe.re[6,1]+zc*coe.re[6,4])  
  
  var.comp<-data.frame(delta=delta,A=Phi[1,1],C=Phi[3,3],B=Phi[2,2],Lambda[1,1],Lambda[3,3],Lambda[2,2])
  colnames(var.comp)<-c("delta","Random(A)","Random(C)","Random(B)","Var(Error A)","Var(Error C)","Var(Error B)")
  
  s2.C2<-data.frame(C2=s2.C2[1],Error=s2.C2[2],Total=sum(s2.C2))
  colnames(s2.C2)<-c("Random(C')","Var(Error)","Var(C')")
  
  HL<-cma.h(dat,delta=delta,A.ik=At,B.ik=Bt,C.ik=Ct,b=b0,u=u,Phi=Phi,Lambda=Lambda,random.indep=random.indep,u.int=u.int)
  
  re<-list(delta=delta.est,Coefficients=coe.re,Cor.comp=cor.comp,Var.comp=var.comp,Var.C2=s2.C2,logLik=LL,HL=HL$h)
  
  return(re)
}

### ---- from macc/macc/R/cma.uni.plot.R ----
cma.uni.plot <-
function(re.cma.sens,re.cma=NULL,delta=NULL,legend.pos="topright",xlab=expression(delta),ylab=expression(hat(AB)),
                       cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
{
  dt<-re.cma.sens$Coefficients[,"delta"]
  AB.p<-re.cma.sens$Coefficients[,"ABp.Estimate"]
  AB.p.ub<-re.cma.sens$Coefficients[,"ABp.UB"]
  AB.p.lb<-re.cma.sens$Coefficients[,"ABp.LB"]
  idx<-sort(AB.p,index.return=TRUE)$ix
  
  #################################################
  # plot delta = 0 result or not
  if(!is.null(re.cma))
  {
    if(re.cma$delta==0)
    {
      plot.delta0<-FALSE
    }
  }
  if(!is.null(delta))
  {
    if(length(which(delta==0))>0)
    {
      plot.delta0<-FALSE
    }
  }
  #################################################
  
  plot(range(dt[idx]),range(c(AB.p.lb,AB.p.ub)),type="n",xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis)
  polygon(c(rev(dt),dt),c(rev(AB.p.ub),AB.p.lb),col="grey80",border=NA)
  abline(v=0)
  abline(h=0)
  if(plot.delta0)
  {
    abline(h=AB.p[which(dt==0)],lty=2,col=2) 
  }
  lines(dt,AB.p,lwd=2)
  lines(dt,AB.p.lb,lty=2,col=8)
  lines(dt,AB.p.ub,lty=2,col=8)
  
  if(!is.null(re.cma))
  {
    points(re.cma$delta,re.cma$Coefficients["ABp",1],pch=16,col=4,cex=0.75)
    lines(rep(re.cma$delta,2),re.cma$Coefficients["ABp",c(3,4)],lty=2,col=4)
    lines(c(re.cma$delta-0.02,re.cma$delta+0.02),rep(re.cma$Coefficients["ABp",3],2),col=4)
    lines(c(re.cma$delta-0.02,re.cma$delta+0.02),rep(re.cma$Coefficients["ABp",4],2),col=4)
    
    if(plot.delta0)
    {
      legend(legend.pos,legend=c(expression(delta==0),
                                 substitute(delta==d,list(d=round(re.cma$delta,digits=3)))),
             lty=2,col=c(2,4),pch=c(NA,16),pt.cex=lgd.pt.cex,bty="n",cex=lgd.cex) 
    }else
    {
      legend(legend.pos,legend=substitute(delta==d,list(d=round(re.cma$delta,digits=3))),
             lty=2,col=4,pch=16,pt.cex=lgd.pt.cex,bty="n",cex=lgd.cex)
    }
  }
  if(!is.null(delta))
  {
    for(j in 1:length(delta))
    {
      idx.tmp<-which.min(abs(dt-delta[j]))
      
      points(dt[idx.tmp],AB.p[idx.tmp],pch=16,col=4,cex=0.75)
      lines(rep(dt[idx.tmp],2),c(AB.p.lb[idx.tmp],AB.p.ub[idx.tmp]),lty=2,col=4)
      lines(c(dt[idx.tmp]-0.02,dt[idx.tmp]+0.02),rep(AB.p.lb[idx.tmp],2),col=4)
      lines(c(dt[idx.tmp]-0.02,dt[idx.tmp]+0.02),rep(AB.p.ub[idx.tmp],2),col=4)
    }
    
    if(plot.delta0)
    {
      legend(legend.pos,legend=expression(delta==0),lty=2,col=2,bty="n",cex=lgd.cex) 
    }
  }
}

### ---- from macc/macc/R/cma.uni.sens.R ----
cma.uni.sens <-
function(dat,delta=seq(-1,1,by=0.01),conf.level=0.95)
{
  Z<-matrix(dat$Z,ncol=1)
  M<-matrix(dat$M,ncol=1)
  R<-matrix(dat$R,ncol=1)
  
  n<-nrow(Z)
  
  Coe<-NULL
  Sp<-NULL
  dt<-NULL
  for(i in 1:length(delta))
  {
    if(abs(delta[i])!=1)
    {
      dt<-c(dt,delta[i])
      re<-cma.uni.delta(dat,delta=delta[i],conf.level=conf.level)
      table1<-re$Coefficients
      Coe<-rbind(Coe,c(t(table1)))
      colnames(Coe)<-paste(rep(rownames(table1),each=ncol(table1)),
                           rep(colnames(table1),nrow(table1)),sep=".")
      Sp<-rbind(Sp,c(re$Sigma[1,1],re$Sigma[2,2],re$Sigma[1,2]))
      colnames(Sp)<-c("sigma1^2","sigma2^2","rho")
    }  
  }
  return(list(Coefficients=cbind(delta=dt,Coe),Covariance=cbind(Sp,delta=dt)))
}

### ---- from macc/macc/R/macc-internal.R ----
.Random.seed <-
c(403L, 176L, -1544986965L, 1967274764L, 1370809523L, -1074270084L, 
1272176128L, 1030583483L, -1010711648L, 1059138151L, -483154787L, 
-1776473706L, -1701613475L, 1435703338L, 947425862L, -1590223903L, 
-1723459570L, -1972862455L, 1179771015L, 957296456L, 1902093711L, 
-1770282344L, 795246644L, 140161391L, 254545204L, -139580333L, 
340226129L, -1804441190L, -269308751L, 1937323550L, 124920034L, 
-1762164355L, 1959266330L, -1988657731L, 757882307L, 1543171908L, 
1724862363L, -651897436L, -1781537528L, -810717901L, 1924213624L, 
1177794447L, 1998184613L, -1716173970L, 2065772549L, 800161266L, 
1295019646L, 1650855449L, 1691246950L, -105684031L, -1481477425L, 
-1934196272L, -1331513945L, 970422288L, 1531325900L, 1382968055L, 
-1870928068L, 1791440843L, 897475817L, -439589136L, 1554042565L, 
-813989729L, -301180813L, 1943254945L, 1756319316L, 155935222L, 
-740163664L, -252403230L, -1340738149L, -1516896647L, -787872547L, 
-729218533L, -408997102L, -1258169720L, 1208143906L, 1550622044L, 
392603017L, 1795007883L, 742986023L, 530348045L, -1337030544L, 
1459780850L, -1125833692L, -259591234L, 88464111L, 45930357L, 
-468697183L, -107122473L, -1344980802L, 1375301300L, -867239530L, 
-1434586648L, 1297476253L, 899279447L, 847375451L, -2048726695L, 
-1493358804L, -1455098706L, 23862632L, -626004870L, -1351406461L, 
1911848497L, -129902427L, -1281305949L, 20038938L, -2071130352L, 
1201539306L, -587269740L, 1044993441L, 1366121539L, 164385007L, 
415537733L, -122823576L, -1399188470L, 1933428220L, 1504572982L, 
885128439L, 1388573709L, -1126656567L, 1102945151L, -2052812634L, 
-67546052L, -673437890L, -1601916640L, -1537756331L, -893183857L, 
-339049149L, 506675601L, 690486532L, 913459078L, 2127487616L, 
1287673682L, 101601195L, 786767273L, 1916743181L, 549664267L, 
-839446366L, -1488557128L, -1615626286L, 631698348L, 1447016345L, 
-2021220197L, -292902857L, -1881794979L, 1185699072L, 839774082L, 
-586122444L, -980581266L, -1457234561L, -1790112603L, -228373071L, 
-1736174489L, -997956690L, 1752239012L, -604534650L, 1340776024L, 
327446509L, 1591241223L, 1881170411L, -2000063159L, 104953948L, 
-1972938754L, 898228856L, -342481238L, -1739987501L, -1808372831L, 
-1344875563L, 484206995L, 1716467562L, -486208320L, 1724427098L, 
-42196060L, -1187597071L, -253567853L, 1715683583L, 1006507093L, 
-1746448904L, 256997658L, 810042828L, 209956518L, -369014969L, 
-410927107L, 1880201241L, -1031188401L, -762054698L, -129622292L, 
1362699630L, 738609616L, 1186211813L, 233046719L, -1525722221L, 
648546369L, -1478406540L, -472288426L, 48030032L, 446945794L, 
-1786564933L, 984687833L, 565407997L, 1889846715L, -23695758L, 
-268839960L, -620916414L, 241704892L, -18549655L, 411729451L, 
-1843383993L, -686672531L, 1653863504L, 1960054354L, -1227926588L, 
-389007394L, -366963313L, 161671509L, -1201947007L, 1743980599L, 
-1745432418L, -1879525868L, 1096665462L, 1200614728L, 847420157L, 
756652791L, -998071685L, -1006029319L, -969777716L, -1124046258L, 
-1698547320L, 1417926490L, -1342256093L, -196039919L, 1197578949L, 
1119019395L, -1653913670L, 1545934256L, 879720778L, -1888688524L, 
2714625L, -1979670429L, -711127473L, -1883521371L, 1984136904L, 
2025253610L, -479942308L, -2012277034L, -192969769L, 1628157101L, 
72655017L, 427776415L, -1070626106L, 70515996L, 1903196830L, 
1824678400L, -809131147L, 2107180079L, -1870659997L, -1413616079L, 
1130040228L, 1996332262L, -449329760L, 672100082L, -1644686261L, 
2005305353L, 1169500845L, 92637227L, -1731598078L, -312756840L, 
-2011030798L, 1577933836L, 638635001L, 118816699L, -463318313L, 
898999101L, 734854112L, 583549410L, 2098851412L, 1968960014L, 
-1981185633L, -917051515L, -1492247407L, -789099065L, 2099886734L, 
1295165828L, 429295078L, -1839658440L, 477470029L, 63026599L, 
-985666055L, 1034782213L, -412286571L, -790210713L, 522176316L, 
1971536828L, 245708480L, -1495735766L, -1652733781L, -1368455381L, 
-1469340449L, -292457171L, -773587510L, -1177436986L, -643392178L, 
1616125308L, 1171587301L, 1188271601L, 662270889L, 750024587L, 
1515454480L, -1624348216L, 1342097132L, -1094454866L, -281848481L, 
444021919L, 388836555L, -1239251767L, -259756666L, -238734678L, 
1898466866L, -1436103608L, -137337999L, 1393647997L, 658401613L, 
1711536575L, -151352572L, -310389532L, -2089394088L, 929191746L, 
-133902029L, -261576429L, 1032224471L, -2135640747L, -1142269502L, 
926580238L, -111724474L, 866638772L, -1746769555L, 376404281L, 
1550141457L, 96767411L, 635368936L, -233557344L, -924110156L, 
219737462L, 978333847L, 700208503L, 496185811L, -1644651423L, 
107870702L, -493681662L, -1105894854L, -922429168L, -820870327L, 
-1280216747L, -1255609563L, -1368077897L, 1696533612L, 629377132L, 
1923183888L, 2110569722L, 1356506779L, 2126563739L, 851611951L, 
1466041213L, 1219396794L, -521905162L, -385657698L, -1428961684L, 
1617923125L, 918804129L, -342218695L, -246068261L, -398543776L, 
1082040344L, 757594332L, -782668642L, -1681351857L, -784883025L, 
-1593670917L, -1499433639L, -1226189610L, -1181073318L, 1227900482L, 
-1838991368L, 491859873L, 1870599597L, 2040601213L, 1484851631L, 
-722547244L, 391967220L, -410294328L, 346659890L, -277302269L, 
-608475997L, -1074471481L, -307935867L, 1310134098L, 1180821790L, 
-556864842L, 1239342916L, -2029011555L, -1258747191L, -1737421439L, 
1166310371L, 1349383512L, -380688496L, 1549212164L, 1430654342L, 
9724135L, -26301785L, -132507805L, 890541905L, -367662498L, -2038110254L, 
-256086870L, 361498464L, 1089814937L, 573851045L, 1511474357L, 
-1090951161L, 1653945436L, -398263780L, -1910610080L, 137109002L, 
2088262667L, -175532213L, 910120255L, 1870047245L, 1098068394L, 
-1361220378L, -220120210L, -461483940L, 1636723653L, -1440538095L, 
1379202441L, -1463330965L, -1559910800L, 1377893416L, 1534508300L, 
-1651567922L, -251616833L, -1144633153L, 349821419L, -1828648727L, 
1934400550L, -878005558L, -428556590L, 1687953832L, 726242129L, 
478040413L, 18371501L, -1944269601L, -37648732L, -1496001468L, 
-180791752L, -303122526L, 1191193299L, -176662134L, -169810239L, 
-597762722L, -1237058559L, -1513666756L, 636185597L, -1840969582L, 
523859335L, 724418934L, -1994372899L, 581086822L, 1636653361L, 
1921905272L, 1616055413L, 613594758L, 1866380111L, -1356369910L, 
477781585L, -1609733978L, -1797299479L, -2048249348L, 1168666909L, 
-1496749302L, -163482385L, -560871642L, 1724315773L, 741262870L, 
1012483601L, 1427775464L, 1254532669L, 674445918L, -1518592017L, 
-1622499238L, 1523383681L, -382490322L, -1269848767L, -681095124L, 
-873016739L, 1758733714L, -1940204105L, 1622881878L, 974569405L, 
268629478L, -1079898479L, 242442856L, 936686837L, -183634202L, 
-1561492465L, -336372726L, -1013092767L, 602317878L, 1329465897L, 
234604476L, 562849101L, -976675318L, 1499070975L, 1725335798L, 
1231181293L, 1918171764L, -1378645235L, -936170079L, -1091894972L, 
-463576294L, -1766647351L, 124856665L, -1171518550L, 1712554880L, 
1751978713L, -1172894071L, -1433477008L, 1286671650L, 264184549L, 
-1754982135L, 95710674L, -385139492L, -163619915L, 1870185545L, 
1407177668L, -1668151686L, 1447306081L, 1349411905L, -372850014L, 
1676704832L, -379807527L, 363229025L, -1811155848L, 1466227106L, 
1172072277L, 1570472273L, -1949090742L, 500073588L, -1857748499L, 
365258545L, 1945910228L, 1026089322L, 1683891369L, -1238099143L, 
1093381786L, 2124637696L, 1575370089L, 1928292617L, -1867487584L, 
678187282L, -499157627L, -1390891687L, 350093346L, 2062754060L, 
-112443899L, 1565557577L, 1847265236L, -1554928694L, 308216145L, 
-1524925839L, 205380002L, -591109184L, -1046384983L, 1172445521L, 
-4067592L, -705998606L, 1310307381L, -549376111L, -1100979830L, 
386374580L, 1757515213L, -314683359L, -1388018524L, -1140266534L, 
-47692119L, 794649817L, 873282890L, -934426336L, -118437447L, 
903780297L, 644541776L, -1279658814L, 691981093L, -1167276151L, 
-1068351374L, -2028926916L, -2108673195L, -1503671831L, -1363787868L, 
563095898L, -1178685919L, 979431585L, -1817590526L, 1201497056L, 
1406170041L, -1700646623L, -2055595208L, 842054754L, 562468373L, 
1540454449L, -1982186742L, -299592204L, -390573075L, 1850186257L, 
571324884L, 129351818L, -68868343L, -1865198087L, 1177654746L, 
1584908448L, 1176587625L, -1188822263L, -112187552L, -1872452398L, 
-1511561787L, -817570247L, -1310399423L, 1341366807L)

### ---- from macc/macc/R/macc.R ----
macc <-
function(dat,model.type=c("single","multilevel","twolevel"),method=c("HL","TS","HL-TS"),
               delta=NULL,interval=c(-0.90,0.90),tol=10e-4,max.itr=500,conf.level=0.95,
               optimizer=c("optimx","bobyqa","Nelder_Mead"),mix.pkg=c("nlme","lme4"),
               random.indep=TRUE,random.var.equal=FALSE,u.int=FALSE,Sigma.update=TRUE,
               var.constraint=TRUE,random.var.update=TRUE,logLik.type=c("logLik","HL"),
               error.indep=TRUE,error.var.equal=FALSE,
               sens.plot=FALSE,sens.interval=seq(-1,1,by=0.01),legend.pos="topright",
               xlab=expression(delta),ylab=expression(hat(AB)),
               cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
{
  if(model.type[1]=="single")
  {
    if(is.null(delta))
    {
      delta<-0
    }
    run.time<-system.time(re<-cma.uni.delta(dat,delta=delta,conf.level=conf.level))
    if(sens.plot)
    {
      re.sens<-cma.uni.sens(dat,delta=sens.interval,conf.level=conf.level)
      cma.uni.plot(re.sens,re,delta=NULL,legend.pos=legend.pos,xlab=xlab,ylab=ylab,
                   cex.lab=cex.lab,cex.axis=cex.axis,lgd.cex=lgd.cex,lgd.pt.cex=1,plot.delta0=plot.delta0,...)
    }
  }else
    if(model.type[1]=="multilevel")
    {
      if(is.null(delta))
      {
        if(method[1]=="TS")
        {
          run.time1<-system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=0,optimizer=optimizer,
                                               mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                               u.int=u.int,Sigma.update=Sigma.update,logLik.type=logLik.type,maximum=TRUE))
          run.time2<-system.time(re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                 random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
          
          run.time<-run.time1+run.time2
        }else
        {
          run.time1<-system.time(re1<-optimize(cma.uni.mix.dhl,interval=interval,dat=dat,tol=tol,max.itr=max.itr,optimizer=optimizer,
                                               mix.pkg=mix.pkg,random.indep=random.indep,random.var.equal=random.var.equal,
                                               u.int=u.int,Sigma.update=Sigma.update,var.constraint=var.constraint,
                                               random.var.update=random.var.update,logLik.type=logLik.type,maximum=TRUE))
          if(method[1]=="HL-TS")
          {
            run.time2<-system.time(re<-cma.uni.mix(dat,delta=re1$maximum,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                   random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
          }
          if(method[1]=="HL")
          {
            run.time2<-system.time(re<-cma.uni.mix.hl(dat,delta=re1$maximum,tol=tol,max.itr=max.itr,alpha=1-conf.level,
                                                      random.indep=random.indep,optimizer=optimizer,mix.pkg=mix.pkg,
                                                      random.var.equal=random.var.equal,u.int=u.int,
                                                      Sigma.update=Sigma.update,var.constraint=var.constraint,
                                                      random.var.update=random.var.update))
          }
          
          run.time<-run.time1+run.time2
        }
      }else
      {
        if(method[1]=="TS")
        {
          run.time<-system.time(re<-cma.uni.mix(dat,delta=delta,conf.level=conf.level,optimizer=optimizer,mix.pkg=mix.pkg,
                                                random.indep=random.indep,random.var.equal=random.var.equal,u.int=u.int))
        }
        if(method[1]=="HL")
        {
          run.time<-system.time(re<-cma.uni.mix.hl(dat,delta=delta,tol=tol,max.itr=max.itr,alpha=1-conf.level,random.indep=random.indep,
                                                   optimizer=optimizer,mix.pkg=mix.pkg,random.var.equal=random.var.equal,u.int=u.int,
                                                   Sigma.update=Sigma.update,var.constraint=var.constraint,
                                                   random.var.update=random.var.update))
        }
      }
    }else
      if(model.type[1]=="twolevel")
      {
        if(is.null(delta))
        {
          if(method[1]=="TS")
          {
            run.time1<-system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=0,tol=tol,
                                                 error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                 var.constraint=var.constraint,maximum=TRUE))
            run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,
                                                    error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                    var.constraint=var.constraint))
            
            run.time<-run.time1+run.time2
          }else
          {
            run.time1<-system.time(re1<-optimize(cma.delta.lm.HL,interval=interval,dat=dat,max.itr=max.itr,tol=tol,
                                                 error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                 var.constraint=var.constraint,maximum=TRUE))
            if(method[1]=="HL")
            {
              run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                      error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                      var.constraint=var.constraint))
            }
            if(method[1]=="HL-TS")
            {
              run.time2<-system.time(re<-cma.delta.lm(dat,delta=re1$maximum,max.itr=0,tol=tol,error.indep=error.indep,
                                                      error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                      var.constraint=var.constraint))
            }
            
            run.time<-run.time1+run.time2
          }
        }else
        {
          if(method[1]=="TS")
          {
            run.time<-system.time(re<-cma.delta.lm(dat,delta=delta,max.itr=0,tol=tol,error.indep=error.indep,
                                                   error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                   var.constraint=var.constraint))
          }
          if(method[1]=="HL")
          {
            run.time<-system.time(re<-cma.delta.lm(dat,delta=delta,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                   error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                   var.constraint=var.constraint))
          }
        }
      }
  
  re$time<-run.time
  
  return(re)
  if(as.numeric(re$Var.comp[1])<1e-5)
  {
    warning("The variance of A's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[2])<1e-5)
  {
    warning("The variance of C's random effect is less than 1e-5.")
  }
  return(re)
  if(as.numeric(re$Var.comp[3])<1e-5)
  {
    warning("The variance of B's random effect is less than 1e-5.")
  }
}

### ---- from macc/macc/R/sim.data.multi.R ----
sim.data.multi <-
function(Z.list,N,K=1,Theta,Sigma,Psi=diag(rep(1,3)),Lambda=diag(rep(1,3)))
{
  n<-matrix(NA,N,K)
  for(i in 1:N)
  {
    for(j in 1:K)
    {
      n[i,j]<-length(Z.list[[i]][[j]])
    }
  }
  
  if(K>1)
  {
    s.Psi<-svd(Psi)
    Psi.root<-s.Psi$u%*%diag(sqrt(s.Psi$d))%*%t(s.Psi$v)
    
    s.Lambda<-svd(Lambda)
    Lambda.root<-s.Lambda$u%*%diag(sqrt(s.Lambda$d))%*%t(s.Lambda$v)
    
    u<-matrix(rnorm(3*N),nrow=N)%*%Psi.root
    alpha<-u[,1]
    beta<-u[,2]
    gamma<-u[,3]
    
    A=B=C<-matrix(NA,N,K)
    epsA=epsB=epsC<-matrix(NA,N,K)
    for(i in 1:N)
    {
      eta<-matrix(rnorm(3*K),nrow=K)%*%Lambda.root
      
      epsA[i,]<-eta[,1]
      epsB[i,]<-eta[,2]
      epsC[i,]<-eta[,3]
    }
    A<-Theta[1,1]+matrix(rep(alpha,K),nrow=N)+epsA
    B<-Theta[2,2]+matrix(rep(beta,K),nrow=N)+epsB
    C<-Theta[1,2]+matrix(rep(gamma,K),nrow=N)+epsC
    
    dat<-list()
    for(i in 1:N)
    {
      dat[[i]]<-list()
      for(j in 1:K)
      {
        dat[[i]][[j]]<-sim.data.single(Z.list[[i]][[j]],Theta=matrix(c(A[i,j],0,C[i,j],B[i,j]),2,2),Sigma)
      }
    }
    
    re<-list(data=dat,A=A,B=B,C=C,type="multilevel")
  }else
  {
    s.Lambda<-svd(Lambda)
    Lambda.root<-s.Lambda$u%*%diag(sqrt(s.Lambda$d))%*%t(s.Lambda$v)
    
    eta<-matrix(rnorm(3*N),nrow=N)%*%Lambda.root
    A<-Theta[1,1]+eta[,1]
    B<-Theta[2,2]+eta[,2]
    C<-Theta[1,2]+eta[,3]
    
    dat<-list()
    for(i in 1:N)
    {
      dat[[i]]<-sim.data.single(Z.list[[i]],Theta=matrix(c(A[i],0,C[i],B[i]),2,2),Sigma)
    }
    
    re<-list(data=dat,A=A,B=B,C=C,type="twolevel")
  }
  
  return(re)
}

### ---- from macc/macc/R/sim.data.single.R ----
sim.data.single <-
function(Z,Theta,Sigma)
{
  n<-length(Z)
  
  s<-svd(Sigma)
  Sigma.root<-s$u%*%diag(sqrt(s$d))%*%t(s$v)
  
  E<-matrix(rnorm(2*n),nrow=n)%*%Sigma.root
  
  M<-Z*Theta[1,1]+E[,1]
  R<-Z*Theta[1,2]+M*Theta[2,2]+E[,2]
  
  re<-data.frame(Z=Z,M=M,R=R)
  return(re)
}

