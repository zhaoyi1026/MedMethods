##############################################################################
# MedMethods method module: gma
# Granger mediation analysis
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from gma/gma/R/cma.delta.ts.arp.error.lm.HL.R ----
cma.delta.ts.arp.error.lm.HL <-
function(dat,delta=0,p=1,max.itr=500,tol=1e-4,error.indep=FALSE,error.var.equal=FALSE,Sigma.update=FALSE,
                                       var.constraint=FALSE)
{
  re<-cma.delta.ts.arp.error.lm(dat,delta,p=p,max.itr,tol,error.indep,error.var.equal,Sigma.update,var.constraint)
  return(re$HL)
}

### ---- from gma/gma/R/cma.delta.ts.arp.error.lm.R ----
cma.delta.ts.arp.error.lm <-
function(dat,delta=0,p=1,max.itr=500,tol=1e-4,error.indep=FALSE,error.var.equal=FALSE,
                                    Sigma.update=FALSE,var.constraint=FALSE)
{
  N<-length(dat)
  K<-1
  
  ##################################################################
  # Estimate A, B, C and C' for each subject
  At=Bt=Ct=C2t<-matrix(NA,N,K)
  sigma1.hat=sigma2.hat<-matrix(NA,N,K)
  Wt<-array(NA,c(N,2*p,2))
  thetat<-matrix(NA,N,6*p+3)
  if(p==1)
  {
    colnames(thetat)<-c("A","phi1","psi11","psi21","C","phi2","psi12","psi22","B")
  }else
  {
    colnames(thetat)<-c("A",paste0(rep(c("phi1","psi11","psi21"),each=p),"_",rep(1:p,3)),
                        "C",paste0(rep(c("phi2","psi12","psi22"),each=p),"_",rep(1:p,3)),
                        "B")
  }
  n<-matrix(NA,N,K)
  for(i in 1:N)
  {
    dd<-dat[[i]]
    n[i,1]<-nrow(dd)
    
    re<-cma.uni.delta.ts.arp.error(dd,delta=delta,p=p,var.asmp=FALSE)
    
    At[i,1]<-re$Coefficients[1,1]
    Bt[i,1]<-re$Coefficients[3,1]
    Ct[i,1]<-re$Coefficients[2,1]
    C2t[i,1]<-re$Coefficients[4,1]
    Wt[i,,]<-c(re$W)
    
    Theta<-re$D
    
    phi<-matrix(NA,p,2)
    colnames(phi)<-c("phi1","phi2")
    phi[,1]<--Wt[i,seq(1,2*p-1,by=2),1]*Theta[1,1]-Wt[i,seq(2,2*p,by=2),1]*Theta[1,2]
    phi[,2]<--Wt[i,seq(1,2*p-1,by=2),2]*Theta[1,1]-Wt[i,seq(2,2*p,by=2),2]*Theta[1,2]
    
    Psi<-matrix(NA,2*p,2)
    Psi[seq(1,2*p-1,by=2),1]<-Wt[i,seq(1,2*p-1,by=2),1]-Wt[i,seq(2,2*p,by=2),1]*Theta[2,2]
    Psi[seq(2,2*p,by=2),1]<-Wt[i,seq(2,2*p,by=2),1]
    Psi[seq(1,2*p-1,by=2),2]<-Wt[i,seq(1,2*p-1,by=2),2]-Wt[i,seq(2,2*p,by=2),2]*Theta[2,2]
    Psi[seq(2,2*p,by=2),2]<-Wt[i,seq(2,2*p,by=2),2]
    
    theta1<-c(Theta[1,1],phi[,1],Psi[seq(1,2*p-1,by=2),1],Psi[seq(2,2*p,by=2),1])
    theta2<-c(Theta[1,2],phi[,2],Psi[seq(1,2*p-1,by=2),2],Psi[seq(2,2*p,by=2),2])
    
    thetat[i,]<-c(theta1,theta2,Bt[i,1])
    
    sigma1.hat[i,1]<-re$Sigma[1,1]
    sigma2.hat[i,1]<-re$Sigma[2,2]
  }
  ##################################################################
  
  ##################################################################
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
      Lambda.confint[1,]<-(intervals(fit.A)[[2]][c(1,3)])^2
      Lambda.confint[2,]<-(intervals(fit.B)[[2]][c(1,3)])^2
      Lambda.confint[3,]<-(intervals(fit.C)[[2]][c(1,3)])^2
    }
  }else
    if(var.constraint==TRUE)
    {
      Lambda.confint<-matrix(NA,3,2)
      colnames(Lambda.confint)<-c("LB","UB")
      rownames(Lambda.confint)<-c("A","B","C")
      Lambda.confint[1,]<-(intervals(fit.A)[[2]][c(1,3)])^2
      Lambda.confint[2,]<-(intervals(fit.B)[[2]][c(1,3)])^2
      Lambda.confint[3,]<-(intervals(fit.C)[[2]][c(1,3)])^2
    }else
    {
      Lambda.confint<-NULL
    }
  ##################################################################
  
  ##################################################################
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
    if(error.indep==TRUE)
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
  ##################################################################
  
  ##################################################################
  # max.itr=0: two-stage approach
  # max.itr>0: full likelihood (h-likelihood) approach
  J1<-cbind(diag(rep(1,3*p+1)),matrix(0,3*p+1,3*p+1),matrix(0,3*p+1,1))
  J2<-cbind(matrix(0,3*p+1,3*p+1),diag(rep(1,3*p+1)),matrix(0,3*p+1,1))
  J3<-matrix(c(rep(0,6*p+2),1),nrow=1)
  J<-rbind(c(1,rep(0,3*p+3*p+1+1)),c(rep(0,6*p+2),1),c(rep(0,3*p+1),1,rep(0,3*p+1)))
  
  diff<-100
  s<-0
  while(s<max.itr&diff>=tol)
  {
    bi.new<-matrix(NA,N,3)
    theta.new<-matrix(NA,N,6*p+3)
    for(i in 1:N)
    {
      dd<-dat[[i]]
      
      Zx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$Z[(p+1-x):(n[i,1]-x)])})
      Mx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$M[(p+1-x):(n[i,1]-x)])})
      Rx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$R[(p+1-x):(n[i,1]-x)])})
      
      Zy<-matrix(dd$Z[-(1:p)],ncol=1)
      My<-matrix(dd$M[-(1:p)],ncol=1)
      Ry<-matrix(dd$R[-(1:p)],ncol=1)
      
      X<-cbind(Zy,Zx,Mx,Rx)
      
      kappa<-delta*sqrt(sigma2.hat[i]/sigma1.hat[i])
      
      S1<-t(X%*%J1)%*%(X%*%J1)/sigma1.hat[i]+t(J)%*%solve(Lambda.hat)%*%J+
        t(My%*%J3+X%*%J2-kappa*X%*%J1)%*%(My%*%J3+X%*%J2-kappa*X%*%J1)/(sigma2.hat[i]*(1-delta^2))
      S2<-t(X%*%J1)%*%My/sigma1.hat[i]+t(My%*%J3+X%*%J2-kappa*X%*%J1)%*%(Ry-kappa*My)/(sigma2.hat[i]*(1-delta^2))+
        t(J)%*%solve(Lambda.hat)%*%b.hat
      theta.new[i,]<-solve(S1)%*%S2
      
      if(Sigma.update==TRUE)
      {
        theta1<-theta.new[i,1:(3*p+1)]
        theta2<-theta.new[i,(3*p+2):(6*p+2)]
        
        e1<-My-X%*%theta1
        e2<-Ry-My*theta.new[i,6*p+3]-X%*%theta2
        
        S<-matrix(NA,2,2)
        S[1,1]<-t(e1)%*%e1
        S[1,2]=S[2,1]<-t(e1)%*%e2
        S[2,2]<-t(e2)%*%e2
        
        sigma1.hat[i,1]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/((n[i,1]-p)*(1-delta^2))
        sigma2.hat[i,1]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/((n[i,1]-p)*(1-delta^2))
      }
      
      bi.new[i,]<-theta.new[i,c(1,6*p+3,3*p+2)]
    }
    b.new<-apply(bi.new,2,mean)
    Lambda.new<-matrix(0,3,3)
    Lambda.tmp<-t(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))%*%(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE))/N
    if(error.var.equal==FALSE)
    {
      if(error.indep==TRUE)
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
      if(error.indep==TRUE)
      {
        lambda2<-t(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))%*%(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))/(3*N)
        diag(Lambda.new)<-rep(lambda2,3)
      }else
      {
        warning("This variance structure is not valide! The errors are assumed to be independent.")
        lambda2<-t(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))%*%(c(bi.new-matrix(rep(b.new,N),ncol=3,byrow=TRUE)))/(3*N)
        diag(Lambda.new)<-rep(lambda2,3)
      }
    }
    
    diff<-max(abs(b.hat-b.new))
    
    b.hat<-b.new
    Lambda.hat<-Lambda.new
    At<-bi.new[,1]
    Bt<-bi.new[,2]
    Ct<-bi.new[,3]
    
    theta1.new<-theta.new[,1:(3*p+1)]
    theta2.new<-theta.new[,(3*p+2):(6*p+2)]
    
    Wt[,seq(2,2*p,by=2),1]<-theta1.new[,(2*p+2):(3*p+1)]
    Wt[,seq(2,2*p,by=2),2]<-theta2.new[,(2*p+2):(3*p+1)]
    Wt[,seq(1,2*p-1,by=2),1]<-theta1.new[,((p+2):(2*p+1))]+theta1.new[,(2*p+2):(3*p+1)]*Bt
    Wt[,seq(1,2*p-1,by=2),2]<-theta2.new[,((p+2):(2*p+1))]+theta2.new[,(2*p+2):(3*p+1)]*Bt
    
    s<-s+1
  }
  ##################################################################
  
  ##################################################################
  # summary results
  HL<-cma.ts.arp.error.lm.h(dat,delta=delta,Ai=At,Bi=Bt,Ci=Ct,Wi=Wt,b=b.hat,Lambda=Lambda.hat,p=p,Sigma.update=Sigma.update)
  if(max.itr==0)
  {
    # re.HL<-HL$h2
    re.HL<-as.numeric(logLik(fit.A))+as.numeric(logLik(fit.B))+as.numeric(logLik(fit.C))
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
  sigma.hat<-cbind(sigma1.hat,sigma2.hat)
  colnames(sigma.hat)<-c("E1","E2")
  ##################################################################
  
  re<-list(delta=delta,Coefficients=coe.re,Lambda=Lambda.hat,Sigma=sigma.hat,W=apply(Wt,c(2,3),mean,na.rm=TRUE),HL=re.HL,
           convergence=(s<max.itr|max.itr==0),var.constraint=Lambda.confint)
  
  return(re)
}

### ---- from gma/gma/R/cma.delta.ts.arp.error.lm.ts.logLik.R ----
cma.delta.ts.arp.error.lm.ts.logLik <-
function(dat,delta=0,p=p,error.indep=FALSE,error.var.equal=FALSE)
{
  re<-cma.delta.ts.arp.error.lm.ts(dat,delta=delta,p=p,error.indep=error.indep,error.var.equal=error.var.equal)
  return(re$logLik.lm)
}

### ---- from gma/gma/R/cma.delta.ts.arp.error.lm.ts.R ----
cma.delta.ts.arp.error.lm.ts <-
function(dat,delta=0,p=1,error.indep=FALSE,error.var.equal=FALSE)
{
  N<-length(dat)
  K<-1
  
  ##############################################
  # Estimate A, B, C and C' for each subject
  At=Bt=Ct=C2t<-matrix(NA,N,K)
  sigma1.hat=sigma2.hat<-matrix(NA,N,K)
  Wi<-array(NA,c(N,2*p,2))
  for(i in 1:N)
  {
    dd<-dat[[i]]
    re<-cma.uni.delta.ts.arp.error(dd,delta=delta,p=p,var.asmp=FALSE)
    
    At[i,1]<-re$Coefficients[1,1]
    Ct[i,1]<-re$Coefficients[2,1]
    Bt[i,1]<-re$Coefficients[3,1]
    C2t[i,1]<-re$Coefficients[4,1]
    
    sigma1.hat[i]<-re$Sigma[1,1]
    sigma2.hat[i]<-re$Sigma[2,2]
    
    Wi[i,,]<-re$W
  }
  Wt<-apply(Wi,c(2,3),mean,na.rm=TRUE)
  ##############################################
  
  ##############################################
  #
  Lambda.hat<-matrix(0,3,3)
  colnames(Lambda.hat)=rownames(Lambda.hat)<-c("A","B","C")
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
      
      ll<-as.numeric(logLik(fit.A)+logLik(fit.B)+logLik(fit.C))
    }else
    {
      fit.b<-lm(cbind(At,Bt,Ct)~1)
      aov.b<-Anova(fit.b)
      b.hat<-as.vector(coef(fit.b))
      Lambda.hat<-aov.b$SSPE/aov.b$error.df
      res<-cbind(At-b.hat[1],Bt-b.hat[1],Ct-b.hat[3])
      
      ll<-log(2*pi)*(-3*N/2)-log(det(Lambda.hat))/2-sum(diag(solve(Lambda.hat)%*%t(res)%*%res))/2
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
      
      ll<-as.numeric(logLik(fit))
    }else
    {
      warning("This variance structure is not valid! The errors are assumed to be independent.")
      bt<-rbind(At,Bt,Ct)
      group<-c(rep("A",length(At)),rep("B",length(Bt)),rep("C",length(Ct)))
      fit<-lm(bt~group)
      b.hat<-as.vector(by(bt[,1],group,mean,na.rm=TRUE))
      diag(Lambda.hat)<-rep((summary(fit)$sigma)^2,3)
      
      ll<-as.numeric(logLik(fit))
    }
  }
  ##############################################
  
  ##############################################
  # summary results
  AB.p<-b.hat[1]*b.hat[2]
  AB.d<-mean(C2t,na.rm=TRUE)-b.hat[3]
  coe.re<-matrix(NA,6,1)
  colnames(coe.re)<-c("Estimate")
  rownames(coe.re)<-c("A","C","B","C2","AB.prod","AB.diff")
  coe.re[,1]<-c(b.hat[1],b.hat[3],b.hat[2],mean(C2t,na.rm=TRUE),AB.p,AB.d)
  sigma.hat<-cbind(sigma1.hat,sigma2.hat)
  colnames(sigma.hat)<-c("E1","E2")
  ##############################################
  
  re<-list(delta=delta,Coefficients=coe.re,Lambda=Lambda.hat,Sigma=sigma.hat,W=Wt,logLik.lm=ll)
  return(re)
}

### ---- from gma/gma/R/cma.ts.arp.error.lm.h.R ----
cma.ts.arp.error.lm.h <-
function(dat,delta=0,Ai,Bi,Ci,Wi,b,Lambda,p=NULL,Sigma.update=FALSE)
{
  N<-length(dat)
  K<-1
  
  if(is.null(p))
  {
    p<-dim(Wi)[2]/2
  }
  
  sigma1.hat=sigma2.hat=n<-matrix(NA,N,K)
  h1=h2<-0
  for(i in 1:N)
  {
    dd<-dat[[i]]
    n[i,1]<-nrow(dd)
    
    ###################################################
    # h1
    Zx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$Z[(p+1-x):(n[i,1]-x)])})
    Mx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$M[(p+1-x):(n[i,1]-x)])})
    Rx<-apply(matrix(1:p,ncol=1),1,function(x){return(dd$R[(p+1-x):(n[i,1]-x)])})
    
    Zy<-matrix(dd$Z[-(1:p)],ncol=1)
    My<-matrix(dd$M[-(1:p)],ncol=1)
    Ry<-matrix(dd$R[-(1:p)],ncol=1)
    
    X<-cbind(Zy,Zx,Mx,Rx)
    
    Theta<-matrix(c(Ai[i],0,Ci[i],Bi[i]),2,2)
    
    phi<-matrix(NA,p,2)
    colnames(phi)<-c("phi1","phi2")
    phi[,1]<--Wi[i,seq(1,2*p-1,by=2),1]*Theta[1,1]-Wi[i,seq(2,2*p,by=2),1]*Theta[1,2]
    phi[,2]<--Wi[i,seq(1,2*p-1,by=2),2]*Theta[1,1]-Wi[i,seq(2,2*p,by=2),2]*Theta[1,2]
    
    Psi<-matrix(NA,2*p,2)
    Psi[seq(1,2*p-1,by=2),1]<-Wi[i,seq(1,2*p-1,by=2),1]-Wi[i,seq(2,2*p,by=2),1]*Theta[2,2]
    Psi[seq(2,2*p,by=2),1]<-Wi[i,seq(2,2*p,by=2),1]
    Psi[seq(1,2*p-1,by=2),2]<-Wi[i,seq(1,2*p-1,by=2),2]-Wi[i,seq(2,2*p,by=2),2]*Theta[2,2]
    Psi[seq(2,2*p,by=2),2]<-Wi[i,seq(2,2*p,by=2),2]
    
    theta1<-c(Theta[1,1],phi[,1],Psi[seq(1,2*p-1,by=2),1],Psi[seq(2,2*p,by=2),1])
    theta2<-c(Theta[1,2],phi[,2],Psi[seq(1,2*p-1,by=2),2],Psi[seq(2,2*p,by=2),2])
    
    if(Sigma.update)
    {
      e1<-My-X%*%theta1
      e2<-Ry-My*Bi[i]-X%*%theta2
      
      S<-matrix(NA,2,2)
      S[1,1]<-t(e1)%*%e1
      S[1,2]=S[2,1]<-t(e1)%*%e2
      S[2,2]<-t(e2)%*%e2
      
      sigma1.hat[i,1]<-(S[1,1]-delta*S[1,2]*sqrt(S[1,1]/S[2,2]))/((n[i,1]-p)*(1-delta^2))
      sigma2.hat[i,1]<-(S[2,2]-delta*S[1,2]*sqrt(S[2,2]/S[1,1]))/((n[i,1]-p)*(1-delta^2))
    }else
    {
      re<-cma.uni.delta.ts.arp.error(dd,delta,p=p,var.asmp=FALSE)
      sigma1.hat[i,1]<-re$Sigma[1,1]
      sigma2.hat[i,1]<-re$Sigma[2,2]
    }
    
    Sigma<-matrix(c(sigma1.hat[i,1],delta*sqrt(sigma1.hat[i,1]*sigma2.hat[i,1]),
                    delta*sqrt(sigma1.hat[i,1]*sigma2.hat[i,1]),sigma2.hat[i,1]),2,2)
    
    h1<-h1+cma.uni.ts.arp.error.ll(dd,Theta,W=matrix(Wi[i,,],ncol=2),Sigma=Sigma,p=p)
    
    b.hat<-c(Ai[i],Bi[i],Ci[i])
    
    h2<-h2-((log(det(Lambda))+t(b.hat-b)%*%solve(Lambda)%*%(b.hat-b))/2)[1,1]
  }
  
  const1<--log(2*pi)*sum(n-p)
  const2<--log(2*pi)*N*3/2
  
  h<-h1+h2
  re<-data.frame(h1=const1+h1,h2=const2+h2,h=const1+const2+h)
  
  return(re)
}

### ---- from gma/gma/R/cma.uni.delta.ts.arp.error.R ----
cma.uni.delta.ts.arp.error <-
function(dat,delta=0,p=1,conf.level=0.95,var.asmp=TRUE)
{
  Z<-matrix(dat$Z,ncol=1)
  M<-matrix(dat$M,ncol=1)
  R<-matrix(dat$R,ncol=1)
  
  n<-nrow(Z)
  
  z.alpha<-qnorm(1-(1-conf.level)/2)
  
  # time series: AR(p) autoregressive errors
  Zx<-apply(matrix(1:p,ncol=1),1,function(x){return(Z[(p+1-x):(n-x)])})
  Mx<-apply(matrix(1:p,ncol=1),1,function(x){return(M[(p+1-x):(n-x)])})
  Rx<-apply(matrix(1:p,ncol=1),1,function(x){return(R[(p+1-x):(n-x)])})
  
  Zy<-matrix(Z[-(1:p),1],ncol=1)
  My<-matrix(M[-(1:p),1],ncol=1)
  Ry<-matrix(R[-(1:p),1],ncol=1)
  
  ###########################################################
  # total effect model
  X<-cbind(Zy,Zx,Mx,Rx)
  fit.M<-lm(My~0+X)
  fit.R2<-lm(Ry~0+X)
  
  beta.M<-coef(fit.M)
  beta.R2<-coef(fit.R2)
  
  Sigma.B.hat<-(t(cbind(My-X%*%beta.M,Ry-X%*%beta.R2))%*%(cbind(My-X%*%beta.M,Ry-X%*%beta.R2)))/(n-1)
  # sigma12.hat<-Sigma.B.hat[1,1]
  # sigma22.hat<-det(Sigma.B.hat)/(Sigma.B.hat[1,1]*(1-delta^2))
  # Sigma.hat<-matrix(c(sigma12.hat,delta*sqrt(sigma12.hat*sigma22.hat),delta*sqrt(sigma12.hat*sigma22.hat),sigma22.hat),2,2)
  
  # Variance of beta.M and beta.R2
  beta.M.var<-Sigma.B.hat[1,1]*solve(t(X)%*%X)
  beta.R2.var<-Sigma.B.hat[2,2]*solve(t(X)%*%X)
  ###########################################################
  
  ###########################################################
  # coefficient estimate given delta
  
  # projection matrices
  P.M<-My%*%solve(t(My)%*%My)%*%t(My)
  O.M<-diag(rep(1,n-p))-P.M
  P.MX<-O.M%*%X%*%solve(t(X)%*%O.M%*%X)%*%t(X)%*%O.M
  P.X<-X%*%solve(t(X)%*%X)%*%t(X)
  O.X<-diag(rep(1,n-p))-P.X
  
  sigma12.hat<-(t(My)%*%O.X%*%My/(n-p))[1,1]
  sigma22.hat<-(t(Ry)%*%(O.M-P.MX)%*%Ry/((n-p)*(1-delta^2)))[1,1]
  Sigma.hat<-matrix(c(sigma12.hat,delta*sqrt(sigma12.hat*sigma22.hat),delta*sqrt(sigma12.hat*sigma22.hat),sigma22.hat),2,2)
  
  kappa.hat<-delta*sqrt(sigma22.hat/sigma12.hat)
  
  theta1.hat<-solve(t(X)%*%X)%*%t(X)%*%My
  theta2.hat<-solve(t(X)%*%O.M%*%X)%*%t(X)%*%O.M%*%Ry+kappa.hat*theta1.hat
  B.hat<-(solve(t(My)%*%My)%*%t(My)%*%(diag(rep(1,n-p))-X%*%solve(t(X)%*%O.M%*%X)%*%t(X)%*%O.M)%*%Ry-kappa.hat)[1,1]
  
  # C'
  C2.hat<-coef(fit.R2)[1]
  
  # A, B and C
  A.hat<-theta1.hat[1,1]
  C.hat<-theta2.hat[1,1]
  ABp.hat<-A.hat*B.hat
  ABd.hat<-C2.hat-C.hat
  
  # transition matrix of error
  W.hat<-matrix(NA,2*p,2)
  W.hat[seq(2,2*p,by=2),c(1,2)]<-cbind(theta1.hat[(2*p+2):(3*p+1)],theta2.hat[(2*p+2):(3*p+1)])
  W.hat[seq(1,2*p-1,by=2),1]<-theta1.hat[(p+2):(2*p+1)]+theta1.hat[(2*p+2):(3*p+1)]*B.hat
  W.hat[seq(1,2*p-1,by=2),2]<-theta2.hat[(p+2):(2*p+1)]+theta2.hat[(2*p+2):(3*p+1)]*B.hat
  
  phi.hat<-cbind(theta1.hat[2:(p+1)],theta2.hat[2:(p+1)])
  colnames(phi.hat)<-c("phi1","phi2")
  
  if(var.asmp)
  {
    J1<-kronecker(diag(rep(1,p)),c(1,0))
    J2<-kronecker(diag(rep(1,p)),c(0,1))
    if(p==1)
    {
      Fm<-t(W.hat)
    }else
    {
      Fm<-rbind(t(W.hat),cbind(diag(rep(1,2*(p-1))),matrix(0,2*(p-1),2))) 
    }
    Xi<-matrix(0,2*p,2*p)
    Xi[1:2,1:2]<-Sigma.hat
    Pi<-matrix(solve(diag(rep(1,(2*p)^2))-kronecker(Fm,Fm))%*%c(Xi),2*p,2*p)
    
    # Xt*Xt
    q<-mean(Z)
    ZZ<-q*diag(rep(1,p))
    ZM<-A.hat*q*diag(rep(1,p))
    ZR<-(C.hat+ABp.hat)*q*diag(rep(1,p))
    MM<-A.hat^2*q*diag(rep(1,p))+t(J1)%*%Pi%*%J1
    MR<-A.hat*(C.hat+ABp.hat)*q*diag(rep(1,p))+B.hat*t(J1)%*%Pi%*%J1+t(J1)%*%Pi%*%J2
    RR<-(C.hat+ABp.hat)^2*q*diag(rep(1,p))+B.hat^2*t(J1)%*%Pi%*%J1+B.hat*t(J1)%*%Pi%*%J2+B.hat*t(J2)%*%Pi%*%J1+t(J1)%*%Pi%*%J1
    
    X.cov<-matrix(NA,3*p+1,3*p+1)
    X.cov[1,]<-c(q,rep(0,p),rep(0,p),rep(0,p))
    X.cov[2:(p+1),]<-cbind(rep(0,p),ZZ,ZM,ZR)
    X.cov[(p+2):(2*p+1),]<-cbind(rep(0,p),t(ZM),MM,MR)
    X.cov[(2*p+2):(3*p+1),]<-cbind(rep(0,p),t(ZR),t(MR),RR)
    
    # Xt*Mt
    psi11.hat<-matrix(W.hat[seq(1,2*p,by=2),1]-W.hat[seq(2,2*p,by=2),1]*B.hat,ncol=1)
    psi21.hat<-matrix(W.hat[seq(2,2*p,by=2),1],ncol=1)
    XM.cov<-matrix(NA,3*p+1,1)
    XM.cov[1,1]<-A.hat*q
    XM.cov[2:(p+1),1]<-ZZ%*%phi.hat[,1]+ZM%*%psi11.hat+ZR%*%psi21.hat
    XM.cov[(p+2):(2*p+1),1]<-t(ZM)%*%phi.hat[,1]+MM%*%psi11.hat+MR%*%psi21.hat
    XM.cov[(2*p+2):(3*p+1),1]<-t(ZR)%*%phi.hat[,1]+t(MR)%*%psi11.hat+RR%*%psi21.hat
    
    M.cov<-A.hat^2*q+Sigma.hat[1,1]+t(phi.hat[,1])%*%(ZZ%*%phi.hat[,1]+ZM%*%psi11.hat+ZR%*%psi21.hat)+
      t(psi11.hat)%*%(t(ZM)%*%phi.hat[,1]+MM%*%psi11.hat+MR%*%psi21.hat)+
      t(psi21.hat)%*%(t(ZR)%*%phi.hat[,1]+t(MR)%*%psi11.hat+RR%*%psi21.hat)
    
    dn1<-sigma12.hat*(1-delta^2)
    dn2<-sigma22.hat*(1-delta^2)
    Fisher.info<-rbind(cbind(X.cov/dn1,-kappa.hat*X.cov/dn2,-kappa.hat*XM.cov/dn2),
                       cbind(-kappa.hat*t(X.cov)/dn2,X.cov/dn2,XM.cov/dn2),
                       cbind(-kappa.hat*t(XM.cov)/dn2,t(XM.cov)/dn2,M.cov/dn2))*(n-p)
    
    # theta.var<-solve(Fisher.info)
    theta.var<-ginv(Fisher.info)
    
    theta1.var<-theta.var[1:(3*p+1),1:(3*p+1)]
    theta2.var<-theta.var[(3*p+2):(6*p+2),(3*p+2):(6*p+2)]
    B.var<-theta.var[6*p+3,6*p+3]
  }else
  {
    theta1.var<-solve(t(X)%*%X/sigma12.hat)
    theta2.var<-solve(t(X)%*%X/(sigma22.hat*(1-delta^2)))
    B.var<-((sigma22.hat*(1-delta^2))/t(M)%*%M)[1,1] 
  }
  
  C2.hat.se<-sqrt(beta.R2.var[1,1])
  
  A.hat.se<-sqrt(theta1.var[1,1])
  B.hat.se<-sqrt(B.var)
  C.hat.se<-sqrt(theta2.var[1,1])
  
  ABp.hat.se<-sqrt(A.hat^2*B.hat.se^2+B.hat^2*A.hat.se^2)
  ABd.hat.se<-sqrt(C2.hat.se^2+C.hat.se^2)
  
  cma.re<-matrix(NA,6,4)
  rownames(cma.re)<-c("A","C","B","C2","AB.p","AB.d")
  colnames(cma.re)<-c("Estimate","SE","LB","UB")
  cma.re[,1]<-c(A.hat,C.hat,B.hat,C2.hat,ABp.hat,ABd.hat)
  cma.re[,2]<-c(A.hat.se,C.hat.se,B.hat.se,C2.hat.se,ABp.hat.se,ABd.hat.se)
  cma.re[,3]<-cma.re[,1]-z.alpha*cma.re[,2]
  cma.re[,4]<-cma.re[,1]+z.alpha*cma.re[,2]
  D.hat<-matrix(c(A.hat,0,C.hat,B.hat),2,2)
  
  re<-list(Coefficients=cma.re,D=D.hat,Sigma=Sigma.hat,delta=delta,W=W.hat)
  
  return(re)
}

### ---- from gma/gma/R/cma.uni.plot.ts.arp.error.R ----
cma.uni.plot.ts.arp.error <-
function(re.cma.sens,re.cma=NULL,delta=NULL,legend.pos="topright",
                                    xlab=expression(delta),ylab=expression(hat(AB)),
                                    cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,...)
{
  dt<-re.cma.sens$coefficients[,"delta"]
  AB.p<-re.cma.sens$coefficients[,"AB.p.Estimate"]
  AB.p.ub<-re.cma.sens$coefficients[,"AB.p.UB"]
  AB.p.lb<-re.cma.sens$coefficients[,"AB.p.LB"]
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
  lines(dt,AB.p.lb,lty=2,lwd=1,col=8)
  lines(dt,AB.p.ub,lty=2,lwd=1,col=8)
  
  if(!is.null(re.cma))
  {
    points(re.cma$delta,re.cma$Coefficients[5,1],pch=16,col=4,cex=0.75)
    lines(rep(re.cma$delta,2),re.cma$Coefficients[5,c(3,4)],lty=2,col=4)
    lines(c(re.cma$delta-0.02,re.cma$delta+0.02),rep(re.cma$Coefficients[5,3],2),col=4)
    lines(c(re.cma$delta-0.02,re.cma$delta+0.02),rep(re.cma$Coefficients[5,4],2),col=4)
    
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
  }else
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

### ---- from gma/gma/R/cma.uni.sens.ts.arp.error.R ----
cma.uni.sens.ts.arp.error <-
function(dat,delta=seq(-1,1,by=0.01),p=1,conf.level=0.95,var.asmp=TRUE)
{
  Coe=Sp=dt=Wmat<-NULL
  for(i in 1:length(delta))
  {
    if(abs(delta[i])!=1)
    {
      # delta
      dt<-c(dt,delta[i])
      
      re<-cma.uni.delta.ts.arp.error(dat,delta=delta[i],p=p,conf.level=conf.level,var.asmp=var.asmp)
      
      # Coefficients and confidence intervals
      Coe<-rbind(Coe,c(t(re$Coefficients)))
      colnames(Coe)<-paste(rep(rownames(re$Coefficients),each=ncol(re$Coefficients)),
                           rep(colnames(re$Coefficients),nrow(re$Coefficients)),sep=".")
      
      # Sigma
      Sp<-rbind(Sp,c(re$Sigma[1,1],re$Sigma[2,2],re$Sigma[1,2]))
      colnames(Sp)<-c("sigma12","sigma22","rho")
      
      # W matrix
      Wmat<-rbind(Wmat,c(re$W))
      colnames(Wmat)<-paste0(rep(c("W11","W21","W12","W22"),p),"_",rep(1:p,each=4))
    }
  }
  return(list(coefficients=cbind(delta=dt,Coe),Sigma=cbind(Sp,delta=dt),W=Wmat))
}

### ---- from gma/gma/R/cma.uni.ts.arp.error.ll.R ----
cma.uni.ts.arp.error.ll <-
function(dat,Theta,W,Sigma,p=NULL)
{
  Z<-matrix(dat$Z,ncol=1)
  M<-matrix(dat$M,ncol=1)
  R<-matrix(dat$R,ncol=1)
  
  n<-nrow(Z)
  
  if(is.null(p))
  {
    p<-nrow(W)/2
  }
  
  # time series: AR(p) autoregressive errors
  Zx<-apply(matrix(1:p,ncol=1),1,function(x){return(Z[(p+1-x):(n-x)])})
  Mx<-apply(matrix(1:p,ncol=1),1,function(x){return(M[(p+1-x):(n-x)])})
  Rx<-apply(matrix(1:p,ncol=1),1,function(x){return(R[(p+1-x):(n-x)])})
  
  Zy<-matrix(Z[-(1:p),1],ncol=1)
  My<-matrix(M[-(1:p),1],ncol=1)
  Ry<-matrix(R[-(1:p),1],ncol=1)
  
  X<-cbind(Zy,Zx,Mx,Rx)
  
  phi<-matrix(NA,p,2)
  colnames(phi)<-c("phi1","phi2")
  phi[,1]<--W[seq(1,2*p-1,by=2),1]*Theta[1,1]-W[seq(2,2*p,by=2),1]*Theta[1,2]
  phi[,2]<--W[seq(1,2*p-1,by=2),2]*Theta[1,1]-W[seq(2,2*p,by=2),2]*Theta[1,2]
  
  Psi<-matrix(NA,2*p,2)
  Psi[seq(1,2*p-1,by=2),1]<-W[seq(1,2*p-1,by=2),1]-W[seq(2,2*p,by=2),1]*Theta[2,2]
  Psi[seq(2,2*p,by=2),1]<-W[seq(2,2*p,by=2),1]
  Psi[seq(1,2*p-1,by=2),2]<-W[seq(1,2*p-1,by=2),2]-W[seq(2,2*p,by=2),2]*Theta[2,2]
  Psi[seq(2,2*p,by=2),2]<-W[seq(2,2*p,by=2),2]
  
  theta1<-c(Theta[1,1],phi[,1],Psi[seq(1,2*p-1,by=2),1],Psi[seq(2,2*p,by=2),1])
  theta2<-c(Theta[1,2],phi[,2],Psi[seq(1,2*p-1,by=2),2],Psi[seq(2,2*p,by=2),2])
  
  sigma12<-Sigma[1,1]
  sigma22<-Sigma[2,2]
  delta<-Sigma[1,2]/sqrt(sigma12*sigma22)
  
  kappa<-delta*sqrt(sigma22/sigma12)
  
  const<--(n-p)*log(sigma12*sigma22*(1-delta^2))/2
  
  ll1<-(-t(My-X%*%theta1)%*%(My-X%*%theta1)/(2*sigma12))[1,1]
  e2<-(Ry-My*Theta[2,2]-X%*%theta2)-kappa*(My-X%*%theta1)
  ll2<-(-t(e2)%*%e2/(2*sigma22*(1-delta^2)))[1,1]
  
  ll<-const+ll1+ll2
  return(ll)
}

### ---- from gma/gma/R/gma-internal.R ----
.Random.seed <-
c(403L, 408L, -1017079477L, 2059942818L, 533275634L, -2135981122L, 
-318880162L, 1906218233L, 271789403L, 1636142624L, -108752961L, 
-826386172L, -807017685L, 1774965448L, 1734265075L, -1871624294L, 
1487897610L, 937052826L, 853921576L, -341460200L, -1164634601L, 
1437304157L, 1612852952L, -1162997786L, -1563473304L, 1098869603L, 
-1048559943L, -938501410L, 1182846804L, 1770204495L, 1454580886L, 
-2120699501L, -1278269938L, 1273989249L, -538092211L, 1819737079L, 
1147527002L, -897038234L, -1015689789L, -57244779L, -894114801L, 
-404953472L, -1373714892L, 808875994L, -2064669204L, -209884387L, 
-245530530L, -756158231L, 1562228800L, -149430823L, 1098194970L, 
1085959493L, -294255444L, -153308469L, 827535419L, 1148711202L, 
-571803302L, 477913361L, -1543167024L, 1153484607L, 1970954655L, 
-580842638L, 858207747L, 863875697L, -200263715L, -1181213131L, 
323365123L, 1445406838L, -2097916381L, -1733319404L, 581836943L, 
-913485552L, 31025823L, -958177018L, 1700028632L, 1339298163L, 
-1867960034L, -1788139154L, -1261850927L, 1965138768L, -544903283L, 
2072782658L, 406051844L, 455257072L, -1717770740L, -405083L, 
132748353L, 506714441L, 1780698757L, -685043775L, -310248357L, 
445303664L, 959347624L, -1864840995L, -2099533385L, -876893543L, 
-1091360572L, 1123973167L, -754291057L, -1635871145L, -1675748358L, 
-729293502L, -1676018925L, 2144694828L, 829358386L, -166976249L, 
499033557L, 1378215453L, 1813921370L, -1771240758L, -1025848713L, 
-153398465L, -1090948886L, 1050236442L, 1077579284L, -1861744247L, 
-1176453328L, 835357989L, 1016437333L, -638990701L, -1495329005L, 
849595291L, -1123948526L, -981997248L, -1333985399L, 161692570L, 
-81881387L, -1120180817L, 840706179L, 194784483L, -1596915350L, 
-386130535L, -1202335441L, 772342578L, -1534093105L, -743382007L, 
-1958377597L, 4364912L, 1398397266L, 432647668L, -1262073806L, 
-493079058L, 747473664L, -1216991566L, 1505056438L, 2070385086L, 
1059527144L, -1320072963L, 932447762L, -156623981L, -1914314700L, 
-976424967L, -767494007L, -1759059655L, -1386913511L, -266425450L, 
1427377358L, -732671329L, -29531580L, -171214501L, -1596403433L, 
-1494332055L, -53791911L, -1507920074L, -2041319661L, -1854329430L, 
1719885959L, 785833685L, 541724291L, -829253750L, 1098130114L, 
755565903L, -98488616L, -597015824L, -1620038756L, -248182914L, 
1885975940L, -210608260L, -445718021L, -1546598369L, -1824077503L, 
918872853L, -1207123863L, 924941772L, -2021168544L, -871486084L, 
-681808504L, -1809629002L, -533664670L, 302344495L, 666281981L, 
-1870811661L, -562283265L, -300711005L, -1629444594L, -1334956882L, 
-361636613L, 1235733588L, -883080460L, -1277828362L, -1329089294L, 
-2103117458L, 1473748905L, 144467920L, -276167268L, 74325391L, 
721653352L, 1202558901L, -783625197L, -82191971L, -1896828453L, 
-1159274776L, 1745106380L, -1183041702L, -450032967L, -1318316184L, 
-402096778L, -821140931L, 1305538207L, -1250998724L, 877111154L, 
-1027984612L, 1463140019L, 657566642L, 1788068480L, -1142436984L, 
1901888706L, 564155723L, 1433233362L, 1177856611L, -1788623908L, 
-176254057L, -1734979474L, 110873241L, -1036603317L, -932906763L, 
1128836612L, -2058543313L, -1295973401L, -625012338L, 2054067452L, 
-1219463303L, 1865506420L, 1391644448L, 1673091543L, 856730154L, 
-1051763778L, 1009732001L, 297893328L, 491896911L, -1600640857L, 
-1041892923L, -136769404L, 284089492L, 1535425356L, 406547049L, 
940287435L, 2121376148L, 742142611L, 1772440345L, 2026419890L, 
-2079122036L, 1924240492L, -1960442607L, -183717876L, -2083531267L, 
225624689L, -1391430099L, 1213371670L, -1367883621L, 1466452228L, 
-210421453L, 1813282217L, 267089364L, -1614910193L, 746211287L, 
-2049256020L, 1952977042L, -702914099L, -365946071L, 233095543L, 
-1940919532L, 1247201687L, -1449719773L, 2086403507L, 1084830660L, 
2102865939L, -939542870L, -1234587943L, 1542526986L, -1930877359L, 
1099405557L, 280643327L, 1269410129L, -174837931L, 748003749L, 
1425540754L, 247854185L, -1644019975L, 579778824L, 468461862L, 
1483074597L, -747796728L, -2137503914L, 4581264L, 785056149L, 
-980612251L, -1082889349L, 1382117141L, 1627155149L, 1695439354L, 
2066465610L, -1391338269L, -961705871L, 1318138535L, 1853871539L, 
1473761223L, 932737309L, -910057009L, -904314107L, -1657290691L, 
13992709L, 1335134781L, -2033392540L, -1820443074L, 1777089993L, 
-962085577L, -154527777L, 883262922L, -621115638L, 1848891308L, 
-675602889L, -1482463844L, -143916843L, -104305480L, -1733326526L, 
-154211342L, -1505734487L, 1727618035L, 1544696703L, 602217013L, 
1730434982L, -1524744544L, -2111702563L, 1157086458L, 961798664L, 
791536410L, -631657849L, 1572412767L, 169162964L, -879681142L, 
781007696L, -1153809667L, -1701118724L, 332084190L, 1782153808L, 
-2002247041L, 1421015915L, 1761933056L, -1073503749L, 1580483678L, 
-699000994L, 670567730L, -2009391557L, -449674447L, 754573981L, 
615935929L, -1633986915L, -584996183L, 878046137L, -1418590024L, 
1624628171L, -996052062L, 1974683007L, -694175399L, -520702295L, 
-76873044L, 1863821843L, 2064454218L, 774363629L, 1254734102L, 
1998491029L, -197759242L, -721873827L, -1207718721L, -1515742375L, 
-1158557621L, -1870459589L, 517257934L, 1202478796L, 679272322L, 
-821249278L, -1278459682L, 1964625487L, 773239161L, 791359603L, 
1556649675L, 561385398L, 105076413L, 1732500437L, 844855726L, 
-1189095688L, -1762848700L, -1336489943L, 919055440L, -137130941L, 
-1144769707L, -266954569L, -681169087L, -1886466401L, -1343412835L, 
1890125515L, 1004628671L, -247044407L, -1105948527L, 899197336L, 
2040365713L, -1968680966L, -1957897965L, -978963642L, -930115235L, 
1915576134L, 1579450025L, 1292569617L, 776573421L, 1045595315L, 
1272154677L, -1330171592L, -783820126L, 443786446L, 1393895400L, 
-1405467447L, 2061478185L, 1618712447L, -1326517162L, 660208674L, 
1626406893L, 205250037L, -638480088L, -1896074497L, -1590269174L, 
52111823L, -614772845L, 355701291L, 1979570122L, -1513758332L, 
-1048020580L, -804576137L, 469281934L, -1793663007L, -1367938771L, 
1917101349L, -914580861L, 2008033583L, -1057439616L, 817666561L, 
-1322627152L, 1621511636L, -1379279580L, -1583132926L, 489523281L, 
-413012889L, -1015258743L, 564864981L, -2001421967L, -1801217327L, 
1061168367L, -1696982072L, -1188800667L, -1972987327L, -1985712253L, 
-373584868L, -1182525949L, -423559010L, 848752323L, -1273163849L, 
-2061957511L, -1504276063L, -1909712383L, -370632754L, -1342657904L, 
-958126926L, -1712565148L, -645017444L, 606763583L, -1702232653L, 
747964293L, -687168448L, -2109833589L, 2056927900L, -1504238682L, 
935651071L, -232356112L, 1136757329L, 1660092650L, 1441360769L, 
1276575166L, 1107301108L, 949680289L, -1542606107L, -1004659032L, 
2109615975L, 1393096519L, 195189128L, -207236302L, -434097359L, 
-1736499824L, -284947286L, -17568617L, -1257604226L, 471423183L, 
-114483957L, 1715284503L, -2120443201L, 2129476694L, -1410562278L, 
-703653630L, 497142094L, 1966127944L, -61876254L, 663113674L, 
-1725321614L, -1543949270L, -1098570836L, -554088780L, -839747798L, 
843238566L, 1792787013L, -915346259L, 1311963051L, -1587250900L, 
283048806L, -416286850L, -23082795L, 255767347L, 639685700L, 
823206846L, -693771016L, -382709675L, 134887853L, -2124607164L, 
1713693765L, -119739154L, -1337733941L, -1134169686L, -2114096532L, 
1337457929L, 553887014L, -1426097671L, -1627492331L, -1504354879L, 
-898719362L, -325441193L, -760891074L, 92987228L, -397987385L, 
-689592649L, 439940408L, -694833608L, -689001447L, -1264280279L, 
2043400604L, -567184355L, -1143712453L, -2137670650L, -1691724271L, 
-271058973L, -227326517L, 1098243593L, -643444494L, 1167917292L, 
-411716532L, -1117561811L, -1094905992L, -1962092056L, 1975184601L, 
522478848L, -833937839L, -1617105033L, -1107098430L, 1585340281L, 
-1269929762L, -74012726L, 729658649L, -1570826429L, 1163029581L, 
860360436L, -194455404L, 1322827579L, 212266455L, 1821674510L, 
1739382562L, 1512223535L, -805547110L, 1042387178L, 1406265206L, 
-1904447900L, -1086851992L, 1870820117L, 1451073841L, 1786900696L, 
532021278L, -1705195574L, 1866801267L, 1840294502L, -122242999L, 
1524002159L, -1585969429L, -1622887518L, 188527134L, -1273770711L, 
290046251L, 1391539505L, -296031941L, 605327891L, 1862598290L, 
995163101L, -260229466L, -2100358436L, -1322569936L, 158198560L, 
1057499036L, -1808541013L, -544963810L, 1046683335L, -1090803001L
)

### ---- from gma/gma/R/gma.R ----
gma <-
function(dat,model.type=c("single","twolevel"),method=c("HL","TS","HL-TS"),delta=NULL,p=1,
              single.var.asmp=TRUE,sens.plot=FALSE,sens.delta=seq(-1,1,by=0.01),legend.pos="topright",
              xlab=expression(delta),ylab=expression(hat(AB)),
              cex.lab=1,cex.axis=1,lgd.cex=1,lgd.pt.cex=1,plot.delta0=TRUE,
              interval=c(-0.9,0.9),tol=1e-4,max.itr=500,conf.level=0.95,error.indep=TRUE,error.var.equal=FALSE,
              Sigma.update=TRUE,var.constraint=TRUE,...)
{
  if(model.type[1]=="single")
  {
    if(is.null(delta)==TRUE)
    {
      delta<-0
    }
    run.time<-system.time(re1<-cma.uni.delta.ts.arp.error(dat,delta=delta,p=p,conf.level=conf.level,var.asmp=single.var.asmp))
    
    single.ll<-cma.uni.ts.arp.error.ll(dat,Theta=re1$D,W=re1$W,Sigma=re1$Sigma,p=p)
    
    re<-re1
    re$LL<-single.ll
    
    # sensitivity plot
    if(sens.plot==TRUE)
    {
      re.cma.sens<-cma.uni.sens.ts.arp.error(dat,delta=sens.delta,p=p,conf.level=conf.level,var.asmp=single.var.asmp)
      cma.uni.plot.ts.arp.error(re.cma.sens,re.cma=re1,delta=delta,legend.pos=legend.pos,
                                xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis,
                                lgd.cex=lgd.cex,lgd.pt.cex=lgd.pt.cex,plot.delta0=plot.delta0)
    }
  }else
    if(model.type[1]=="twolevel")
    {
      if(is.null(delta)==TRUE)
      {
        if(method[1]=="TS")
        {
          t1<-system.time(re1<-optimize(cma.delta.ts.arp.error.lm.HL,interval=interval,dat=dat,p=p,max.itr=0,tol=tol,
                                        error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                        var.constraint=var.constraint,maximum=TRUE))
          t2<-system.time(re<-cma.delta.ts.arp.error.lm(dat,delta=re1$maximum,p=p,max.itr=0,tol=tol,error.indep=error.indep,
                                                        error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                        var.constraint=var.constraint))
          
          run.time<-t1+t2
        }else
        {
          t1<-system.time(re1<-optimize(cma.delta.ts.arp.error.lm.HL,interval=interval,dat=dat,p=p,max.itr=max.itr,tol=tol,
                                        error.indep=error.indep,error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                        var.constraint=var.constraint,maximum=TRUE))
          if(method[1]=="HL")
          {
            t2<-system.time(re<-cma.delta.ts.arp.error.lm(dat,delta=re1$maximum,p=p,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                          error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                          var.constraint=var.constraint))
          }
          if(method[1]=="HL-TS")
          {
            t2<-system.time(re<-cma.delta.ts.arp.error.lm(dat,delta=re1$maximum,p=p,max.itr=0,tol=tol,error.indep=error.indep,
                                                          error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                          var.constraint=var.constraint))
          }
          
          run.time<-t1+t2
        }
      }else
      {
        if(method[1]=="TS")
        {
          run.time<-system.time(re<-cma.delta.ts.arp.error.lm(dat,delta=delta,p=p,max.itr=0,tol=tol,error.indep=error.indep,
                                                              error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                              var.constraint=var.constraint))
        }
        if(method[1]=="HL")
        {
          run.time<-system.time(re<-cma.delta.ts.arp.error.lm(dat,delta=delta,p=p,max.itr=max.itr,tol=tol,error.indep=error.indep,
                                                              error.var.equal=error.var.equal,Sigma.update=Sigma.update,
                                                              var.constraint=var.constraint))
        }
      }
    }
  
  re$time<-run.time
  
  return(re)
}

### ---- from gma/gma/R/sim.data.ts.single.R ----
sim.data.ts.single <-
function(n,Z,A,B,C,Sigma,W,Delta=NULL,p=NULL,nburn=100)
{
  if(is.null(Delta)==TRUE)
  {
    Delta<-Sigma
  }
  
  if(is.null(p)==TRUE)
  {
    p<-nrow(W)/2
  }
  
  # n subjects and nburn burning samples
  nt<-n+nburn
  
  # For AR(p) model, need to generate p points first
  Ep<-matrix(NA,p,2)
  colnames(Ep)<-c("E1","E2")
  # initial error
  s0<-svd(Delta)
  u<-rnorm(2)
  Ep[1,]<-s0$u%*%diag(sqrt(s0$d))%*%t(s0$v)%*%u
  s<-svd(Sigma)
  
  if(p>1)
  {
    for(j in 2:p)
    {
      u<-rnorm(2)
      e<-s$u%*%diag(sqrt(s$d))%*%t(s$v)%*%u
      
      Ep[j,]<-t(W[1:(2*(j-1)),])%*%c(t(Ep[(j-1):1,]))+e
    } 
  }
  
  E<-matrix(0,nt+p,2)
  E[1:p,]<-Ep
  for(j in (p+1):nrow(E))
  {
    u<-rnorm(2)
    e<-s$u%*%diag(sqrt(s$d))%*%t(s$v)%*%u
    
    E[j,]<-t(W)%*%c(t(E[(j-1):(j-p),]))+e
  }
  
  E1<-E[(nburn+p+1):(nrow(E)),1]
  E2<-E[(nburn+p+1):(nrow(E)),2]
  
  M<-Z*A+E1
  R<-Z*C+M*B+E2
  
  re1<-data.frame(Z=Z,M=M,R=R)
  re2<-data.frame(E1=E1,E2=E2)
  re<-list(data=re1,error=re2)
  return(re)
}

### ---- from gma/gma/R/sim.data.ts.two.R ----
sim.data.ts.two <-
function(Z.list,N,theta,Sigma,W,Delta=NULL,p=NULL,Lambda=diag(rep(1,3)),nburn=100)
{
  n<-rep(NA,N)
  for(i in 1:N)
  {
    n[i]<-length(Z.list[[i]])
  }
  
  s.Lambda<-svd(Lambda)
  Lambda.root<-s.Lambda$u%*%diag(sqrt(s.Lambda$d))%*%t(s.Lambda$v)
  
  eta<-matrix(rnorm(3*N),nrow=N)%*%Lambda.root
  A<-theta[1]+eta[,1]
  B<-theta[2]+eta[,2]
  C<-theta[3]+eta[,3]
  
  dat<-list()
  error<-list()
  for(i in 1:N)
  {
    re.tmp<-sim.data.ts.single(n[i],Z.list[[i]],A[i],B[i],C[i],Sigma,W,Delta,p,nburn)
    dat[[i]]<-re.tmp$data
    error[[i]]<-re.tmp$error
  }
  
  re<-list(data=dat,error=error,A=A,B=B,C=C,type="twolevel")
  return(re)
}

