##############################################################################
# MedMethods method module: spcma
# Sparse principal component based mediation analysis
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from spcma/spcma-master/R/BC.CI.R ----
BC.CI <-
function(theta,sims,conf.level=0.95) 
{
  low <- (1 - conf.level)/2
  high <- 1 - low
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- (1/6) * (sum(U^2))^{3/2}
  a <- top/under
  lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower2 <- lower <- quantile(theta, lower.inv)
  upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper2 <- upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
}

### ---- from spcma/spcma-master/R/deCor.R ----
deCor <-
function(X)
{
  n<-nrow(X)
  q<-ncol(X)
  
  if(q>1)
  {
    Xnew<-matrix(NA,n,q)
    Xnew[,1]<-X[,1]
    for(j in 2:q)
    {
      xtmp<-X[,1:(j-1)]
      fit<-lm(X[,j]~xtmp)
      Xnew[,j]<-fit$residuals+fit$coefficients[1]
    }
  }else
  {
    Xnew<-X
  }
  return(Xnew)
}

### ---- from spcma/spcma-master/R/deCorM.X.R ----
deCorM.X <-
function(M.tilde,X)
{
  n<-nrow(M.tilde)
  q<-ncol(M.tilde)
  
  if(q>1)
  {
    Mnew<-matrix(NA,n,q)
    Mnew[,1]<-M.tilde[,1]
    for(j in 2:q)
    {
      xtmp<-M.tilde[,1:(j-1)]
      fit<-lm(M.tilde[,j]~X+xtmp)
      Mnew[,j]<-fit$residuals+cbind(rep(1,n),X)%*%fit$coefficients[c(1,2)]
    }
  }else
  {
    Mnew<-M.tilde
  }
  return(Mnew)
}

### ---- from spcma/spcma-master/R/mcma_BK.R ----
mcma_BK <-
function(X,M,Y,sims=1000,boot=TRUE,boot.ci.type=c("bca","perc"),conf.level=0.95,p.adj.method=c("BH","bonferroni","BY"))
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  alpha<-matrix(NA,p,5)
  colnames(alpha)<-c("Estimate","pvalue","LB","UB","adjpv")
  rownames(alpha)<-colnames(M)
  beta=gamma=IE<-alpha
  
  # boot==TRUE
  IE.sims<-matrix(NA,sims,p)
  TE.sims<-matrix(NA,sims,p)
  DE.sims<-rep(NA,sims)
  # boot==FALSE
  IE.se<-rep(NA,p)
  for(j in 1:p)
  {
    dat.tmp<-data.frame(X=X,M=M[,j],Y=Y)
    
    fit.m<-lm(M~X,data=dat.tmp)
    fit.y<-lm(Y~X+M,data=dat.tmp)
    
    alpha[j,1:4]<-c(coef(fit.m)[2],summary(fit.m)$coefficients[2,4],confint(fit.m,level=conf.level)[2,])
    gamma[j,1:4]<-c(coef(fit.y)[2],summary(fit.y)$coefficients[2,4],confint(fit.y,level=conf.level)[2,])
    beta[j,1:4]<-c(coef(fit.y)[3],summary(fit.y)$coefficients[3,4],confint(fit.y,level=conf.level)[3,])
    
    if(boot)
    {
      re.med<-mediate(fit.m,fit.y,treat="X",mediator="M",sims=sims,boot=boot,boot.ci.type=boot.ci.type[1],conf.level=conf.level)
      
      IE[j,1:4]<-c(re.med$d1,re.med$d1.p,re.med$d1.ci)
      gamma[j,1:4]<-c(re.med$z1,re.med$z1.p,re.med$z1.ci) 
      
      IE.sims[,j]<-re.med$d1.sims
      TE.sims[,j]<-re.med$d1.sims+re.med$z1.sims
      IE.se[j]<-sd(re.med$d1.sims)
    }else
    {
      IE[j,1]<-alpha[j,1]*beta[j,1]
      IE.se[j]<-sqrt((alpha[j,1]*summary(fit.y)$coefficients[3,2])^2+(summary(fit.m)$coefficients[2,2]*beta[j,1])^2)
      IE[j,2]<-2*pnorm(abs(IE[j,1]/IE.se[j]),lower.tail=FALSE)
      IE[j,c(3,4)]<-c(IE[j,1]-IE.se[j]*qnorm((1-conf.level)/2,lower.tail=FALSE),IE[j,1]+IE.se[j]*qnorm((1-conf.level)/2,lower.tail=FALSE))
    }
  }
  alpha[,5]<-p.adjust(alpha[,2],method=p.adj.method[1])
  beta[,5]<-p.adjust(beta[,2],method=p.adj.method[1])
  gamma[,5]<-p.adjust(gamma[,2],method=p.adj.method[1])
  IE[,5]<-p.adjust(IE[,2],method=p.adj.method[1])
  
  IE.total<-matrix(NA,1,4)
  colnames(IE.total)<-c("Estiamte","pvalue","LB","UB")
  rownames(IE.total)<-"Total"
  DE<-matrix(NA,1,4)
  colnames(DE)<-c("Estimate","pvalue","LB","UB")
  rownames(DE)<-"DE"
  IE.total[1,1]<-sum(IE[,1])
  IE.total.se<-sqrt(sum(IE.se^2))
  IE.total[1,2]<-2*pnorm(abs(IE.total[1,1]/IE.total.se),lower.tail=FALSE)
  if(boot)
  {
    DE.sims<-apply(TE.sims,1,mean,na.rm=TRUE)-apply(IE.sims,1,sum,na.rm=TRUE)
    DE[1,1]<-mean(DE.sims,na.rm=TRUE)
    DE[1,2]<-2*pnorm(abs(DE[1,1]/sd(DE.sims,na.rm=TRUE)),lower.tail=FALSE)
    
    if(boot.ci.type[1]=="bca")
    {
      IE.total[1,c(3,4)]<-BC.CI(apply(IE.sims,1,sum),sims=sims,conf.level=conf.level)
      DE[1,c(3,4)]<-BC.CI(DE.sims,sims=sims,conf.level=conf.level)
    }else
    {
      IE.total[1,c(3,4)]<-quantile(apply(IE.sims,1,sum),probs=c((1-conf.level)/2,1-(1-conf.level)/2))
      DE[1,c(3,4)]<-quantile(DE.sims,probs=c((1-conf.level)/2,1-(1-conf.level)/2))
    }
  }else
  {
    IE.total[1,c(3,4)]<-c(IE.total[1,1]-IE.total.se*qnorm((1-conf.level)/2,lower.tail=FALSE),IE.total[1,1]+IE.total.se*qnorm((1-conf.level)/2,lower.tail=FALSE))
    
    # total effect model
    fit<-lm(Y~X)
    TE.se<-sqrt(summary(fit)$cov.unscaled[2,2]*(summary(fit)$sigma)^2)
    # direct effect
    DE[1,1]<-fit$coefficients[2]-IE.total[1,1]
    DE.se<-sqrt(IE.total.se^2+TE.se^2)
    DE[1,2]<-2*pnorm(abs(DE[1,1]/DE.se),lower.tail=FALSE)
    DE[1,c(3,4)]<-c(DE[1,1]-DE.se*qnorm((1-conf.level)/2,lower.tail=FALSE),DE[1,1]+DE.se*qnorm((1-conf.level)/2,lower.tail=FALSE))
  }
  
  re<-list(IE=IE,DE=DE,alpha=alpha,beta=beta,gamma=gamma,IE.total=IE.total)
  return(re)
}

### ---- from spcma/spcma-master/R/mcma_PCA.R ----
mcma_PCA <-
function(X,M,Y,adaptive=FALSE,var.per=0.8,n.pc=NULL,boot=TRUE,sims=1000,boot.ci.type=c("bca","perc"),conf.level=0.95,p.adj.method=c("BH","bonferroni","BY"))
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  # PCA on M~X residuals
  fit.m<-lm(M~X)
  Sigma.m<-cov(fit.m$residuals)
  svd.m<-svd(Sigma.m)
  if(adaptive)
  {
    n.pc<-sum(cumsum(svd.m$d)/sum(svd.m$d)<var.per)+1
  }else
  {
    if(is.null(n.pc))
    {
      n.pc<-p
    }
  }
  U<-svd.m$u[,1:n.pc]
  colnames(U)<-paste0("PC",1:n.pc)
  rownames(U)<-colnames(M)
  M.pc<-M%*%U
  colnames(M.pc)<-paste0("PC",1:n.pc)
  
  # run marginal mediation on PCs
  re.pc<-mcma_BK(X,M.pc,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
  re.pc$U<-U
  re.pc$var.per<-cumsum(svd.m$d[1:n.pc])/sum(svd.m$d)
  
  return(re.pc)
}

### ---- from spcma/spcma-master/R/plot_spcma.R ----
plot_spcma <-
function(object,plot.coef=c("alpha","beta","IE"),cex.lab=1,cex.axis=1,pt.cex=1,...)
{
  plot.idx<-which(names(object)==plot.coef)
  
  out<-as.matrix(object[[plot.idx]][,c(1,3,4)])
  colnames(out)<-c("Estimate","LB","UB")
  
  K<-nrow(out)
  
  sig<-as.numeric(out[,2]*out[,3]>0)
  neg<-as.numeric(out[,1]<0)*sig*3
  pos<-as.numeric(out[,1]>0)*sig*2
  pt.type<-(1-sig)+pos+neg
  col.tmp<-c(1,2,4)[pt.type]
  pt.tmp<-c(19,17,15)[pt.type]
  
  par(mar=c(5,5,3,3))
  plot(range(1-0.5/K,K+0.5/K),range(out[,c(2,3)],na.rm=TRUE),type="n",xaxt="n",xlab="",ylab=plot.coef,cex.lab=cex.lab,cex.axis=cex.axis)
  axis(side=1,at=1:K,labels=rownames(out),cex.axis=cex.axis)
  abline(h=0,lty=2,lwd=2,col=8)
  points(1:K,out[,1],pch=pt.tmp,col=col.tmp,cex=pt.cex)
  for(j in 1:K)
  {
    lines(rep(j,2),out[j,c(2,3)],lty=2,lwd=1,col=col.tmp[j])
    lines(c(j-0.3/K,j+0.3/K),rep(out[j,2],2),lty=1,lwd=1,col=col.tmp[j])
    lines(c(j-0.3/K,j+0.3/K),rep(out[j,3],2),lty=1,lwd=1,col=col.tmp[j])
  }
}

### ---- from spcma/spcma-master/R/R2.flasso.R ----
R2.flasso <-
function(E,U,D=NULL,gamma=0,eps=1e-4,maxsteps=2000,per.jump=0.7)
{
  p<-ncol(E)
  n.pc<-ncol(U)
  
  # PC
  E.pc<-E%*%U
  
  # % of variance
  var.total<-sum(diag(cov(E)))
  var.per=var.per.ind<-rep(NA,n.pc)
  for(j in 1:n.pc)
  {
    var.per.ind[j]<-var(E.pc[,j])/var.total
    var.per[j]<-sum(diag(cov(matrix(E.pc[,1:j],ncol=j))))/var.total
  }
  
  lambda.est<-rep(NA,n.pc)
  V<-matrix(NA,p,n.pc)
  var.per.new<-rep(NA,n.pc)
  # first PC lambda choice
  if(is.null(D))
  {
    out.tmp<-fusedlasso1d(y=E.pc[,1],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)  
  }else
  {
    out.tmp<-fusedlasso(y=E.pc[,1],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
  }
  var.per.tmp<-rep(NA,length(out.tmp$lambda))
  for(k in 1:length(out.tmp$lambda))
  {
    var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
  }
  var.per.diff.tmp<-var.per.tmp[2:length(var.per.tmp)]-var.per.tmp[1:(length(var.per.tmp)-1)]
  lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
  # lambda.idx.tmp<-min(which(abs(var.per.tmp-var.per.ind[1])<var.diff))
  lambda.est[1]<-out.tmp$lambda[lambda.idx.tmp]
  V[,1]<-out.tmp$beta[,lambda.idx.tmp]
  var.per.new[1]<-var.per.tmp[lambda.idx.tmp]
  
  for(j in 2:n.pc)
  {
    Etmp<-deCor(cbind(E%*%V[,1:(j-1)],E.pc[,j]))
    
    if(is.null(D))
    {
      out.tmp<-fusedlasso1d(y=Etmp[,j],X=E,gamma=gamma,eps=eps,maxsteps=maxsteps)  
    }else
    {
      out.tmp<-fusedlasso(y=Etmp[,j],X=E,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps)
    }
    
    var.per.tmp=var.per.tol.tmp<-rep(NA,length(out.tmp$lambda))
    for(k in 1:length(out.tmp$lambda))
    {
      var.per.tmp[k]<-var(out.tmp$fit[,k])/var.total
      dtmp<-deCor(cbind(E%*%V[,1:(j-1)],out.tmp$fit[,k]))
      var.per.tol.tmp[k]<-sum(diag(cov(dtmp)))/var.total
    }
    var.per.diff.tmp<-var.per.tol.tmp[2:length(var.per.tol.tmp)]-var.per.tol.tmp[1:(length(var.per.tol.tmp)-1)]
    lambda.idx.tmp<-max(which(var.per.diff.tmp>quantile(var.per.diff.tmp,probs=per.jump)))+1
    lambda.est[j]<-out.tmp$lambda[lambda.idx.tmp]
    V[,j]<-out.tmp$beta[,lambda.idx.tmp]
    var.per.new[j]<-var.per.tol.tmp[lambda.idx.tmp]
  }
  
  re<-list(lambda=lambda.est,V=V,var.per=var.per.new)
  return(re)
}

### ---- from spcma/spcma-master/R/SPCA.R ----
SPCA <-
function(X,M,adaptive=FALSE,var.per=0.8,n.pc=NULL,D=NULL,gamma=0,eps=1e-4,trace=TRUE,maxsteps=2000,lambda.tune=c("R2"),per.jump=0.7)
{
  n<-nrow(M)
  p<-ncol(M)
  
  if(is.null(colnames(M)))
  {
    colnames(M)<-paste0("M",1:p)
  }
  
  #==================================================
  # PCA on M~X residuals
  fit.m<-lm(M~X)
  E.m<-fit.m$residuals
  Sigma.m<-cov(E.m)
  svd.m<-svd(Sigma.m)
  if(adaptive)
  {
    n.pc<-sum(cumsum(svd.m$d)/sum(svd.m$d)<var.per)+1
  }else
  {
    if(is.null(n.pc))
    {
      n.pc<-p
    }
  }
  U<-svd.m$u[,1:n.pc]
  colnames(U)<-paste0("PC",1:n.pc)
  rownames(U)<-colnames(M)
  M.pc<-M%*%U
  colnames(M.pc)<-paste0("PC",1:n.pc)
  
  var.pc<-cumsum(svd.m$d[1:n.pc])/sum(svd.m$d)
  #==================================================
  
  
  #==================================================
  # Sparse PCA
  if(lambda.tune[1]=="R2")
  {
    # lambda chosen by adjusted total variance
    re.SPCA<-R2.flasso(E.m,U,D=D,gamma=gamma,eps=eps,maxsteps=maxsteps,per.jump=per.jump)
    V<-re.SPCA$V
    SPCA.var.per.cum<-re.SPCA$var.per
  }
  W<-apply(V,2,function(x){return(x/sqrt(sum(x^2)))})
  colnames(V)=colnames(W)<-paste0("PC",1:n.pc)
  rownames(V)=rownames(W)<-colnames(M)
  #==================================================
  
  # U: original loading matrix
  # V: sparsify loading matrix
  # W: sparsify loading matrix with l2-norm 1
  re<-list(U=U,V=V,W=W,var.pc=var.pc,var.spc=SPCA.var.per.cum)
  return(re)
}

### ---- from spcma/spcma-master/R/spcma-internal.R ----
.Random.seed <-
c(403L, 248L, 538071540L, 1176134445L, 2082784759L, 323466785L, 
-1179084428L, 1830411057L, 1706295580L, 1110474386L, 1339855211L, 
1585077755L, -170447494L, -952512894L, -122957877L, 633849249L, 
432221362L, 662259443L, -1435115883L, 12894496L, -232709384L, 
-1297662106L, -1735291474L, -131816034L, 6581051L, 1732690847L, 
-657118474L, 765149223L, 724178103L, -525154710L, 1047529867L, 
676529505L, 91739612L, -1127129816L, -213760579L, -275065960L, 
1680776327L, 1808679707L, -2106833994L, -850141157L, -905081316L, 
1719380743L, -1038350206L, -1936189551L, 1893414233L, 588685935L, 
1919925819L, -325408983L, 1670172108L, 101105528L, 1278652119L, 
435076140L, 1314766660L, 924868423L, 928004938L, 158823681L, 
206669412L, 1598619609L, 1259335545L, 1774957547L, 1483567182L, 
1353494775L, -1753689613L, -1938080736L, -1066435566L, 2087153273L, 
53157131L, 503627342L, -1330551436L, -1688777243L, 205887093L, 
1527807940L, -1419975209L, 1807974067L, 1752200916L, -977631203L, 
-116797827L, -1389955931L, -1685307790L, 933201236L, 191840473L, 
1964864906L, -1454156422L, 893742670L, 1860645402L, 694828895L, 
1744168099L, 1514881170L, 38987062L, 787264814L, 1793520973L, 
-62606897L, -1998941542L, -1085889480L, 111191082L, 791456695L, 
1992435453L, 1646632836L, 1251027280L, -1457845358L, -2071062218L, 
1322489803L, -226793130L, 1112346633L, 1951535391L, -1537489985L, 
-2052785841L, -548573942L, -1906552757L, -2129464602L, 1587784768L, 
552730002L, -1246934706L, 449292110L, -1824516560L, -1747617400L, 
-155269215L, -1433546060L, 1693582942L, 1036541281L, 1477912545L, 
-1765582174L, -688690940L, 2005618672L, 656538446L, -181440832L, 
-1484168974L, -757835148L, -461561136L, -476724688L, -1584597305L, 
-346338775L, -1866596378L, 1933346041L, 544962344L, 1865633851L, 
-1567727937L, 2106359601L, 1356699719L, -1891953678L, 1275785504L, 
1140802495L, -1263474636L, -1319066498L, -617395534L, 1745431654L, 
1146186614L, 26496608L, -308004774L, 1075866070L, -529085555L, 
-1590161765L, 1870678487L, 1502443032L, 1150643871L, -742902336L, 
-1971736457L, -1326404003L, -1507719490L, 1531139969L, -851984896L, 
2079350301L, -1565681410L, -155058829L, 1950279729L, 706924804L, 
121167440L, 2095414058L, 549844794L, -1343793727L, -1181703509L, 
-1143906319L, 95235037L, -334359551L, 32869335L, -144576273L, 
-1252341916L, 1854121816L, 1755886549L, 822356589L, 195731962L, 
-1533132911L, -189270605L, -1425288456L, 1300412129L, 1956227146L, 
1479593054L, -1663350663L, 1230548574L, -1209257150L, 2111769312L, 
-172717470L, 322231116L, -1294952002L, 1098692334L, -1387443004L, 
1731640158L, 614886006L, 1073409917L, 1477459505L, 1228470789L, 
1115593006L, 176058398L, 1546939324L, -860428786L, 1485384224L, 
-496838515L, -1420322348L, -990566425L, 171062082L, -455610800L, 
692922553L, 1198079513L, 1420127451L, 1692852003L, 781673440L, 
-315063617L, -1268653032L, -1604125628L, 1692597249L, 1960699071L, 
-2123619649L, 50669491L, 242619931L, -1240447704L, 1873827496L, 
-1162702982L, -539772494L, 1931268741L, -1306380423L, 1462100075L, 
1342036389L, -666812672L, -1892649518L, 1453163853L, -2028431682L, 
-537152332L, 2137182957L, 1842731883L, 441199136L, -475785272L, 
1091411157L, 2084592457L, -847646447L, -530480191L, 588262337L, 
1212066286L, -211378602L, -1013130941L, -791843768L, -241598403L, 
-8609623L, -1372983082L, 1278046117L, -672689420L, -755359248L, 
1242252224L, -49299708L, -2051699592L, 184668555L, -336662249L, 
1080732679L, 480617321L, -1557075590L, -11178943L, 835644761L, 
1359039182L, -1285113545L, 1043312607L, -85220025L, 548006636L, 
-1889086200L, -1077683985L, 431867113L, 644290476L, 1637022453L, 
-970784877L, 1386923606L, 476644299L, -1046266568L, 932505612L, 
1007592304L, -855467754L, 654815830L, -1427470603L, -2133679564L, 
575645745L, -712927236L, -1774049783L, -1203340042L, 965275096L, 
-686441010L, 1346242251L, -982680186L, -528026452L, 339753197L, 
-428001657L, -381794375L, -169061728L, 1551762310L, 707467000L, 
-1688197804L, -989918302L, 241311696L, -1237138615L, 170374059L, 
-1541721020L, -869816914L, -1654979358L, 777797423L, 54536234L, 
48712509L, -613331588L, -1801075260L, -1926958535L, 906390603L, 
-1132236896L, -103796852L, -1043254452L, 1403593900L, -775642436L, 
375518001L, -380051053L, 591069931L, 83770570L, -676980075L, 
689818398L, -1609368059L, 954794652L, -1969394180L, 1254077510L, 
840785907L, 82387434L, 855592615L, -1967989769L, 1075391698L, 
1168771949L, -1502989393L, 583366897L, -1968974206L, 1524507151L, 
17438409L, -646537606L, -287763427L, -378312654L, 818765225L, 
-978828002L, 612119292L, -2085313403L, -961468452L, -1051813186L, 
-237590278L, -124259256L, 1256640773L, 624226070L, -775425900L, 
-1819387058L, 15360305L, -471503820L, -2033401701L, 837525466L, 
574370849L, -1803874502L, 1651433510L, 1272817487L, -412682688L, 
-1491404709L, -1248665297L, -1012437485L, 784364637L, -240183790L, 
-1080079117L, 1922198096L, 1681508726L, -401733320L, 325068966L, 
-1306288225L, 1382069441L, -651434992L, -409534646L, -1636482750L, 
-1066055518L, -131384063L, 276909441L, 280997209L, 1185633815L, 
-1082864496L, 1869538773L, 931097448L, 872492226L, -829875528L, 
-1662150129L, -390951345L, 1862920803L, 1136089931L, 1927585456L, 
1729561581L, 844226584L, -939958104L, -795014212L, 540602054L, 
-808560971L, 138876158L, -530575617L, 1019046677L, 913705255L, 
-538514551L, 1178002945L, 29221803L, -1609447908L, 2099837081L, 
482402413L, 702405990L, -1973997468L, -1630773450L, -158163215L, 
1655561075L, -1944762559L, -1621740089L, 1671393844L, 775808302L, 
2014214586L, 3998534L, -1521966409L, -563401013L, 997541035L, 
513422494L, 1771352371L, 1283310762L, 1906925214L, 151475696L, 
736719149L, -665220628L, -1651483141L, 1473845133L, -1722442114L, 
-1197500083L, 1366547510L, -1457627424L, 60363316L, 43200650L, 
1861272758L, -1706639334L, 1016753833L, 1923784892L, -557600598L, 
1521686733L, 671020988L, 1744195615L, 71956347L, 796198548L, 
-302783775L, 235280073L, -1433323968L, 1191525728L, -717570372L, 
-1176125630L, -1203431831L, -1543622758L, 1872299350L, 1184896106L, 
1286839954L, 173196302L, -413588450L, 201086111L, 697257052L, 
1484493637L, -414316515L, -1775394688L, 954677334L, 594320049L, 
-955525291L, 1598672457L, -1135163280L, -1799742970L, -43150359L, 
237914351L, -72558049L, 694336676L, 670526018L, 831522868L, -1905100177L, 
-528236426L, -1406652385L, 1867219958L, -803859258L, 1768855496L, 
1672871835L, -1242983817L, 1323412685L, 1329629805L, -387813855L, 
-1307044951L, -1368660480L, 340166664L, -55505174L, 1140413508L, 
27201207L, 453677945L, 521886433L, -665810098L, 4016197L, 342195301L, 
-1732817932L, -1378503565L, 1490831349L, 541820992L, -2126884625L, 
317934374L, -672590686L, -1267642253L, -55500243L, 597941928L, 
2062513336L, -490354280L, -904771404L, 336930117L, 966664642L, 
1599417221L, -472229919L, -1696184014L, -191232520L, -878881733L, 
1313233761L, -1035780184L, 2118318956L, 682531065L, 77210165L, 
-613297490L, 321030558L, -1772874959L, 1872184247L, -311588594L, 
1050962588L, -425660236L, 511810016L, 1090127110L, 1400446584L, 
-2028388481L, 956040449L, -2109919429L, -1706166872L, -2071969963L, 
78487528L, -47887899L, -428950572L, -1614535033L, 1676221079L, 
-2129134970L, -528217916L, 1031115392L, -343195780L, 596109906L, 
1802124788L, -301338401L, -352813148L, -500754138L, -1277810585L, 
-651729821L, 1341558406L, 1866971493L, -1739178316L, -1682253330L, 
-1700658498L, -1646568746L, -1363656034L, 1407158120L, 671236624L, 
1813923547L, -959356990L, 1909508525L, -1247467881L, 1141950823L, 
-679719823L, 1763159131L, 1386474448L, -2003915603L, 1865031196L, 
-2055647960L, -391084412L, 546288158L, 896368834L, -691657771L, 
1434236115L, -391079951L, 398088221L, -267506693L, -138495343L, 
-1345049265L, -247771344L, -872146733L, 1789655399L, 855332282L, 
-1064822893L, 1229620174L, 1137390922L, -914128058L, 778092622L, 
1831197316L, -217215778L, -1901299705L, -1185478006L, 402442674L, 
1022787585L, 720327019L, 318499644L, -1039228870L, -261101621L, 
-1744511304L, 555198114L, -1795170958L, 1945987476L, -810289240L, 
-1324627958L, 1803094076L, -595763103L, 2120256470L, 751927146L, 
1056607119L, 1930820760L, -1543258749L, -1805147226L, -1883024042L, 
1991622797L, 747424685L, -1745860526L)

### ---- from spcma/spcma-master/R/spcma.R ----
spcma <-
function(X,M,Y,adaptive=FALSE,var.per=0.8,n.pc=NULL,D=NULL,gamma=0,eps=1e-4,maxsteps=2000,per.jump=0.7,
                     boot=TRUE,sims=1000,boot.ci.type=c("bca","perc"),conf.level=0.95,
                     p.adj.method=c("BH","bonferroni","BY"),PC.run=TRUE)
{
  n<-nrow(M)
  p<-ncol(M)
  
  re.SPCA<-SPCA(X,M,adaptive=adaptive,var.per=var.per,n.pc=n.pc,D=D,gamma=gamma,eps=eps,trace=FALSE,maxsteps=maxsteps,lambda.tune="R2",per.jump=per.jump)
  n.pc<-ncol(re.SPCA$U)
  
  #==================================================
  # run marginal mediaiton on sparse PCs
  M.spc<-M%*%re.SPCA$W
  M.spc.dCor<-deCorM.X(M.spc,X)
  re.spc<-mcma_BK(X,M.spc.dCor,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
  re.spc$W<-re.SPCA$W
  re.spc$var.per<-re.SPCA$var.spc
  #==================================================
  
  if(PC.run)
  {
    #==================================================
    # run marginal mediation on PCs
    M.pc<-M%*%re.SPCA$U
    re.pc<-mcma_BK(X,M.pc,Y,sims=sims,boot=boot,boot.ci.type=boot.ci.type,conf.level=conf.level,p.adj.method=p.adj.method)
    re.pc$U<-re.SPCA$U
    re.pc$var.per<-re.SPCA$var.pc
    #==================================================
    
    re<-list(PCA=re.pc,SPCA=re.spc)
  }else
  {
    re<-list(SPCA=re.spc)
  }
  
  return(re)
}

