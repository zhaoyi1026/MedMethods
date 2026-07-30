##############################################################################
# MedMethods method module: cfma
# Causal functional mediation analysis
#
# Assembled by tools/build_medpkg.R from the original method sources.
# Sourced into a private environment at .onLoad (see R/zzz.R), so internal
# helper names may safely collide with those of other method modules.
# Do not edit by hand -- edit the source files and re-run the build script.
##############################################################################

### ---- from cfma/cfma-master/R/BC.CI.R ----
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

### ---- from cfma/cfma-master/R/cfma-internal.R ----
.Random.seed <-
c(403L, 360L, 1132518362L, 861990609L, -788278369L, -1163499403L, 
207948770L, 1510631621L, -1777181360L, 1314363300L, -694110329L, 
551023851L, 1385710750L, -1906451013L, -26091222L, 1447842144L, 
-1012060895L, 907695500L, 1400340326L, -253539393L, 1411737125L, 
1559534860L, 1764170235L, -1429599369L, 1939961739L, -1110354225L, 
-1415906346L, 673948497L, 1078888975L, 228414461L, 1923357435L, 
-553313098L, -2140495858L, 2004480967L, 1829422236L, -1365053415L, 
1470506860L, -18744450L, -1619672052L, 1379098890L, -1253906429L, 
1298489912L, -304356180L, -289146654L, -1195148306L, -949495766L, 
1333680023L, 572465889L, 226456753L, 1182495779L, -1376014771L, 
-1921465601L, 1152172087L, 1259582935L, -1886539752L, 198209311L, 
1072776922L, -1662741975L, 419400346L, -553078174L, 1686708454L, 
303688778L, -1358505629L, 917134745L, 1976243960L, 1664761244L, 
1245153047L, 563531475L, 891853036L, 589566061L, 1975064195L, 
-1227668591L, -936025415L, 630723535L, -55197193L, 1956310701L, 
-1684398075L, 1957684047L, -891358227L, 1349146354L, 1010320786L, 
-1921968124L, 2125031051L, 1037567433L, -132147654L, 638801870L, 
-1340663680L, -2007041461L, -1070782590L, -456115621L, -796668427L, 
888318457L, -632746157L, 1746623991L, 942053078L, 594426102L, 
1474113269L, -2031217363L, -1083050216L, 1692300315L, 1430042455L, 
-696423042L, 1097904277L, 898783504L, -395857674L, 612077208L, 
383296871L, 395594459L, 208169413L, 1563884215L, -528112394L, 
1875924872L, -626163811L, 764143946L, 1584618874L, 1388731811L, 
-857847353L, -1700218489L, -866187259L, 113720102L, -1390980537L, 
-1698263802L, 877811673L, -1694725936L, 536531435L, 1927996262L, 
1632510030L, 358627649L, 1207034176L, -851124435L, -842121803L, 
-1762116878L, -2966813L, -812784288L, -272579776L, -1841929859L, 
-1099590619L, -1645785246L, -803017629L, -946290084L, 868844573L, 
-1213623000L, -393011841L, 1409553822L, -967374111L, 280013470L, 
-2031368195L, 790245349L, -1652381307L, 1788521360L, -1861339489L, 
991252145L, -74387923L, 545548400L, -1389053950L, 682391587L, 
407443766L, 1353757976L, -318394186L, 1475321112L, 1556880580L, 
-1814008520L, -1880136204L, -1679481494L, 1718298823L, -1031126656L, 
-2068098037L, -2066955009L, -750103390L, 853694464L, 1498938411L, 
1272380378L, -1995962432L, 280298975L, -32699401L, 907811281L, 
-605135704L, -1499840175L, -1211847079L, 171056456L, 2008598963L, 
-1409289645L, 523211392L, -984733713L, 616751751L, 867163436L, 
-1315088100L, 1235139753L, 188054416L, -297911498L, -1415866495L, 
-1377907199L, 733700952L, 466120511L, 1779540584L, -326839108L, 
-643884590L, -915947622L, -2035944524L, -1995665584L, -588787349L, 
-1441991467L, -1401808896L, -2029889463L, -572956037L, -969469510L, 
-798819466L, -750623650L, -1585707825L, -1481857050L, -420268518L, 
362579887L, 1278499928L, 32773050L, -974961059L, -511979482L, 
-517879215L, -1551350319L, 737981257L, 258896073L, 732699983L, 
448876024L, 2125452654L, -715558475L, -320821664L, -1576298910L, 
-763549250L, -823241412L, 597271142L, 1396046714L, -752499765L, 
1818865484L, -1292453143L, 899050296L, -908126363L, -1388021139L, 
1502620753L, -1266784264L, 1642218970L, -1638318815L, -1964399261L, 
-1011863122L, -675100004L, 1223308265L, 968448608L, 2075177603L, 
-279582648L, -657616298L, -396462198L, -307650963L, -1069352794L, 
-2111968462L, 435337292L, 518736050L, 237449109L, -1648201877L, 
1810731116L, 1284939378L, 227288196L, -1564102992L, 951451168L, 
1869528808L, -1571647021L, 1157527433L, -1950656126L, -811124345L, 
-700373320L, -2039685532L, -844175057L, -1917941968L, -1193181872L, 
684137833L, 303392514L, 497204783L, 1207813122L, 1227790956L, 
-1740930001L, -1089619523L, -1703366730L, 1126261398L, 298050432L, 
-486248773L, -1419358311L, 268634294L, -1051357027L, 535293173L, 
269740641L, -1842562765L, -2010431765L, 961815475L, -1954325755L, 
-796406658L, 1955141636L, -1855599441L, 1639596017L, 1605277764L, 
-1865549768L, -335905228L, -1971423507L, -542553009L, -509743255L, 
-950623693L, 313623999L, -880117750L, -1366330246L, 1288472201L, 
-1426632767L, 1358220041L, 1009614162L, -541054725L, 798008276L, 
304175801L, 1768901676L, 476886361L, 1791003659L, 1034414786L, 
185137916L, 1016320999L, 1322184158L, -1673041301L, 1906051081L, 
-117608698L, 1298782629L, -1771559834L, -1670774844L, 975282160L, 
448707331L, -814120071L, -380540311L, 865728128L, 396748078L, 
-725243586L, 1763228333L, 1796689575L, -1362715216L, -1115312439L, 
618961223L, 1446809436L, 1292247054L, -1405505848L, -1395183636L, 
711542749L, -1042908586L, 494991906L, -311900162L, 568996910L, 
152280061L, 1711221995L, 662440914L, 93690619L, -916581236L, 
1715592456L, -638269503L, -589335206L, 1368821098L, 1176567731L, 
-205035268L, 473900086L, 1393099385L, -1042219353L, -2116354209L, 
1476191241L, -1281651323L, -834648258L, 2113389014L, -220632741L, 
-420497935L, 2089785925L, -2127457202L, -1045737215L, -42068869L, 
43867670L, 2034775208L, 821943522L, 343372937L, 1685805696L, 
-644639271L, -2116737980L, -2118514108L, -992327614L, -798866562L, 
1749396414L, -1866780086L, 759158593L, 826865066L, -1809113784L, 
1206653359L, 1051063989L, 359645148L, 1527586183L, 1316884993L, 
167381547L, -1883903728L, 493393562L, -1513170093L, -670310309L, 
-2126448801L, -1017495496L, 673451714L, 2126438402L, 1009139574L, 
403398355L, 2132148348L, 980603139L, -964203027L, 1289659531L, 
735886098L, -792662145L, 1338045751L, -1037763639L, -1304779677L, 
1260581533L, -746585923L, 979746817L, 1543619051L, 158109563L, 
-961372664L, -1071604820L, -2126706726L, 309245633L, 415523552L, 
738477438L, 150038873L, -25428360L, -325609448L, 1597564517L, 
-2011313171L, 171072605L, 1465983038L, -1150218247L, 1815732214L, 
1681394540L, 1147231209L, 1135888984L, -1199207L, -1909072273L, 
615173550L, -1552726651L, 1999617650L, 2093816649L, 2147422288L, 
-4731156L, -325872120L, 105943285L, 653860852L, 206598164L, -570221212L, 
-858071643L, 193724328L, -1333940112L, -1854746869L, -1476461102L, 
72019748L, 1969257662L, -1308966848L, -1739837953L, -1271610067L, 
998700476L, -619692210L, 1483779728L, 1556699011L, -218275984L, 
-829581331L, 1573401378L, -1453112029L, 372173993L, 1604432876L, 
1880770939L, 1186502239L, -134704388L, -1151891004L, 606183928L, 
-100420817L, -928852197L, 1415907523L, -1446405908L, 460722364L, 
-2007389251L, 351303196L, -1991814034L, 1623374763L, -1859291937L, 
2029234940L, -1550541182L, -633893458L, -1222338801L, 993842601L, 
-171438217L, -1784152229L, -1500177960L, 364098585L, -1928456365L, 
793484269L, 1334568914L, 1636088279L, -1680891690L, 956804845L, 
1442459052L, -920085794L, 1240266550L, 553001206L, -394368566L, 
1972074217L, -1394721452L, 1821785436L, 1690548766L, -1861859538L, 
-641425401L, 1117956775L, 315982103L, 71585282L, -827734504L, 
69596246L, 1043099093L, 43902617L, 19494225L, -919137059L, -1939731214L, 
-53138900L, 2036316389L, -1225582805L, -1857174600L, -1140365447L, 
82808181L, -1560984960L, -604141348L, 1930755637L, -722880767L, 
1053486769L, 1971519368L, 1306990021L, -1536417809L, -620513043L, 
-1072425446L, 1820959885L, -578764225L, 733412217L, 458272590L, 
287688382L, -1756259848L, 1536074420L, -254771346L, 1757801998L, 
-26945343L, -1732122659L, -1695653538L, -501453234L, 20240172L, 
-534193275L, -820292917L, 60796139L, 959601573L, -191877067L, 
-1845940047L, -2118132780L, -1794643706L, -378927337L, 1275515944L, 
635214249L, -1727421914L, 100574632L, 1643882659L, 1655542273L, 
-2060087605L, 823894380L, -1265318962L, 1703997402L, -1207081333L, 
256673605L, -1603286746L, 2101589351L, 1939846641L, 643979664L, 
-534237209L, 311169101L, 217213903L, -677501453L, 364361832L, 
2075521663L, -1834102584L, 1443100896L, 1785310422L, 975381966L, 
1143338391L, -39441045L, 1076314107L, -1753473673L, 1553165157L, 
464255752L, -105241267L, 545022675L, -342370295L, 1079039099L, 
-1613409648L, 16552657L, 2010843977L, -564079709L, 1319839015L, 
163750384L, -598009134L, -1634502247L, -111172241L, -2007556833L, 
745289584L, -158942852L, -351485074L, -2086456841L, -1014326301L, 
379238002L, -932450797L, -1924313975L, 1290823697L, 1931295247L, 
-1940672213L, -305652307L, -1572117408L, 2082413139L, 971161356L, 
-955001853L, -1013522102L, -1838059624L, -67396908L, 1365335583L, 
539068482L, 1124868714L, 165962543L)

### ---- from cfma/cfma-master/R/design.mat.R ----
design.mat <-
function(X,intercept=TRUE)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  # design matrix
  q0<-dim(X)[3]
  if(is.na(q0))
  {
    if(intercept)
    {
      q<-2
      
      W<-array(NA,c(N,ntp,q))
      W[,,1]<-matrix(1,nrow=N,ncol=ntp)
      W[,,2]<-X
    }else
    {
      q<-1
      
      W<-array(NA,c(N,ntp,q))
      W[,,1]<-X
    }
  }else
    if(intercept)
    {
      q<-q0+1
      
      W<-array(NA,c(N,ntp,q))
      W[,,1]<-matrix(1,nrow=N,ncol=ntp)
      W[,,2:q]<-X
    }else
    {
      q<-q0
      W<-X
    }
  
  return(list(q=q,W=W))
}

### ---- from cfma/cfma-master/R/designD.mat.R ----
designD.mat <-
function(X,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,delta.grid=1)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  # design matrix
  X.design<-design.mat(X,intercept=intercept)
  q<-X.design$q
  W<-X.design$W
  
  # time grids
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # lambda
  if(length(lambda1)==1)
  {
    lambda1<-rep(lambda1,q)
  }else
    if(length(lambda1)>q)
    {
      lambda1<-lambda1[1:q]
    }else
    {
      lambda1<-rep(lambda1[1],q)
    }
  if(length(lambda2)==1)
  {
    lambda2<-rep(lambda2,q)
  }else
    if(length(lambda2)>q)
    {
      lambda2<-lambda2[1:q]
    }else
    {
      lambda2<-rep(lambda2[1],q)
    }
  
  
  K1<-nbasis1*q
  K2<-nbasis2*q
  K<-(nbasis1*nbasis2)*q
  ######################################################
  U<-matrix(0,K,K)
  V<-matrix(0,K,K)
  D<-array(0,c(N,K,ntp))
  for(j in 1:q)
  {
    U1tmp<-apply(array(apply(basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    U2tmp<-apply(array(apply(Ld2.basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    U[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda1[j]*kronecker(U1tmp,U2tmp)
    
    V1tmp<-apply(array(apply(Ld2.basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    V2tmp<-apply(array(apply(basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    V[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda2[j]*kronecker(V1tmp,V2tmp)
    
    
    Ws<-array(NA,c(N,nbasis1,ntp))
    Wstmp<-array(NA,c(N,nbasis1,ntp))
    for(i in 1:ntp)
    {
      Wstmp[,,i]<-matrix(W[,i,j],nrow=N)%*%matrix(basis1[i,],ncol=nbasis1)
      
      rtmp<-max(i-delta.grid,1)
      
      Ws[,,i]<-apply(array(Wstmp[,,rtmp:i],c(N,nbasis1,length(rtmp:i))),c(1,2),int.func,timeinv=c(timegrids[rtmp],timegrids[i]),timegrids=timegrids[rtmp:i])
      
      D[,((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),i]<-kronecker(t(basis2[i,]),Ws[,,i])
    }
  }
  
  return(list(D=D,U=U,V=V))
}

### ---- from cfma/cfma-master/R/FDA.concurrent.CV.R ----
FDA.concurrent.CV <-
function(X,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,
                            lambda=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),nfolds=5,verbose=TRUE)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  idx.tmp<-sample(1:N,N,replace=FALSE)
  cv.idx<-split(idx.tmp,sort(idx.tmp%%nfolds))
  
  mse.fit<-matrix(NA,N,length(lambda))
  
  for(kk in 1:length(cv.idx))
  {
    idx.ts<-cv.idx[[kk]]
    idx.tr<-unlist(cv.idx[-kk])
    
    Ytmp<-Y[idx.tr,]
    Yts<-matrix(Y[idx.ts,],nrow=length(idx.ts))
    if(length(dim(X))==2)
    {
      Xtmp<-X[idx.tr,]
      Xts<-matrix(X[idx.ts,],nrow=length(idx.ts))
    }
    if(length(dim(X))==3)
    {
      Xtmp<-X[idx.tr,,]
      Xts<-array(X[idx.ts,,],c(length(idx.ts),dim(X)[2],dim(X)[3]))
    }
    
    for(ii in 1:length(lambda))
    {
      re<-FDA.concurrent(Xtmp,Ytmp,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda[ii])
      
      Xts.design<-design.mat(Xts,intercept=intercept)
      yfit<-t(apply(Xts.design$W,1,function(x){return(apply(x*t(re$gamma.curve),1,sum))}))
      
      mse.fit[idx.ts,ii]<-apply((yfit-Yts)^2,1,int.func,timeinv=timeinv,timegrids=timegrids)
    }
    
    if(verbose)
    {
      print(paste0("Fold ",kk))
    }
  }
  
  lambda.idx<-which.min(apply(mse.fit,2,mean))[1]
  lambda.est<-lambda[lambda.idx]
  re.fit<-FDA.concurrent(X,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda.est)
  re.fit$lambda<-lambda.est
  re.fit$mse<-mse.fit
  
  return(re.fit)
}

### ---- from cfma/cfma-master/R/FDA.concurrent.R ----
FDA.concurrent <-
function(X,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,lambda=0.01)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  # design matrix
  X.design<-design.mat(X,intercept=intercept)
  q<-X.design$q
  W<-X.design$W
  
  # time grids
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # basis functions
  if(is.null(basis))
  {
    if(basis.type[1]=="fourier")
    {
      basis<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis)
      
      Ld2.basis<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis)
    }
  }else
  {
    nbasis<-ncol(basis)
  }
  
  # lambda
  if(length(lambda)==1)
  {
    lambda<-rep(lambda,q)
  }else
    if(length(lambda)>q)
    {
      lambda<-lambda[1:q]
    }else
    {
      lambda<-rep(lambda[1],q)
    }
  
  K<-nbasis*q
  
  U<-matrix(0,K,K)
  Theta<-array(0,c(q,K,ntp))
  for(j in 1:q)
  {
    Utmp<-lambda[j]*apply(array(apply(Ld2.basis,1,function(x){return(x%*%t(x))}),c(nbasis,nbasis,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    U[((j-1)*nbasis+1):(j*nbasis),((j-1)*nbasis+1):(j*nbasis)]<-Utmp
    
    Theta[j,((j-1)*nbasis+1):(j*nbasis),]<-t(basis)
  }
  
  # W %*% Theta matrix
  mat1<-array(NA,c(K,K,ntp))
  mat2<-array(NA,c(K,1,ntp))
  for(s in 1:ntp)
  {
    dtmp<-matrix(W[,s,],ncol=q)%*%matrix(Theta[,,s],nrow=q)
    mat1[,,s]<-t(dtmp)%*%dtmp
    mat2[,,s]<-t(matrix(Theta[,,s],nrow=q))%*%t(matrix(W[,s,],ncol=q))%*%Y[,s]
  }
  V1<-apply(mat1,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  V2<-apply(mat2,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  g<-solve(V1+U)%*%V2
  
  # gamma.est<-matrix(apply(Theta,3,function(x){return(x%*%g)}),ncol=q)
  gamma.est<-matrix(NA,ntp,q)
  for(j in 1:q)
  {
    gamma.est[,j]<-basis%*%g[((j-1)*nbasis+1):(j*nbasis)]
  }
  if(intercept)
  {
    colnames(gamma.est)<-c("Intercept",paste0("X",1:(q-1))) 
  }else
  {
    colnames(gamma.est)<-paste0("X",1:q)
  }
  yfit<-t(apply(W,1,function(x){return(apply(x*gamma.est,1,sum))}))
  
  re<-list(coefficients=c(g),basis=basis,gamma.curve=t(gamma.est),fitted=yfit)
  
  return(re)
}

### ---- from cfma/cfma-master/R/FDA.historical.CV.R ----
FDA.historical.CV <-
function(X,Y,delta.grid=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,
                            lambda1=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),lambda2=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),nfolds=5,verbose=TRUE)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  idx.tmp<-sample(1:N,N,replace=FALSE)
  cv.idx<-split(idx.tmp,sort(idx.tmp%%nfolds))
  
  mse.fit<-matrix(NA,N,min(length(lambda1),length(lambda2)))
  
  for(kk in 1:length(cv.idx))
  {
    idx.ts<-cv.idx[[kk]]
    idx.tr<-unlist(cv.idx[-kk])
    
    Ytmp<-Y[idx.tr,]
    Yts<-matrix(Y[idx.ts,],nrow=length(idx.ts))
    if(length(dim(X))==2)
    {
      Xtmp<-X[idx.tr,]
      Xts<-matrix(X[idx.ts,],nrow=length(idx.ts))
    }
    if(length(dim(X))==3)
    {
      Xtmp<-X[idx.tr,,]
      Xts<-array(X[idx.ts,,],c(length(idx.ts),dim(X)[2],dim(X)[3]))
    }
    
    for(ii in 1:min(length(lambda1),length(lambda2)))
    {
      re<-FDA.historical(Xtmp,Ytmp,delta.grid=delta.grid,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                         nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1[ii],lambda2=lambda2[ii])
      
      XtsD.mat<-designD.mat(Xts,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,nbasis1=nbasis1,nbasis2=nbasis2,
                            timeinv=timeinv,timegrids=timegrids,delta.grid=delta.grid)
      yfit<-apply(XtsD.mat$D,c(1,3),function(x){return(t(x)%*%c(re$coef.vec))})
      
      mse.fit[idx.ts,ii]<-apply((yfit-Yts)^2,1,int.func,timeinv=timeinv,timegrids=timegrids)
    }
    
    if(verbose)
    {
      print(paste0("Fold ",kk))
    }
  }
  
  lambda.idx<-which.min(apply(mse.fit,2,mean))[1]
  lambda.est<-lambda1[lambda.idx]
  re.fit<-FDA.historical(X,Y,delta.grid=delta.grid,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                         nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda.est,lambda2=lambda.est)
  re.fit$lambda1<-lambda.est
  re.fit$lambda2<-lambda.est
  re.fit$mse<-mse.fit
  
  return(re.fit)
}

### ---- from cfma/cfma-master/R/FDA.historical.R ----
FDA.historical <-
function(X,Y,delta.grid=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),
                         nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,lambda1=0.01,lambda2=0.01)
{
  N<-dim(X)[1]             # # of subject
  ntp<-dim(X)[2]           # # of time points
  
  # design matrix
  X.design<-design.mat(X,intercept=intercept)
  q<-X.design$q
  W<-X.design$W
  
  # time grids
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # lambda
  if(length(lambda1)==1)
  {
    lambda1<-rep(lambda1,q)
  }else
    if(length(lambda1)>q)
    {
      lambda1<-lambda1[1:q]
    }else
    {
      lambda1<-rep(lambda1[1],q)
    }
  if(length(lambda2)==1)
  {
    lambda2<-rep(lambda2,q)
  }else
    if(length(lambda2)>q)
    {
      lambda2<-lambda2[1:q]
    }else
    {
      lambda2<-rep(lambda2[1],q)
    }
  
  
  K1<-nbasis1*q
  K2<-nbasis2*q
  K<-(nbasis1*nbasis2)*q
  ######################################################
  U<-matrix(0,K,K)
  V<-matrix(0,K,K)
  D<-array(0,c(N,K,ntp))
  for(j in 1:q)
  {
    U1tmp<-apply(array(apply(basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    U2tmp<-apply(array(apply(Ld2.basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    U[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda1[j]*kronecker(U1tmp,U2tmp)
    
    V1tmp<-apply(array(apply(Ld2.basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    V2tmp<-apply(array(apply(basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
    V[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda2[j]*kronecker(V1tmp,V2tmp)
    
    
    Ws<-array(NA,c(N,nbasis1,ntp))
    Wstmp<-array(NA,c(N,nbasis1,ntp))
    for(i in 1:ntp)
    {
      Wstmp[,,i]<-matrix(W[,i,j],nrow=N)%*%matrix(basis1[i,],ncol=nbasis1)
      
      rtmp<-max(i-delta.grid,1)
      
      Ws[,,i]<-apply(array(Wstmp[,,rtmp:i],c(N,nbasis1,length(rtmp:i))),c(1,2),int.func,timeinv=c(timegrids[rtmp],timegrids[i]),timegrids=timegrids[rtmp:i])
      
      D[,((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),i]<-kronecker(t(basis2[i,]),Ws[,,i])
    }
  }
  mat1<-array(NA,c(K,K,ntp))
  mat2<-array(NA,c(K,1,ntp))
  for(i in 1:ntp)
  {
    mat1[,,i]<-t(D[,,i])%*%D[,,i]
    mat2[,,i]<-t(D[,,i])%*%Y[,i]
  }
  Q1<-apply(mat1,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  Q2<-apply(mat2,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  g<-solve(Q1+U+V)%*%Q2
  Gmat<-matrix(0,K1,K2)
  for(j in 1:q)
  {
    Gmat[((j-1)*nbasis1+1):(j*nbasis1),((j-1)*nbasis2+1):(j*nbasis2)]<-matrix(g[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))],nbasis1,nbasis2)
  }
  ######################################################
  
  gamma.est<-array(NA,c(ntp,ntp,q))
  for(j in 1:q)
  {
    for(u in 1:ntp)
    {
      for(s in 1:ntp)
      {
        gamma.est[u,s,j]<-t(basis1[u,])%*%Gmat[((j-1)*nbasis1+1):(j*nbasis1),((j-1)*nbasis2+1):(j*nbasis2)]%*%basis2[s,]
      }
    }
  }
  if(intercept)
  {
    dimnames(gamma.est)[[3]]<-c("Intercept",paste0("X",1:(q-1))) 
  }else
  {
    dimnames(gamma.est)<-list(NULL)
    dimnames(gamma.est)[[3]]<-paste0("X",1:q)
  }
  yfit<-apply(D,c(1,3),function(x){return(t(x)%*%g)})
  
  re<-list(coefficients=Gmat,coef.vec=g,basis1=basis1,basis2=basis2,gamma.curve=gamma.est,fitted=yfit)
  
  return(re)
}

### ---- from cfma/cfma-master/R/FDA.historical2.CV.R ----
FDA.historical2.CV <-
function(X1,X2,Y,delta.grid1=1,delta.grid2=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,
                             lambda1=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),lambda2=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),nfolds=5,verbose=TRUE)
{
  N<-dim(X1)[1]             # # of subject
  ntp<-dim(X1)[2]           # # of time points
  
  idx.tmp<-sample(1:N,N,replace=FALSE)
  cv.idx<-split(idx.tmp,sort(idx.tmp%%nfolds))
  
  mse.fit<-matrix(NA,N,min(length(lambda1),length(lambda2)))
  
  for(kk in 1:length(cv.idx))
  {
    idx.ts<-cv.idx[[kk]]
    idx.tr<-unlist(cv.idx[-kk])
    
    Ytmp<-Y[idx.tr,]
    Yts<-matrix(Y[idx.ts,],nrow=length(idx.ts))
    if(length(dim(X1))==2)
    {
      X1tmp<-X1[idx.tr,]
      X1ts<-matrix(X1[idx.ts,],nrow=length(idx.ts))
    }
    if(length(dim(X1))==3)
    {
      X1tmp<-X1[idx.tr,,]
      X1ts<-array(X1[idx.ts,,],c(length(idx.ts),dim(X1)[2],dim(X1)[3]))
    }
    if(length(dim(X2))==2)
    {
      X2tmp<-X2[idx.tr,]
      X2ts<-matrix(X2[idx.ts,],nrow=length(idx.ts))
    }
    if(length(dim(X2))==3)
    {
      X2tmp<-X2[idx.tr,,]
      X2ts<-array(X2[idx.ts,,],c(length(idx.ts),dim(X2)[2],dim(X2)[3]))
    }
    
    for(ii in 1:min(length(lambda1),length(lambda2)))
    {
      re<-FDA.historical2(X1tmp,X2tmp,Ytmp,delta.grid1=delta.grid1,delta.grid2=delta.grid2,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                          nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1[ii],lambda2=lambda2[ii])
      
      X1tsD.mat<-designD.mat(X1ts,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,nbasis1=nbasis1,nbasis2=nbasis2,
                             timeinv=timeinv,timegrids=timegrids,delta.grid=delta.grid1)
      X2tsD.mat<-designD.mat(X2ts,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,nbasis1=nbasis1,nbasis2=nbasis2,
                             timeinv=timeinv,timegrids=timegrids,delta.grid=delta.grid2)
      D1<-X1tsD.mat$D
      D2<-X2tsD.mat$D
      D<-array(NA,c(dim(D1)[1],dim(D1)[2]+dim(D2)[2],dim(D1)[3]))
      D[,1:(dim(D1)[2]),]<-D1
      D[,(dim(D1)[2]+1):(dim(D1)[2]+dim(D2)[2]),]<-D2
      
      yfit<-apply(D,c(1,3),function(x){return(t(x)%*%c(re$coef.vec))})
      
      mse.fit[idx.ts,ii]<-apply((yfit-Yts)^2,1,int.func,timeinv=timeinv,timegrids=timegrids)
    }
    
    if(verbose)
    {
      print(paste0("Fold ",kk))
    }
  }
  
  lambda.idx<-which.min(apply(mse.fit,2,mean))[1]
  lambda.est<-lambda1[lambda.idx]
  re.fit<-FDA.historical2(X1,X2,Y,delta.grid1=delta.grid1,delta.grid2=delta.grid2,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                          nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda.est,lambda2=lambda.est)
  re.fit$lambda1<-lambda.est
  re.fit$lambda2<-lambda.est
  re.fit$mse<-mse.fit
  
  return(re.fit)
}

### ---- from cfma/cfma-master/R/FDA.historical2.R ----
FDA.historical2 <-
function(X1,X2,Y,delta.grid1=1,delta.grid2=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),
                          nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,lambda1=0.01,lambda2=0.01)
{
  N<-dim(X1)[1]             # # of subject
  ntp<-dim(X1)[2]           # # of time points
  
  # design matrix
  X.dim<-function(X)
  {
    q0<-dim(X)[3]
    if(is.na(q0))
    {
      if(intercept)
      {
        q<-2
        
        W<-array(NA,c(N,ntp,q))
        W[,,1]<-matrix(1,nrow=N,ncol=ntp)
        W[,,2]<-X
      }else
      {
        q<-1
        
        W<-array(NA,c(N,ntp,q))
        W[,,1]<-X
      }
    }else
      if(intercept)
      {
        q<-q0+1
        
        W<-array(NA,c(N,ntp,q))
        W[,,1]<-matrix(1,nrow=N,ncol=ntp)
        W[,,2:q]<-X
      }else
      {
        q<-q0
        W<-X
      }
    return(list(q=q,W=W))
  }
  X.dim.re1<-X.dim(X1)
  X.dim.re2<-X.dim(X2)
  
  q1<-X.dim.re1$q
  W1<-X.dim.re1$W
  q2<-X.dim.re2$q
  W2<-X.dim.re2$W
  
  q<-q1+q2
  
  # time grids
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # lambda
  if(length(lambda1)==1)
  {
    lambda1<-rep(lambda1,q)
  }else
    if(length(lambda1)>q)
    {
      lambda1<-lambda1[1:q]
    }else
    {
      lambda1<-rep(lambda1[1],q)
    }
  if(length(lambda2)==1)
  {
    lambda2<-rep(lambda2,q)
  }else
    if(length(lambda2)>q)
    {
      lambda2<-lambda2[1:q]
    }else
    {
      lambda2<-rep(lambda2[1],q)
    }
  
  
  K1<-c(nbasis1*q1,nbasis1*q2)
  K2<-c(nbasis2*q1,nbasis2*q2)
  K<-c((nbasis1*nbasis2)*q1,(nbasis1*nbasis2)*q2)
  ######################################################
  D.func<-function(W,K1,K2,K,q,delta.grid)
  {
    U<-matrix(0,K,K)
    V<-matrix(0,K,K)
    D<-array(0,c(N,K,ntp))
    for(j in 1:q)
    {
      U1tmp<-apply(array(apply(basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
      U2tmp<-apply(array(apply(Ld2.basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
      U[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda1[j]*kronecker(U1tmp,U2tmp)
      
      V1tmp<-apply(array(apply(Ld2.basis2,1,function(x){return(x%*%t(x))}),c(nbasis2,nbasis2,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
      V2tmp<-apply(array(apply(basis1,1,function(x){return(x%*%t(x))}),c(nbasis1,nbasis1,ntp)),c(1,2),function(x){return(int.func(x,timeinv=timeinv,timegrids=timegrids))})
      V[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))]<-lambda2[j]*kronecker(V1tmp,V2tmp)
      
      
      Ws<-array(NA,c(N,nbasis1,ntp))
      Wstmp<-array(NA,c(N,nbasis1,ntp))
      for(i in 1:ntp)
      {
        Wstmp[,,i]<-matrix(W[,i,j],nrow=N)%*%matrix(basis1[i,],ncol=nbasis1)
        
        rtmp<-max(i-delta.grid,1)
        
        Ws[,,i]<-apply(array(Wstmp[,,rtmp:i],c(N,nbasis1,length(rtmp:i))),c(1,2),int.func,timeinv=c(timegrids[rtmp],timegrids[i]),timegrids=timegrids[rtmp:i])
        
        D[,((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2)),i]<-kronecker(t(basis2[i,]),Ws[,,i])
      }
    }
    
    return(list(D=D,U=U,V=V))
  }
  re.D1<-D.func(W1,K1[1],K2[1],K[1],q1,delta.grid1)
  re.D2<-D.func(W2,K1[2],K2[2],K[2],q2,delta.grid2)
  
  D1<-re.D1$D
  D2<-re.D2$D
  
  D<-array(NA,c(dim(D1)[1],dim(D1)[2]+dim(D2)[2],dim(D1)[3]))
  D[,1:(dim(D1)[2]),]<-D1
  D[,(dim(D1)[2]+1):(dim(D1)[2]+dim(D2)[2]),]<-D2
  
  mat1<-array(NA,c(sum(K),sum(K),ntp))
  mat2<-array(NA,c(sum(K),1,ntp))
  for(i in 1:ntp)
  {
    mat1[,,i]<-t(D[,,i])%*%D[,,i]
    mat2[,,i]<-t(D[,,i])%*%Y[,i]
  }
  Q1<-apply(mat1,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  Q2<-apply(mat2,c(1,2),int.func,timeinv=timeinv,timegrids=timegrids)
  
  U1<-re.D1$U
  V1<-re.D1$V
  U2<-re.D2$U
  V2<-re.D2$V
  
  U=V<-matrix(0,sum(K),sum(K))
  U[1:K[1],1:K[1]]<-U1
  U[(K[1]+1):sum(K),(K[1]+1):sum(K)]<-U2
  V[1:K[1],1:K[1]]<-V1
  V[(K[1]+1):sum(K),(K[1]+1):sum(K)]<-V2
  
  g<-solve(Q1+U+V)%*%Q2
  g1<-g[1:K[1]]
  g2<-g[(K[1]+1):sum(K)]
  
  Gmat.func<-function(K1,K2,q,g)
  {
    Gmat<-matrix(0,K1,K2)
    for(j in 1:q)
    {
      Gmat[((j-1)*nbasis1+1):(j*nbasis1),((j-1)*nbasis2+1):(j*nbasis2)]<-matrix(g[((j-1)*(nbasis1*nbasis2)+1):(j*(nbasis1*nbasis2))],nbasis1,nbasis2)
    }
    
    return(Gmat)
  }
  Gmat1<-Gmat.func(K1[1],K2[1],q1,g1)
  Gmat2<-Gmat.func(K1[2],K2[2],q2,g2)
  
  Gmat<-matrix(0,dim(Gmat1)[1]+dim(Gmat2)[1],dim(Gmat1)[2]+dim(Gmat2)[2])
  Gmat[1:(dim(Gmat1)[1]),1:(dim(Gmat1)[2])]<-Gmat1
  Gmat[(dim(Gmat1)[1]+1):(dim(Gmat1)[1]+dim(Gmat2)[1]),(dim(Gmat1)[2]+1):(dim(Gmat1)[2]+dim(Gmat2)[2])]<-Gmat2
  ######################################################
  gamma.func<-function(Gmat,q)
  {
    gamma.est<-array(NA,c(ntp,ntp,q))
    for(j in 1:q)
    {
      for(u in 1:ntp)
      {
        for(s in 1:ntp)
        {
          gamma.est[u,s,j]<-t(basis1[u,])%*%Gmat[((j-1)*nbasis1+1):(j*nbasis1),((j-1)*nbasis2+1):(j*nbasis2)]%*%basis2[s,]
        }
      }
    }
    
    return(gamma.est)
  }
  gamma.est1<-gamma.func(Gmat1,q1)
  gamma.est2<-gamma.func(Gmat2,q2)
  
  gamma.est<-array(NA,c(ntp,ntp,q))
  gamma.est[,,1:q1]<-gamma.est1
  gamma.est[,,(q1+1):q]<-gamma.est2
  
  if(intercept)
  {
    dimnames(gamma.est)[[3]]<-c("Intercept",paste0("X",1:(q-1))) 
  }else
  {
    dimnames(gamma.est)<-list(NULL)
    dimnames(gamma.est)[[3]]<-paste0("X",1:q)
  }
  yfit<-apply(D,c(1,3),function(x){return(t(x)%*%g)})
  
  re<-list(coefficients=Gmat,coef.vec=g,coef1=Gmat1,coef2=Gmat2,basis1=basis1,basis2=basis2,gamma.curve=gamma.est,fitted=yfit)
  
  return(re)
}

### ---- from cfma/cfma-master/R/FMA.concurrent.boot.R ----
FMA.concurrent.boot <-
function(Z,M,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,lambda.m=0.01,lambda.y=0.01,
                              sims=1000,boot=TRUE,boot.ci.type=c("bca","perc"),conf.level=0.95,verbose=TRUE)
{
  N<-nrow(Z)             # # of subject
  ntp<-ncol(Z)           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # M model
  fit.m<-FDA.concurrent(Z,M,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda.m)
  # Y model
  Xtmp<-array(NA,c(N,ntp,2))
  Xtmp[,,1]<-Z
  Xtmp[,,2]<-M
  fit.y<-FDA.concurrent(Xtmp,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda.y)
  
  if(boot)
  {
    coef.alpha=coef.beta=coef.gamma=coef.IE<-matrix(NA,sims,ncol(fit.m$basis))
    c.alpha=c.beta=c.gamma=c.IE<-matrix(NA,sims,ntp)
    for(b in 1:sims)
    {
      idx.tmp<-sample(1:N,N,replace=TRUE)
      
      Ztmp<-Z[idx.tmp,]
      Mtmp<-M[idx.tmp,]
      Ytmp<-Y[idx.tmp,]
      
      re.tmp<-FMA.concurrent(Ztmp,Mtmp,Ytmp,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,
                             lambda.m=lambda.m,lambda.y=lambda.y)
      
      coef.alpha[b,]<-re.tmp$M$coefficients["Z",]
      coef.gamma[b,]<-re.tmp$Y$coefficients["Z",]
      coef.beta[b,]<-re.tmp$Y$coefficients["M",]
      
      coef.IE[b,]<-re.tmp$IE$coefficients
      
      c.alpha[b,]<-re.tmp$M$curve["Z",]
      c.gamma[b,]<-re.tmp$Y$curve["Z",]
      c.beta[b,]<-re.tmp$Y$curve["M",]
      c.IE[b,]<-re.tmp$IE$curve
      
      if(verbose)
      {
        print(paste0("Bootstrap sample ",b))
      }
    }
    se.alpha<-apply(coef.alpha,2,sd,na.rm=TRUE)
    se.gamma<-apply(coef.gamma,2,sd,na.rm=TRUE)
    se.beta<-apply(coef.beta,2,sd,na.rm=TRUE)
    se.IE<-apply(coef.IE,2,sd,na.rm=TRUE)
    
    se.c.alpha<-apply(c.alpha,2,sd,na.rm=TRUE)
    se.c.gamma<-apply(c.gamma,2,sd,na.rm=TRUE)
    se.c.beta<-apply(c.beta,2,sd,na.rm=TRUE)
    se.c.IE<-apply(c.IE,2,sd,na.rm=TRUE)
    
    if(boot.ci.type[1]=="bca")
    {
      ci.alpha<-apply(coef.alpha,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.gamma<-apply(coef.gamma,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.beta<-apply(coef.beta,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.IE<-apply(coef.IE,2,BC.CI,sims=sims,conf.level=conf.level)
      
      ci.c.alpha<-apply(c.alpha,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.gamma<-apply(c.gamma,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.beta<-apply(c.beta,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.IE<-apply(c.IE,2,BC.CI,sims=sims,conf.level=conf.level)
    }
    if(boot.ci.type[1]=="perc")
    {
      ci.alpha<-apply(coef.alpha,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.gamma<-apply(coef.gamma,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.beta<-apply(coef.beta,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.IE<-apply(coef.IE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      
      ci.c.alpha<-apply(c.alpha,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.gamma<-apply(c.gamma,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.beta<-apply(c.beta,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.IE<-apply(c.IE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
    }
    
    re.alpha<-rbind(apply(coef.alpha,2,mean,na.rm=TRUE),se.alpha,ci.alpha)
    re.gamma<-rbind(apply(coef.gamma,2,mean,na.rm=TRUE),se.gamma,ci.gamma)
    re.beta<-rbind(apply(coef.beta,2,mean,na.rm=TRUE),se.beta,ci.beta)
    re.IE<-rbind(apply(coef.IE,2,mean,na.rm=TRUE),se.IE,ci.IE)
    rownames(re.alpha)=rownames(re.gamma)=rownames(re.beta)=rownames(re.IE)<-c("Estimate","SE","LB","UB")
    colnames(re.alpha)=colnames(re.gamma)=colnames(re.beta)=colnames(re.IE)<-paste0("basis",1:ncol(fit.m$basis))
    
    curve.alpha<-rbind(apply(c.alpha,2,mean,na.rm=TRUE),se.c.alpha,ci.c.alpha)
    curve.gamma<-rbind(apply(c.gamma,2,mean,na.rm=TRUE),se.c.gamma,ci.c.gamma)
    curve.beta<-rbind(apply(c.beta,2,mean,na.rm=TRUE),se.c.beta,ci.c.beta)
    curve.IE<-rbind(apply(c.IE,2,mean,na.rm=TRUE),se.c.IE,ci.c.IE)
    rownames(curve.alpha)=rownames(curve.gamma)=rownames(curve.beta)=rownames(curve.IE)<-c("Estimate","SE","LB","UB")
    
    re<-list(alpha=list(coefficients=re.alpha,curve=curve.alpha),gamma=list(coefficients=re.gamma,curve=curve.gamma),beta=list(coefficients=re.beta,curve=curve.beta),
             IE=list(coefficients=re.IE,curve=curve.IE),DE=list(coefficients=re.gamma,curve=curve.gamma))
    
    return(re)
  }else
  {
    return(FMA.concurrent(Z,M,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda.m=lambda.m,lambda.y=lambda.y))
  }
}

### ---- from cfma/cfma-master/R/FMA.concurrent.CV.boot.R ----
FMA.concurrent.CV.boot <-
function(Z,M,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,
                                 lambda=10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1]),nfolds=5,
                                 sims=1000,boot=TRUE,boot.ci.type=c("bca","perc"),conf.level=0.95,verbose=TRUE)
{
  N<-nrow(Z)             # # of subject
  ntp<-ncol(Z)           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # M model
  fit.m<-FDA.concurrent.CV(Z,M,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda,nfolds=nfolds,verbose=FALSE)
  lambda.m<-fit.m$lambda
  # Y model
  Xtmp<-array(NA,c(N,ntp,2))
  Xtmp[,,1]<-Z
  Xtmp[,,2]<-M
  fit.y<-FDA.concurrent.CV(Xtmp,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda,nfolds=nfolds,verbose=FALSE)
  lambda.y<-fit.y$lambda
  
  if(boot)
  {
    coef.alpha=coef.beta=coef.gamma=coef.IE<-matrix(NA,sims,ncol(fit.m$basis))
    c.alpha=c.beta=c.gamma=c.IE<-matrix(NA,sims,ntp)
    for(b in 1:sims)
    {
      idx.tmp<-sample(1:N,N,replace=TRUE)
      
      Ztmp<-Z[idx.tmp,]
      Mtmp<-M[idx.tmp,]
      Ytmp<-Y[idx.tmp,]
      
      re.tmp<-FMA.concurrent(Ztmp,Mtmp,Ytmp,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,
                             lambda.m=lambda.m,lambda.y=lambda.y)
      
      coef.alpha[b,]<-re.tmp$M$coefficients["Z",]
      coef.gamma[b,]<-re.tmp$Y$coefficients["Z",]
      coef.beta[b,]<-re.tmp$Y$coefficients["M",]
      
      coef.IE[b,]<-re.tmp$IE$coefficients
      
      c.alpha[b,]<-re.tmp$M$curve["Z",]
      c.gamma[b,]<-re.tmp$Y$curve["Z",]
      c.beta[b,]<-re.tmp$Y$curve["M",]
      c.IE[b,]<-re.tmp$IE$curve
      
      if(verbose)
      {
        print(paste0("Bootstrap sample ",b))
      }
    }
    se.alpha<-apply(coef.alpha,2,sd,na.rm=TRUE)
    se.gamma<-apply(coef.gamma,2,sd,na.rm=TRUE)
    se.beta<-apply(coef.beta,2,sd,na.rm=TRUE)
    se.IE<-apply(coef.IE,2,sd,na.rm=TRUE)
    
    se.c.alpha<-apply(c.alpha,2,sd,na.rm=TRUE)
    se.c.gamma<-apply(c.gamma,2,sd,na.rm=TRUE)
    se.c.beta<-apply(c.beta,2,sd,na.rm=TRUE)
    se.c.IE<-apply(c.IE,2,sd,na.rm=TRUE)
    
    if(boot.ci.type[1]=="bca")
    {
      ci.alpha<-apply(coef.alpha,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.gamma<-apply(coef.gamma,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.beta<-apply(coef.beta,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.IE<-apply(coef.IE,2,BC.CI,sims=sims,conf.level=conf.level)
      
      ci.c.alpha<-apply(c.alpha,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.gamma<-apply(c.gamma,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.beta<-apply(c.beta,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.IE<-apply(c.IE,2,BC.CI,sims=sims,conf.level=conf.level)
    }
    if(boot.ci.type[1]=="perc")
    {
      ci.alpha<-apply(coef.alpha,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.gamma<-apply(coef.gamma,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.beta<-apply(coef.beta,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.IE<-apply(coef.IE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      
      ci.c.alpha<-apply(c.alpha,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.gamma<-apply(c.gamma,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.beta<-apply(c.beta,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.IE<-apply(c.IE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
    }
    
    re.alpha<-rbind(apply(coef.alpha,2,mean,na.rm=TRUE),se.alpha,ci.alpha)
    re.gamma<-rbind(apply(coef.gamma,2,mean,na.rm=TRUE),se.gamma,ci.gamma)
    re.beta<-rbind(apply(coef.beta,2,mean,na.rm=TRUE),se.beta,ci.beta)
    re.IE<-rbind(apply(coef.IE,2,mean,na.rm=TRUE),se.IE,ci.IE)
    rownames(re.alpha)=rownames(re.gamma)=rownames(re.beta)=rownames(re.IE)<-c("Estimate","SE","LB","UB")
    colnames(re.alpha)=colnames(re.gamma)=colnames(re.beta)=colnames(re.IE)<-paste0("basis",1:ncol(fit.m$basis))
    
    curve.alpha<-rbind(apply(c.alpha,2,mean,na.rm=TRUE),se.c.alpha,ci.c.alpha)
    curve.gamma<-rbind(apply(c.gamma,2,mean,na.rm=TRUE),se.c.gamma,ci.c.gamma)
    curve.beta<-rbind(apply(c.beta,2,mean,na.rm=TRUE),se.c.beta,ci.c.beta)
    curve.IE<-rbind(apply(c.IE,2,mean,na.rm=TRUE),se.c.IE,ci.c.IE)
    rownames(curve.alpha)=rownames(curve.gamma)=rownames(curve.beta)=rownames(curve.IE)<-c("Estimate","SE","LB","UB")
    
    re<-list(alpha=list(coefficients=re.alpha,curve=curve.alpha),gamma=list(coefficients=re.gamma,curve=curve.gamma),beta=list(coefficients=re.beta,curve=curve.beta),
             IE=list(coefficients=re.IE,curve=curve.IE),DE=list(coefficients=re.gamma,curve=curve.gamma))
    
    return(re)
  }else
  {
    return(FMA.concurrent(Z,M,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda.m=lambda.m,lambda.y=lambda.y))
  }
}

### ---- from cfma/cfma-master/R/FMA.concurrent.CV.R ----
FMA.concurrent.CV <-
function(Z,M,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,
                            lambda=NULL,nfolds=5)
{
  N<-nrow(Z)             # # of subject
  ntp<-ncol(Z)           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }

  if(is.null(lambda))
  {
    lambda<-10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1])
  }
  
  # M model
  fit.m<-FDA.concurrent.CV(Z,M,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda,nfolds=nfolds,verbose=FALSE)
  lambda.m<-fit.m$lambda
  # Y model
  Xtmp<-array(NA,c(N,ntp,2))
  Xtmp[,,1]<-Z
  Xtmp[,,2]<-M
  fit.y<-FDA.concurrent.CV(Xtmp,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda,nfolds=nfolds,verbose=FALSE)
  lambda.y<-fit.y$lambda
  
  if(intercept)
  {
    coef.inter.m<-fit.m$coefficients[1:ncol(fit.m$basis)]
    coef.inter.y<-fit.y$coefficients[1:ncol(fit.m$basis)]
    curve.inter.m<-fit.m$gamma.curve[1,]
    curve.inter.y<-fit.y$gamma.curve[2,]
    
    coef.alpha<-fit.m$coefficients[-(1:ncol(fit.m$basis))]
    coef.gamma<-fit.y$coefficients[(ncol(fit.y$basis)+1):(2*ncol(fit.y$basis))]
    coef.beta<-fit.y$coefficients[(2*ncol(fit.y$basis)+1):(3*ncol(fit.y$basis))]
    
    curve.alpha<-fit.m$gamma.curve[2,]
    curve.gamma<-fit.y$gamma.curve[2,]
    curve.beta<-fit.y$gamma.curve[3,]
    
    coef.IE<-matrix(coef.alpha*coef.beta,nrow=1)
    rownames(coef.IE)<-"IE"
    colnames(coef.IE)<-paste0("basis",1:ncol(fit.m$basis))
    curve.IE<-curve.alpha*curve.beta
    
    # M model
    coef.m<-rbind(coef.inter.m,coef.alpha)
    rownames(coef.m)<-c("Intercept","Z")
    colnames(coef.m)<-paste0("basis",1:ncol(fit.m$basis))
    curve.m<-rbind(curve.inter.m,curve.alpha)
    rownames(curve.m)<-c("Intercept","Z")
    
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda=lambda.m)
    
    # Y model
    coef.y<-rbind(coef.inter.y,coef.gamma,coef.beta)
    rownames(coef.y)<-c("Intercept","Z","M")
    colnames(coef.y)<-paste0("basis",1:ncol(fit.y$basis))
    curve.y<-rbind(curve.inter.y,curve.gamma,curve.beta)
    rownames(curve.y)<-c("Intercept","Z","M")
    
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda=lambda.y)
    
    re.IE<-list(coefficients=coef.IE,curve=curve.IE)
    
    re.DE<-list(coefficients=coef.gamma,curve=curve.gamma)
    
    re<-list(basis=basis,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }else
  {
    coef.alpha<-fit.m$coefficients
    coef.gamma<-fit.y$coefficients[1:ncol(fit.y$basis)]
    coef.beta<-fit.y$coefficients[(ncol(fit.y$basis)+1):(2*ncol(fit.y$basis))]
    
    curve.alpha<-fit.m$gamma.curve[1,]
    curve.gamma<-fit.y$gamma.curve[1,]
    curve.beta<-fit.y$gamma.curve[2,]
    
    coef.IE<-matrix(coef.alpha*coef.beta,nrow=1)
    rownames(coef.IE)<-"IE"
    colnames(coef.IE)<-paste0("basis",1:ncol(fit.m$basis))
    curve.IE<-curve.alpha*curve.beta
    
    # M model
    coef.m<-rbind(coef.alpha)
    rownames(coef.m)<-c("Z")
    colnames(coef.m)<-paste0("basis",1:ncol(fit.m$basis))
    curve.m<-matrix(curve.alpha,nrow=1)
    rownames(curve.m)<-"Z"
    
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda=lambda.m)
    
    # Y model
    coef.y<-rbind(coef.gamma,coef.beta)
    rownames(coef.y)<-c("Z","M")
    colnames(coef.y)<-paste0("basis",1:ncol(fit.y$basis))
    curve.y<-rbind(curve.gamma,curve.beta)
    rownames(curve.y)<-c("Z","M")
    
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda=lambda.y)
    
    re.IE<-list(coefficients=coef.IE,curve=curve.IE)
    
    re.DE<-list(coefficients=coef.gamma,curve=curve.gamma)
    
    re<-list(basis=fit.m$basis,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }
  
  return(re)
}

### ---- from cfma/cfma-master/R/FMA.concurrent.R ----
FMA.concurrent <-
function(Z,M,Y,intercept=TRUE,basis=NULL,Ld2.basis=NULL,basis.type=c("fourier"),nbasis=3,timeinv=c(0,1),timegrids=NULL,lambda.m=0.01,lambda.y=0.01)
{
  N<-nrow(Z)             # # of subject
  ntp<-ncol(Z)           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # M model
  fit.m<-FDA.concurrent(Z,M,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda.m)
  # Y model
  Xtmp<-array(NA,c(N,ntp,2))
  Xtmp[,,1]<-Z
  Xtmp[,,2]<-M
  fit.y<-FDA.concurrent(Xtmp,Y,intercept=intercept,basis=basis,Ld2.basis=Ld2.basis,basis.type=basis.type,nbasis=nbasis,timeinv=timeinv,timegrids=timegrids,lambda=lambda.y)
  
  if(intercept)
  {
    coef.inter.m<-fit.m$coefficients[1:ncol(fit.m$basis)]
    coef.inter.y<-fit.y$coefficients[1:ncol(fit.m$basis)]
    curve.inter.m<-fit.m$gamma.curve[1,]
    curve.inter.y<-fit.y$gamma.curve[2,]
    
    coef.alpha<-fit.m$coefficients[-(1:ncol(fit.m$basis))]
    coef.gamma<-fit.y$coefficients[(ncol(fit.y$basis)+1):(2*ncol(fit.y$basis))]
    coef.beta<-fit.y$coefficients[(2*ncol(fit.y$basis)+1):(3*ncol(fit.y$basis))]
    
    curve.alpha<-fit.m$gamma.curve[2,]
    curve.gamma<-fit.y$gamma.curve[2,]
    curve.beta<-fit.y$gamma.curve[3,]
    
    coef.IE<-matrix(coef.alpha*coef.beta,nrow=1)
    rownames(coef.IE)<-"IE"
    colnames(coef.IE)<-paste0("basis",1:ncol(fit.m$basis))
    curve.IE<-curve.alpha*curve.beta
    
    # M model
    coef.m<-rbind(coef.inter.m,coef.alpha)
    rownames(coef.m)<-c("Intercept","Z")
    colnames(coef.m)<-paste0("basis",1:ncol(fit.m$basis))
    curve.m<-rbind(curve.inter.m,curve.alpha)
    rownames(curve.m)<-c("Intercept","Z")
    
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda=lambda.m)
    
    # Y model
    coef.y<-rbind(coef.inter.y,coef.gamma,coef.beta)
    rownames(coef.y)<-c("Intercept","Z","M")
    colnames(coef.y)<-paste0("basis",1:ncol(fit.y$basis))
    curve.y<-rbind(curve.inter.y,curve.gamma,curve.beta)
    rownames(curve.y)<-c("Intercept","Z","M")
    
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda=lambda.y)
    
    re.IE<-list(coefficients=coef.IE,curve=curve.IE)
    
    re.DE<-list(coefficients=coef.gamma,curve=curve.gamma)
    
    re<-list(basis=basis,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }else
  {
    coef.alpha<-fit.m$coefficients
    coef.gamma<-fit.y$coefficients[1:ncol(fit.y$basis)]
    coef.beta<-fit.y$coefficients[(ncol(fit.y$basis)+1):(2*ncol(fit.y$basis))]
    
    curve.alpha<-fit.m$gamma.curve[1,]
    curve.gamma<-fit.y$gamma.curve[1,]
    curve.beta<-fit.y$gamma.curve[2,]
    
    coef.IE<-matrix(coef.alpha*coef.beta,nrow=1)
    rownames(coef.IE)<-"IE"
    colnames(coef.IE)<-paste0("basis",1:ncol(fit.m$basis))
    curve.IE<-curve.alpha*curve.beta
    
    # M model
    coef.m<-rbind(coef.alpha)
    rownames(coef.m)<-c("Z")
    colnames(coef.m)<-paste0("basis",1:ncol(fit.m$basis))
    curve.m<-matrix(curve.alpha,nrow=1)
    rownames(curve.m)<-"Z"
    
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda=lambda.m)
    
    # Y model
    coef.y<-rbind(coef.gamma,coef.beta)
    rownames(coef.y)<-c("Z","M")
    colnames(coef.y)<-paste0("basis",1:ncol(fit.y$basis))
    curve.y<-rbind(curve.gamma,curve.beta)
    rownames(curve.y)<-c("Z","M")
    
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda=lambda.y)
    
    re.IE<-list(coefficients=coef.IE,curve=curve.IE)
    
    re.DE<-list(coefficients=coef.gamma,curve=curve.gamma)
    
    re<-list(basis=fit.m$basis,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }
  
  return(re)
}

### ---- from cfma/cfma-master/R/FMA.historical.boot.R ----
FMA.historical.boot <-
function(Z,M,Y,delta.grid1=1,delta.grid2=1,delta.grid3=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),
                               nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,lambda1.m=0.01,lambda2.m=0.01,lambda1.y=0.01,lambda2.y=0.01,
                               sims=1000,boot=TRUE,boot.ci.type=c("bca","perc"),conf.level=0.95,verbose=TRUE)
{
  # delta.grid1: M~Z time interval
  # delta.grid2: Y~Z time interval
  # delta.grid3: Y~M time interval
  
  N<-dim(Z)[1]             # # of subject
  ntp<-dim(Z)[2]           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # M model
  fit.m<-FDA.historical(Z,M,delta.grid=delta.grid1,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                        nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1.m,lambda2=lambda2.m)
  # Y model
  fit.y<-FDA.historical2(X1=Z,X2=M,Y,delta.grid1=delta.grid2,delta.grid2=delta.grid3,intercept=intercept,
                         basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                         nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1.y,lambda2=lambda2.y)
  
  if(boot)
  {
    coef.alpha=coef.gamma=coef.beta<-array(NA,c(ncol(fit.m$basis1),ncol(fit.m$basis2),sims))
    c.alpha=c.gamma=c.beta<-array(NA,c(ntp,ntp,sims))
    c.IE=c.DE<-matrix(NA,sims,ntp)
    for(b in 1:sims)
    {
      idx.tmp<-sample(1:N,N,replace=TRUE)
      
      Ztmp<-Z[idx.tmp,]
      Mtmp<-M[idx.tmp,]
      Ytmp<-Y[idx.tmp,]
      
      re.tmp<-FMA.historical(Ztmp,Mtmp,Ytmp,delta.grid1=delta.grid1,delta.grid2=delta.grid2,delta.grid3=delta.grid3,intercept=intercept,
                             basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,nbasis1=nbasis1,nbasis2=nbasis2,
                             timeinv=timeinv,timegrids=timegrids,lambda1.m=lambda1.m,lambda2.m=lambda2.m,lambda1.y=lambda1.y,lambda2.y=lambda2.y)
      
      coef.alpha[,,b]<-re.tmp$M$coefficients$alpha
      coef.gamma[,,b]<-re.tmp$Y$coefficients$gamma
      coef.beta[,,b]<-re.tmp$Y$coefficients$beta
      
      c.alpha[,,b]<-re.tmp$M$curve$alpha
      c.gamma[,,b]<-re.tmp$Y$curve$gamma
      c.beta[,,b]<-re.tmp$Y$curve$beta
      
      c.IE[b,]<-re.tmp$IE$curve
      c.DE[b,]<-re.tmp$DE$curve
      
      if(verbose)
      {
        print(paste0("Bootstrap sample ",b))
      }
    }
    se.alpha<-apply(coef.alpha,c(1,2),sd,na.rm=TRUE)
    se.gamma<-apply(coef.gamma,c(1,2),sd,na.rm=TRUE)
    se.beta<-apply(coef.beta,c(1,2),sd,na.rm=TRUE)
    
    se.c.alpha<-apply(c.alpha,c(1,2),sd,na.rm=TRUE)
    se.c.gamma<-apply(c.gamma,c(1,2),sd,na.rm=TRUE)
    se.c.beta<-apply(c.beta,c(1,2),sd,na.rm=TRUE)
    se.c.IE<-apply(c.IE,2,sd,na.rm=TRUE)
    se.c.DE<-apply(c.DE,2,sd,na.rm=TRUE)
    
    if(boot.ci.type[1]=="bca")
    {
      ci.alpha<-apply(coef.alpha,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      ci.gamma<-apply(coef.gamma,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      ci.beta<-apply(coef.beta,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      
      ci.c.alpha<-apply(c.alpha,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      ci.c.gamma<-apply(c.gamma,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      ci.c.beta<-apply(c.beta,c(1,2),BC.CI,sims=sims,conf.level=conf.level)
      ci.c.IE<-apply(c.IE,2,BC.CI,sims=sims,conf.level=conf.level)
      ci.c.DE<-apply(c.DE,2,BC.CI,sims=sims,conf.level=conf.level)
    }
    if(boot.ci.type[1]=="perc")
    {
      ci.alpha<-apply(coef.alpha,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.gamma<-apply(coef.gamma,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.beta<-apply(coef.beta,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      
      ci.c.alpha<-apply(c.alpha,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.gamma<-apply(c.gamma,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.beta<-apply(c.beta,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.IE<-apply(c.IE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      ci.c.DE<-apply(c.DE,2,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
    }
    
    re.alpha<-list(Estimate=apply(coef.alpha,c(1,2),mean,na.rm=TRUE),SE=se.alpha,LB=ci.alpha[1,,],UB=ci.alpha[2,,])
    re.gamma<-list(Estimate=apply(coef.gamma,c(1,2),mean,na.rm=TRUE),SE=se.gamma,LB=ci.gamma[1,,],UB=ci.gamma[2,,])
    re.beta<-list(Estimate=apply(coef.beta,c(1,2),mean,na.rm=TRUE),SE=se.beta,LB=ci.beta[1,,],UB=ci.beta[2,,])
    
    curve.alpha<-list(Estimate=apply(c.alpha,c(1,2),mean,na.rm=TRUE),SE=se.c.alpha,LB=ci.c.alpha[1,,],UB=ci.c.alpha[2,,])
    curve.gamma<-list(Estimate=apply(c.gamma,c(1,2),mean,na.rm=TRUE),SE=se.c.gamma,LB=ci.c.gamma[1,,],UB=ci.c.gamma[2,,])
    curve.beta<-list(Estimate=apply(c.beta,c(1,2),mean,na.rm=TRUE),SE=se.c.beta,LB=ci.c.beta[1,,],UB=ci.c.beta[2,,])
    
    curve.IE<-rbind(apply(c.IE,2,mean,na.rm=TRUE),se.c.IE,ci.c.IE)
    curve.DE<-rbind(apply(c.DE,2,mean,na.rm=TRUE),se.c.DE,ci.c.DE)
    rownames(curve.IE)=rownames(curve.DE)<-c("Estimate","SE","LB","UB")
    
    re<-list(alpha=list(coefficients=re.alpha,curve=curve.alpha),gamma=list(coefficients=re.gamma,curve=curve.gamma),
             beta=list(coefficients=re.beta,curve=curve.beta),IE=list(curve=curve.IE),DE=list(curve=curve.DE))
    
    return(re)
  }else
  {
    return(FMA.historical(Z,M,Y,delta.grid1=delta.grid1,delta.grid2=delta.grid2,delta.grid3=delta.grid3,intercept=intercept,
                          basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,nbasis1=nbasis1,nbasis2=nbasis2,
                          timeinv=timeinv,timegrids=timegrids,lambda1.m=lambda1.m,lambda2.m=lambda2.m,lambda1.y=lambda1.y,lambda2.y=lambda2.y))
  }
}

### ---- from cfma/cfma-master/R/FMA.historical.CV.R ----
FMA.historical.CV <-
function(Z,M,Y,delta.grid1=1,delta.grid2=1,delta.grid3=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),
                             nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,
                             lambda1=NULL,lambda2=NULL,nfolds=5)
{
  # delta.grid1: M~Z time interval
  # delta.grid2: Y~Z time interval
  # delta.grid3: Y~M time interval
  
  N<-dim(Z)[1]             # # of subject
  ntp<-dim(Z)[2]           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }

  if(is.null(lambda1))
  {
    lambda1<-10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1])
  }
  if(is.null(lambda2))
  {
    lambda2<-10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1])
  }
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # M model
  fit.m<-FDA.historical.CV(Z,M,delta.grid=delta.grid1,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                           nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1,lambda2=lambda2,nfolds=nfolds,verbose=FALSE)
  lambda1.m<-fit.m$lambda1
  lambda2.m<-fit.m$lambda2
  # Y model
  fit.y<-FDA.historical2.CV(X1=Z,X2=M,Y,delta.grid1=delta.grid2,delta.grid2=delta.grid3,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                            nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1,lambda2=lambda2,nfolds=nfolds,verbose=FALSE)
  lambda1.y<-fit.y$lambda1
  lambda2.y<-fit.y$lambda2
  
  if(intercept)
  {
    # M model
    coef.inter.m<-fit.m$coefficients[1:nbasis1,1:nbasis2]
    curve.inter.m<-fit.m$gamma.curve[,,1]
    
    coef.alpha<-fit.m$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.alpha<-fit.m$gamma.curve[,,2]
    
    # Y model
    coef.inter.y<-fit.y$coefficients[1:nbasis1,1:nbasis2]
    curve.inter.y<-fit.y$gamma.curve[,,1]
    
    coef.gamma<-fit.y$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.gamma<-fit.y$gamma.curve[,,2]
    
    coef.beta<-fit.y$coefficients[(2*nbasis1+1):(3*nbasis1),(2*nbasis2+1):(3*nbasis2)]
    curve.beta<-fit.y$gamma.curve[,,3]
    
    # IE and DE
    curve.IE<-rep(NA,ntp)
    curve.DE<-rep(NA,ntp)
    alpha.int<-rep(NA,ntp)
    for(i in 1:ntp)
    {
      rtmp1<-max(i-delta.grid1,1)
      rtmp2<-max(i-delta.grid2,1)
      rtmp3<-max(i-delta.grid3,1)
      
      alpha.int[i]<-int.func(curve.alpha[rtmp1:i,i],timeinv=c(timegrids[rtmp1],timegrids[i]),timegrids=timegrids[rtmp1:i])
      
      curve.IE[i]<-int.func(alpha.int[rtmp3:i]*curve.beta[rtmp3:i,i],timeinv=c(timegrids[rtmp3],timegrids[i]),timegrids=timegrids[rtmp3:i])
      
      curve.DE[i]<-int.func(curve.gamma[rtmp2:i,i],timeinv=c(timegrids[rtmp2],timegrids[i]),timegrids=timegrids[rtmp2:i])
    }
    
    coef.m<-list(Intercept=coef.inter.m,alpha=coef.alpha)
    curve.m<-list(Intercept=curve.inter.m,alpha=curve.alpha)
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda1=lambda1.m,lambda2=lambda2.m)
    
    coef.y<-list(Intercept=coef.inter.y,gamma=coef.gamma,beta=coef.beta)
    curve.y<-list(Intercept=curve.inter.y,gamma=curve.gamma,beta=curve.beta)
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda1=lambda1.y,lambda2=lambda2.y)
    
    re.IE<-list(curve=curve.IE)
    re.DE<-list(curve=curve.DE)
    
    re<-list(basis1=basis1,basis2=basis2,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }else
  {
    # M model
    coef.alpha<-fit.m$coefficients[1:nbasis1,1:nbasis2]
    curve.alpha<-fit.m$gamma.curve[,,1]
    
    # Y model
    coef.gamma<-fit.y$coefficients[1:nbasis1,1:nbasis2]
    curve.gamma<-fit.y$gamma.curve[,,1]
    
    coef.beta<-fit.y$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.beta<-fit.y$gamma.curve[,,2]
    
    # IE and DE
    curve.IE<-rep(NA,ntp)
    curve.DE<-rep(NA,ntp)
    alpha.int<-rep(NA,ntp)
    for(i in 1:ntp)
    {
      rtmp1<-max(i-delta.grid1,1)
      rtmp2<-max(i-delta.grid2,1)
      rtmp3<-max(i-delta.grid3,1)
      
      alpha.int[i]<-int.func(curve.alpha[rtmp1:i,i],timeinv=c(timegrids[rtmp1],timegrids[i]),timegrids=timegrids[rtmp1:i])
      
      curve.IE[i]<-int.func(alpha.int[rtmp3:i]*curve.beta[rtmp3:i,i],timeinv=c(timegrids[rtmp3],timegrids[i]),timegrids=timegrids[rtmp3:i])
      
      curve.DE[i]<-int.func(curve.gamma[rtmp2:i,i],timeinv=c(timegrids[rtmp2],timegrids[i]),timegrids=timegrids[rtmp2:i])
    }
    
    coef.m<-list(alpha=coef.alpha)
    curve.m<-list(alpha=curve.alpha)
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda1=lambda1.m,lambda2=lambda2.m)
    
    coef.y<-list(gamma=coef.gamma,beta=coef.beta)
    curve.y<-list(gamma=curve.gamma,beta=curve.beta)
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda1=lambda1.y,lambda2=lambda2.y)
    
    re.IE<-list(curve=curve.IE)
    re.DE<-list(curve=curve.DE)
    
    re<-list(basis1=basis1,basis2=basis2,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }
  
  return(re)
}

### ---- from cfma/cfma-master/R/FMA.historical.R ----
FMA.historical <-
function(Z,M,Y,delta.grid1=1,delta.grid2=1,delta.grid3=1,intercept=TRUE,basis1=NULL,Ld2.basis1=NULL,basis2=NULL,Ld2.basis2=NULL,basis.type=c("fourier"),
                          nbasis1=3,nbasis2=3,timeinv=c(0,1),timegrids=NULL,lambda1.m=0.01,lambda2.m=0.01,lambda1.y=0.01,lambda2.y=0.01)
{
  # delta.grid1: M~Z time interval
  # delta.grid2: Y~Z time interval
  # delta.grid3: Y~M time interval
  
  N<-dim(Z)[1]             # # of subject
  ntp<-dim(Z)[2]           # # of time points
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  # basis functions
  if(is.null(basis1))
  {
    if(basis.type[1]=="fourier")
    {
      basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
      
      Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
    }
  }else
  {
    nbasis1<-ncol(basis1)
  }
  if(is.null(basis2))
  {
    if(basis.type[1]=="fourier")
    {
      basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
      
      Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)
    }
  }else
  {
    nbasis2<-ncol(basis2)
  }
  
  # M model
  fit.m<-FDA.historical(Z,M,delta.grid=delta.grid1,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                        nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1.m,lambda2=lambda2.m)
  # Y model
  fit.y<-FDA.historical2(X1=Z,X2=M,Y,delta.grid1=delta.grid2,delta.grid2=delta.grid3,intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,basis.type=basis.type,
                         nbasis1=nbasis1,nbasis2=nbasis2,timeinv=timeinv,timegrids=timegrids,lambda1=lambda1.y,lambda2=lambda2.y)
  
  if(intercept)
  {
    # M model
    coef.inter.m<-fit.m$coefficients[1:nbasis1,1:nbasis2]
    curve.inter.m<-fit.m$gamma.curve[,,1]
    
    coef.alpha<-fit.m$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.alpha<-fit.m$gamma.curve[,,2]
    
    # Y model
    coef.inter.y<-fit.y$coefficients[1:nbasis1,1:nbasis2]
    curve.inter.y<-fit.y$gamma.curve[,,1]
    
    coef.gamma<-fit.y$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.gamma<-fit.y$gamma.curve[,,2]
    
    coef.beta<-fit.y$coefficients[(2*nbasis1+1):(3*nbasis1),(2*nbasis2+1):(3*nbasis2)]
    curve.beta<-fit.y$gamma.curve[,,3]
    
    # IE and DE
    curve.IE<-rep(NA,ntp)
    curve.DE<-rep(NA,ntp)
    alpha.int<-rep(NA,ntp)
    for(i in 1:ntp)
    {
      rtmp1<-max(i-delta.grid1,1)
      rtmp2<-max(i-delta.grid2,1)
      rtmp3<-max(i-delta.grid3,1)
      
      alpha.int[i]<-int.func(curve.alpha[rtmp1:i,i],timeinv=c(timegrids[rtmp1],timegrids[i]),timegrids=timegrids[rtmp1:i])
      
      curve.IE[i]<-int.func(alpha.int[rtmp3:i]*curve.beta[rtmp3:i,i],timeinv=c(timegrids[rtmp3],timegrids[i]),timegrids=timegrids[rtmp3:i])
      
      curve.DE[i]<-int.func(curve.gamma[rtmp2:i,i],timeinv=c(timegrids[rtmp2],timegrids[i]),timegrids=timegrids[rtmp2:i])
    }
    
    coef.m<-list(Intercept=coef.inter.m,alpha=coef.alpha)
    curve.m<-list(Intercept=curve.inter.m,alpha=curve.alpha)
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda1=lambda1.m,lambda2=lambda2.m)
    
    coef.y<-list(Intercept=coef.inter.y,gamma=coef.gamma,beta=coef.beta)
    curve.y<-list(Intercept=curve.inter.y,gamma=curve.gamma,beta=curve.beta)
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda1=lambda1.y,lambda2=lambda2.y)
    
    re.IE<-list(curve=curve.IE)
    re.DE<-list(curve=curve.DE)
    
    re<-list(basis1=basis1,basis2=basis2,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }else
  {
    # M model
    coef.alpha<-fit.m$coefficients[1:nbasis1,1:nbasis2]
    curve.alpha<-fit.m$gamma.curve[,,1]
    
    # Y model
    coef.gamma<-fit.y$coefficients[1:nbasis1,1:nbasis2]
    curve.gamma<-fit.y$gamma.curve[,,1]
    
    coef.beta<-fit.y$coefficients[(nbasis1+1):(2*nbasis1),(nbasis2+1):(2*nbasis2)]
    curve.beta<-fit.y$gamma.curve[,,2]
    
    # IE and DE
    curve.IE<-rep(NA,ntp)
    curve.DE<-rep(NA,ntp)
    alpha.int<-rep(NA,ntp)
    for(i in 1:ntp)
    {
      rtmp1<-max(i-delta.grid1,1)
      rtmp2<-max(i-delta.grid2,1)
      rtmp3<-max(i-delta.grid3,1)
      
      alpha.int[i]<-int.func(curve.alpha[rtmp1:i,i],timeinv=c(timegrids[rtmp1],timegrids[i]),timegrids=timegrids[rtmp1:i])
      
      curve.IE[i]<-int.func(alpha.int[rtmp3:i]*curve.beta[rtmp3:i,i],timeinv=c(timegrids[rtmp3],timegrids[i]),timegrids=timegrids[rtmp3:i])
      
      curve.DE[i]<-int.func(curve.gamma[rtmp2:i,i],timeinv=c(timegrids[rtmp2],timegrids[i]),timegrids=timegrids[rtmp2:i])
    }
    
    coef.m<-list(alpha=coef.alpha)
    curve.m<-list(alpha=curve.alpha)
    re.m<-list(coefficients=coef.m,curve=curve.m,fitted=fit.m$fitted,lambda1=lambda1.m,lambda2=lambda2.m)
    
    coef.y<-list(gamma=coef.gamma,beta=coef.beta)
    curve.y<-list(gamma=curve.gamma,beta=curve.beta)
    re.y<-list(coefficients=coef.y,curve=curve.y,fitted=fit.y$fitted,lambda1=lambda1.y,lambda2=lambda2.y)
    
    re.IE<-list(curve=curve.IE)
    re.DE<-list(curve=curve.DE)
    
    re<-list(basis1=basis1,basis2=basis2,M=re.m,Y=re.y,IE=re.IE,DE=re.DE)
  }
  
  return(re)
}

### ---- from cfma/cfma-master/R/fourier.basis.R ----
fourier.basis <-
function(timeinv=c(0,1),ntp,nbasis=3)
{
  if(nbasis%%2==0)
  {
    nbasis<-nbasis-1
  }
  
  timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  
  r<-timeinv[2]-timeinv[1]
  
  basis<-matrix(NA,ntp,nbasis)
  basis[,1]<-rep(1/sqrt(r),ntp)
  for(j in 1:floor(nbasis/2))
  {
    basis[,j*2]<-sin(j*(2*pi*(timegrids-timeinv[1])/r))*sqrt(2/r)
    basis[,j*2+1]<-cos(j*(2*pi*(timegrids-timeinv[1])/r))*sqrt(2/r)
  }
  
  return(basis)
}

### ---- from cfma/cfma-master/R/int.func.R ----
int.func <-
function(x,timeinv=c(0,1),timegrids=NULL)
{
  ntp<-length(x)
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  if(timeinv[2]<=timeinv[1])
  {
    int<-0
  }else
  {
    if(length(timegrids)==ntp)
    {
      int<-0
      for(i in 1:(ntp-1))
      {
        int<-int+(x[i]+x[i+1])*(timegrids[i+1]-timegrids[i])/2
      }
    }else
    {
      stop("Error!")
    }
  }
  
  return(int)
}

### ---- from cfma/cfma-master/R/Ld2.fourier.R ----
Ld2.fourier <-
function(timeinv=c(0,1),ntp,nbasis=3)
{
  if(nbasis%%2==0)
  {
    nbasis<-nbasis-1
  }
  
  timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  
  r<-timeinv[2]-timeinv[1]
  
  db<-matrix(NA,ntp,nbasis)
  db[,1]<-rep(0,ntp)
  for(j in 1:floor(nbasis/2))
  {
    db[,j*2]<--sqrt(2/r)*(2*pi*j/r)^2*sin(2*pi*j*(timegrids-timeinv[1])/r)
    db[,j*2+1]<--sqrt(2/r)*(2*pi*j/r)^2*cos(2*pi*j*(timegrids-timeinv[1])/r)
  }
  
  return(db)
}

