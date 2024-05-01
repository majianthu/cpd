library(copent)
library(changepoint)
library(ecp)
library(mnormt)
library(copula)
library(rid)
library(CptNonPar)
library(npwbs)
library(MFT)
library(jcp)
library(InspectChangepoint)
library(hdbinseg)
library(changepoint.np)
library(changepoint.geo)
library(mosum)
library(SNSeg)
library(offlineChange)
library(IDetect)
library(wbs)
library(breakfast)
library(mscp)

cpd <- function(x,thd=0.13,n=30){
  result = {}
  x = as.matrix(x)
  len1 = dim(x)[1]
  stat1 = rep(0,len1-1)
  for(i in 2:(len1-2)){
    s0 = as.matrix(x[1:i,])
    s1 = as.matrix(x[(i+1):len1,])
    stat1[i] = tst(s0,s1,n)
  }
  if(max(stat1)>thd){
    result$stats = stat1
    result$maxstat = max(stat1)
    result$pos = which(stat1 == max(stat1))+1
  }
  
  return(result)
}

mcpd <- function(x,maxp=5,thd=0.13,minseglen=10,n=30){
  mresult = {}
  x = as.matrix(x)
  len1 = dim(x)[1]
  bisegs = matrix(c(1,len1-1),1,2)
  k = 1
  for(i in 1:maxp){
    if(i>dim(bisegs)[1]){ break } 
    ri = cpd(x[bisegs[i,1]:bisegs[i,2],],thd,n)
    if(!is.null(ri)){
      ri$pos = ri$pos + bisegs[i,1] - 1
      mresult$maxstat[k] = ri$maxstat
      mresult$pos[k] = ri$pos
      k = k + 1
      if((ri$pos-bisegs[i,1])>minseglen) { bisegs = rbind(bisegs,c(bisegs[i,1],ri$pos-1)) }
      if((bisegs[i,2]-ri$pos+1)>minseglen) { bisegs = rbind(bisegs,c(ri$pos,bisegs[i,2])) }
    }
  }
  return(mresult)
}

### real data

## Nile
data("Nile")
result1 = cpd(Nile)

x11();
par(mfrow = c(2,1))
plot(Nile, type = "l", main = "Nile")
lines(x=1870+c(1,result1$pos),y=rep(mean(Nile[1:result1$pos]),2),col='red')
lines(x=1870+c(result1$pos,length(Nile)),y=rep(mean(Nile[(result1$pos+1):length(Nile)]),2),col='red')
plot(result1$stats, type = "l", xaxt = "n", xlab = "Time", ylab = "statistic", main = "Statistics")
axis(1, at = seq(9,90,by = 20), labels = seq(1880,1960,by = 20))
lines(x=c(result1$pos,result1$pos)-1, y=c(0,result1$maxstat),col='red')

### univariate multiple change points
# case 1: mean
x=c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,10,1),rnorm(50,3,1))
mresult1 = mcpd(x)
cpt1 = cpt.mean(x, method = "BinSeg", Q = 5)
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
npwbs1 = detectChanges(x)
mft1 = MFT.mean(x)$CP[,1]
jcp1 = jcp(x)$SFA[[1]]
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)
cptnp1 = cpt.np(x)
mosum1 = mosum(x,G=20)
snseg1 = SNSeg_Uni(x, paras_to_test = "mean")$est_cp
idpcm1 = ID_pcm(x)$cpt
sbs1 = changepoints(sbs(x))$cpt.th[[1]]
wbs1 = changepoints(wbs(x))$cpt.th[[1]]
breakfast1 = breakfast(x)$cptmodel.list
mscp1 = mscp(x)

# case 2: mean-var
x=c(rnorm(50,0,1),rnorm(50,5,3),rnorm(50,10,1),rnorm(50,3,10))
mresult1 = mcpd(x)
cpt1 = cpt.meanvar(x, method = "BinSeg", Q = 5)
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
npwbs1 = detectChanges(x)
jcp1 = jcp(x)$SFA[[1]]
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)
cptnp1 = cpt.np(x)
mosum1 = mosum(x,G=40)
snseg1 = SNSeg_Uni(x, paras_to_test = c("mean","variance"))$est_cp

# case 3: var
x=c(rnorm(50,0,1),rnorm(50,0,10),rnorm(50,0,5),rnorm(50,0,1))
mresult1 = mcpd(x)
cpt1 = cpt.var(x, method = "BinSeg", Q = 5)
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
npwbs1 = detectChanges(x)
jcp1 = jcp(x)$SFA[[1]]
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)
cptnp1 = cpt.np(x)
mosum1 = mosum(x,G=50)
snseg1 = SNSeg_Uni(x, paras_to_test = "variance")$est_cp

### multivariate multiple change points
rho1 = 0.2; rho2 = 0.8; rho3 = 0.1; rho4 = 0.9
n1 = 50

# case 1: mean
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(10,10),matrix(c(1,rho1,rho1,1),2,2))
x3 = rmnorm(n1,c(5,5),matrix(c(1,rho1,rho1,1),2,2))
x4 = rmnorm(n1,c(1,0),matrix(c(1,rho1,rho1,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x)
kcp1 = kcpa(x,5,500)
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
dcbs1 = dcbs.alg(t(x))
sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
snseg1 = SNSeg_Multi(x, paras_to_test = "mean")$est_cp
offline1 = ChangePoints(x)$change_point

# case 2: mean-var
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(10,10),matrix(c(1,rho2,rho2,1),2,2))
x3 = rmnorm(n1,c(5,5),matrix(c(1,rho3,rho3,1),2,2))
x4 = rmnorm(n1,c(1,0),matrix(c(1,rho4,rho4,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x)
kcp1 = kcpa(x,5,200)
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
dcbs1 = dcbs.alg(t(x))
sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
snseg1 = SNSeg_Multi(x, paras_to_test = "mean")$est_cp
offline1 = ChangePoints(x)$change_point

# case 3: var
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(0,0),matrix(c(1,rho2,rho2,1),2,2))
x3 = rmnorm(n1,c(0,0),matrix(c(1,rho3,rho3,1),2,2))
x4 = rmnorm(n1,c(0,0),matrix(c(1,rho4,rho4,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x,thd = 0.05)
kcp1 = kcpa(x,5,200)
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
dcbs1 = dcbs.alg(t(x))
sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
snseg1 = SNSeg_Multi(x, paras_to_test = "covariance")$est_cp
offline1 = ChangePoints(x)$change_point

# case 4: with non-normality
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
mv.cop <- mvdc(frankCopula(0.9), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 0.5)))
x2 <- rMvdc(n1, mv.cop)
x3 = rmnorm(n1,c(0,0),matrix(c(1,rho3,rho3,1),2,2))
mv.cop <- mvdc(normalCopula(0.3), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 0.5)))
x4 <- rMvdc(n1, mv.cop)
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x)
kcp1 = kcpa(x,5,500)
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
dcbs1 = dcbs.alg(t(x))
sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
