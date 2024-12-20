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
library(L2hdchange)
library(fpop)
library(HDCD)
library(HDDchangepoint) # devtools::install_github("rezadrikvandi/HDDchangepoint")
library(HDcpDetect)
library(decp)
library(cpss)

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

# performance metric
pm <- function(x,ref,thd = 5){
  r1 = 0
  if(length(x)>1){
  for(i in 1:length(ref)){
    for(j in 1:length(x))
      if(abs(x[j]-ref[i])<=thd){
        r1 = r1 + 1
        break
      }
  }
  }else if(length(x)==1 & is.numeric(x)){
    for(i in 1:length(ref)){
      if(abs(x-ref[i])<=thd){
        r1 = r1 + 1
        break
      }
    }
  }
  r1
}

## for simulation
k1 = 10
ref1 = c(50,100,150)

### univariate multiple change points
n1 = 50
# case 1: mean
p1 = rep(0,26)
for(k in 1:k1){
x=c(rnorm(n1,0,1),rnorm(n1,5,1),rnorm(n1,10,1),rnorm(n1,3,1))
mresult1 = copent::mcpd(x)
#mresult1
p1[1] = p1[1] + pm(mresult1$pos-1,ref1) / k1

cpt1 = cpt.mean(x, method = "BinSeg", Q = 5)
p1[2] = p1[2] + pm(cpt1@cpts,ref1) / k1
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
p1[4] = p1[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p1[5] = p1[5] + pm(npmojo1,ref1) / k1
npwbs1 = detectChanges(x)
p1[6] = p1[6] + pm(npwbs1,ref1) / k1
mft1 = MFT.mean(x)$CP[,1]
p1[7] = p1[7] + pm(mft1,ref1) / k1
jcp1 = jcp(x)$SFA[[1]]
p1[8] = p1[8] + pm(jcp1,ref1) / k1
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)$changepoints[,1]
p1[9] = p1[9] + pm(inspect1,ref1) / k1
cptnp1 = cpt.np(x)
p1[11] = p1[11] + pm(cptnp1@cpts,ref1) / k1
mosum1 = mosum(x,G=20)
p1[13] = p1[13] + pm(mosum1$cpts,ref1) / k1
snseg1 = SNSeg_Uni(x, paras_to_test = "mean")$est_cp
p1[14] = p1[14] + pm(snseg1,ref1) / k1
idpcm1 = ID_pcm(x)$cpt
p1[16] = p1[16] + pm(idpcm1,ref1) / k1
#sbs1 = changepoints(sbs(x))$cpt.th[[1]]
wbs1 = changepoints(wbs(x))$cpt.th[[1]]
p1[17] = p1[17] + pm(wbs1,ref1) / k1
breakfast1 = breakfast(x)$cptmodel.list
p1[18] = p1[18] + pm(breakfast1[[6]]$cpts,ref1) / k1
mscp1 = mscp(x)
p1[19] = p1[19] + pm(mscp1$cp,ref1) / k1
fpop1 = Fpop(x,n1*4)$t.est
p1[21] = p1[21] + pm(fpop1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p1[22] = p1[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
cpss1 = cpss.mean(x)@cps
p1[26] = p1[26] + pm(cpss1,ref1) / k1
}

# case 2: mean-var
p2 = rep(0,26)
for(k in 1:k1){
x=c(rnorm(n1,0,1),rnorm(n1,5,3),rnorm(n1,10,1),rnorm(n1,3,10))
mresult1 = mcpd(x)
mresult1
p2[1] = p2[1] + pm(mresult1$pos-1,ref1) / k1
cpt1 = cpt.meanvar(x, method = "BinSeg", Q = 5)
p2[2] = p2[2] + pm(cpt1@cpts,ref1) / k1
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
p2[4] = p2[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p2[5] = p2[5] + pm(npmojo1,ref1) / k1
npwbs1 = detectChanges(x)
p2[6] = p2[6] + pm(npwbs1,ref1) / k1
jcp1 = jcp(x)$SFA[[1]]
p2[8] = p2[8] + pm(jcp1,ref1) / k1
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)
p2[9] = p2[9] + pm(inspect1$changepoints,ref1) / k1
cptnp1 = cpt.np(x)
p2[11] = p2[11] + pm(cptnp1@cpts,ref1) / k1
mosum1 = mosum(x,G=40)
p2[13] = p2[13] + pm(mosum1$cpts,ref1) / k1
snseg1 = SNSeg_Uni(x, paras_to_test = c("mean","variance"))$est_cp
p2[14] = p2[14] + pm(snseg1,ref1) / k1
fpop1 = Fpop(x,n1*4)$t.est
p2[21] = p2[21] + pm(fpop1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p2[22] = p2[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
cpss1 = cpss.meanvar(x)@cps
p2[26] = p2[26] + pm(cpss1,ref1) / k1
}

# case 3: var
p3 = rep(0,26)
for(k in 1:k1){
x=c(rnorm(n1,0,1),rnorm(n1,0,10),rnorm(n1,0,5),rnorm(n1,0,1))
mresult1 = mcpd(x,thd = 0.05)
p3[1] = p3[1] + pm(mresult1$pos-1,ref1) / k1
cpt1 = cpt.var(x, method = "BinSeg", Q = 5)
p3[2] = p3[2] + pm(cpt1@cpts,ref1) / k1
d = rid(x, M = 1000, tau = "clustering")
rid1 = localization(x,d$Good_Intervals)
p3[4] = p3[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p3[5] = p3[5] + pm(npmojo1,ref1) / k1
npwbs1 = detectChanges(x)
p3[6] = p3[6] + pm(npwbs1,ref1) / k1
jcp1 = jcp(x)$SFA[[1]]
p3[8] = p3[8] + pm(jcp1,ref1) / k1
thd1 = compute.threshold(200,1)
inspect1 = inspect(x,thd1)
p3[9] = p3[9] + pm(inspect1$changepoints,ref1) / k1
cptnp1 = cpt.np(x)
p3[11] = p3[11] + pm(cptnp1@cpts,ref1) / k1
mosum1 = mosum(x,G=50)
p3[13] = p3[13] + pm(mosum1$cpts,ref1) / k1
snseg1 = SNSeg_Uni(x, paras_to_test = "variance")$est_cp
p3[14] = p3[14] + pm(snseg1,ref1) / k1
fpop1 = Fpop(x,n1*4)$t.est
p3[21] = p3[21] + pm(fpop1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p3[22] = p3[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
cpss1 = cpss.var(x)@cps
p3[26] = p3[26] + pm(cpss1,ref1) / k1
}

### multivariate multiple change points
rho1 = 0.2; rho2 = 0.8; rho3 = 0.1; rho4 = 0.9
n1 = 50

# case 1: mean
p4 = rep(0,26)
for(k in 1:k1){
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(10,10),matrix(c(1,rho1,rho1,1),2,2))
x3 = rmnorm(n1,c(5,5),matrix(c(1,rho1,rho1,1),2,2))
x4 = rmnorm(n1,c(1,1),matrix(c(1,rho1,rho1,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x)
p4[1] = p4[1] + pm(mresult2$pos-1,ref1) / k1
kcp1 = kcpa(x,5,500)
p4[3] = p4[3] + pm(kcp1,ref1) / k1
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
p4[4] = p4[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p4[5] = p4[5] + pm(npmojo1,ref1) / k1
#hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
p4[24] = p4[24] + pm(hdcp2,ref1) / k1
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
p4[9] = p4[9] + pm(inspect1$changepoints,ref1) / k1
dcbs1 = dcbs.alg(t(x))
p4[10] = p4[10] + pm(dcbs1$ecp,ref1) / k1
#sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
p4[12] = p4[12] + pm(geomcp1@ang.cpts,ref1) / k1
snseg1 = SNSeg_Multi(x, paras_to_test = "mean")$est_cp
p4[14] = p4[14] + pm(snseg1,ref1) / k1
offline1 = ChangePoints(x)$change_point
p4[15] = p4[15] + pm(offline1,ref1) / k1
ts_no_nbd1 = ts_hdchange(t(x))
hdchange1 = hdchange(ts_no_nbd1)$time_stamps
p4[20] = p4[20] + pm(hdchange1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p4[22] = p4[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
# hdd1 = multiple_changepoint_detection(x,npermut = 10)$Detected
# p4[23] = p4[23] + pm(hdd1,ref1) / k1
decp1 = decp(x,num_simulations = 5000,num_iterations = 50)
p4[25] = p4[25] + pm(decp1$ordered_change_points,ref1) / k1
}

# case 2: mean-var
p5 = rep(0,26)
for(k in 1:k1){
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(10,10),matrix(c(1,rho2,rho2,1),2,2))
x3 = rmnorm(n1,c(5,5),matrix(c(1,rho3,rho3,1),2,2))
x4 = rmnorm(n1,c(1,1),matrix(c(1,rho4,rho4,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x)
p5[1] = p5[1] + pm(mresult2$pos-1,ref1) / k1
kcp1 = kcpa(x,5,500)
p5[3] = p5[3] + pm(kcp1,ref1) / k1
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
p5[4] = p5[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p5[5] = p5[5] + pm(npmojo1,ref1) / k1
#hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
p5[24] = p5[24] + pm(hdcp2,ref1) / k1
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
p5[9] = p5[9] + pm(inspect1$changepoints,ref1) / k1
dcbs1 = dcbs.alg(t(x))
p5[10] = p5[10] + pm(dcbs1$ecp,ref1) / k1
#sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
p5[12] = p5[12] + pm(geomcp1@ang.cpts,ref1) / k1
snseg1 = SNSeg_Multi(x, paras_to_test = "mean")$est_cp
p5[14] = p5[14] + pm(snseg1,ref1) / k1
offline1 = ChangePoints(x)$change_point
p5[15] = p5[15] + pm(offline1,ref1) / k1
ts_no_nbd1 = ts_hdchange(t(x))
hdchange1 = hdchange(ts_no_nbd1)$time_stamps
p5[20] = p5[20] + pm(hdchange1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p5[22] = p5[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
# hdd1 = multiple_changepoint_detection(x,npermut = 10)$Detected
# p5[23] = p5[23] + pm(hdd1,ref1) / k1
decp1 = decp(x,num_simulations = 5000,num_iterations = 50)
p5[25] = p5[25] + pm(decp1$ordered_change_points,ref1) / k1
}

# case 3: var
p6 = rep(0,26)
for(k in 1:k1){
x1 = rmnorm(n1,c(0,0),matrix(c(1,rho1,rho1,1),2,2))
x2 = rmnorm(n1,c(0,0),matrix(c(1,rho2,rho2,1),2,2))
x3 = rmnorm(n1,c(0,0),matrix(c(1,rho3,rho3,1),2,2))
x4 = rmnorm(n1,c(0,0),matrix(c(1,rho4,rho4,1),2,2))
x = rbind(x1,x2,x3,x4)
mresult2 = mcpd(x,thd = 0.05)
p6[1] = p6[1] + pm(mresult2$pos-1,ref1) / k1
kcp1 = kcpa(x,5,500)
p6[3] = p6[3] + pm(kcp1,ref1) / k1
d = rid(t(x), M = 1000, tau = "clustering")
rid1 = localization(x, d$Good_Intervals)
p6[4] = p6[4] + pm(rid1,ref1) / k1
npmojo1 = np.mojo.multilag(x,30)$cpts[,1]
p6[5] = p6[5] + pm(npmojo1,ref1) / k1
#hdcp1 = binary.segmentation(x)
hdcp2 = wild.binary.segmentation(x)
p6[24] = p6[24] + pm(hdcp2,ref1) / k1
thd1 = compute.threshold(200,2)
inspect1 = inspect(t(x),thd1)
p6[9] = p6[9] + pm(inspect1$changepoints,ref1) / k1
dcbs1 = dcbs.alg(t(x))
p6[10] = p6[10] + pm(dcbs1$ecp,ref1) / k1
#sbs1 = sbs.alg(t(x))
geomcp1 = geomcp(x)
p6[12] = p6[12] + pm(geomcp1@ang.cpts,ref1) / k1
snseg1 = SNSeg_Multi(x, paras_to_test = "mean")$est_cp
p6[14] = p6[14] + pm(snseg1,ref1) / k1
offline1 = ChangePoints(x)$change_point
p6[15] = p6[15] + pm(offline1,ref1) / k1
ts_no_nbd1 = ts_hdchange(t(x))
hdchange1 = hdchange(ts_no_nbd1)$time_stamps
p6[20] = p6[20] + pm(hdchange1,ref1) / k1
esac1 = ESAC(t(x))$changepoints
p6[22] = p6[22] + pm(esac1,ref1) / k1
#pilliat1 = Pilliat(t(x))$changepoints
# hdd1 = multiple_changepoint_detection(x,npermut = 10)$Detected
# p6[23] = p6[23] + pm(hdd1,ref1) / k1
decp1 = decp(x,num_simulations = 5000,num_iterations = 50)
p6[25] = p6[25] + pm(decp1$ordered_change_points,ref1) / k1
}

save(p1,p2,p3,p4,p5,p6, file = "~/Rworks/cpd/sim_results_p16.RData")

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

## plot performance measures of cpt algs

cptname1 = c("copent","changepoint","ecp","rid","CptNonPar","npwbs","MFT","jcp","InspectChangepoint","hdbinseg",
             "changepoint.np","changepoint.geo","mosum","SNSeg","offlineChange","IDetect","wbs","breakfast",
             "mscp","L2hdchange","fpop","HDCD","HDDchangepoint","HDcpDetect","decp","cpss")
col1 = rep('grey',26)
idx1 = c(1,2,4:9,11,13,14,16:19,21,22,26); # mean
idx1 = c(1,2,4:6,8,9,11,13,14,21,22,26) # mean-var, var
idx1 = c(1,3:5,9,10,12,14,15,20,22,24:25) # multivariate
col1[idx1] = 'red'

x11(width = 11, height = 2)
par(mar = c(3,4,1,1))
at1 = barplot(p5, xaxt = "n", ylim = c(-2,3), ylab = "points")
text(x = at1, y = -0.2, srt = 60, adj = 1, xpd = TRUE, labels = cptname1, col = col1, cex = 0.65)
