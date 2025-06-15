library(copent)

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
