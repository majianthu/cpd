library(copent)
library(SMPracticals)

data("coal")
h1 = hist(floor(coal$date), breaks = 1850:1962)
disasters = h1$counts
years = h1$breaks[2:113]

cp1 = cpd(disasters)

x11();
par(mfrow = c(2,1))
plot(years,disasters, type = "l", main = "British Coal Mine Disasters")
lines(x=c(years[1],years[cp1$pos-1]),y=rep(mean(disasters[1:(cp1$pos-1)]),2),col='red')
lines(x=c(years[cp1$pos],years[length(disasters)]),y=rep(mean(disasters[cp1$pos:length(disasters)]),2),col='red')
plot(years,c(0,cp1$stats,0), type = "l", xlab = "Time", ylab = "statistic", main = "Statistics")
lines(x = c(years[cp1$pos],years[cp1$pos]), y = c(0, max(disasters)), col = 'red')

