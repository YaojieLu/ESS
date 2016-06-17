
options(digits=2)
library(plotBy)
source("Normal plants/Bootstrap predict_nls.R")
source("Normal plants/Bootstrap functions.R")
data <- read.csv("Normal plants/Derived variables.csv")
dataNPP <- read.csv("Normal plants/NPP.csv")
colnames(dataNPP) <- c("MAP", "TNPP")
MAP <- dataNPP$MAP
TNPP <- dataNPP$TNPP

nls1 <- nls(TNPP ~ a*MAP^b/exp(c*MAP), start=list(a=0.5346608789, b=1.0787801564, c= 0.0003562646), data=dataNPP)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")
p1 <- data.frame(x=p$x, pred=p$pred, lwr=p$lwr, upr=p$upr)
p2 <- subset(p1, x<=2500)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 3.5), mfrow=c(1,1))
plot(MAP, TNPP, panel.first={
  addpoly(p2$x, p2$lwr, p2$upr)
  lines(p2$x, p2$pred)},
  type="n", xlab=NA, ylab=NA, axes=FALSE,
  xlim=c(0, 3000), ylim=c(0, 1420))

axis(4, ylim=c(0, 1200), pos=3000, lwd=2, at=c(0, 500, 1000))
mtext(expression(TNPP~(gC~m^-2~yr^-1)),side=4,line=2.4, cex=1.3)
rect(0, 0, 100, 1420, col="white", border=NA)

par(new=TRUE)
# average A
Cols <- c("red","darkgreen","blue")
plotBy(data$averA ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 15),
       xlab="MAP (mm)", ylab=NA,
       xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols)

## average m
#plotBy(data$averm ~ data$MAP | k, data=data,
#       type='l', legend=FALSE, legendwhere="topleft",
#       xlim=c(0, 3000), ylim=c(0, 15),
#       xlab="MAP (mm)", ylab=NA,
#       xaxt="n", yaxt="n",
#       cex.lab=1.3, col=Cols, lty=2, add=T)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(bar(italic(A))~(mu*mol~m^-2~s^-1)), side=2, line=1.8, cex=1.3)#*" or "*bar(italic(m))
legend("topleft", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Del~Grosso~italic(et~al.)~2008),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd=c(2, 2, 2, 2))
