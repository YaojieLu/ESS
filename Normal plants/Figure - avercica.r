
library(plotBy)
source("Normal plants/Bootstrap predict_nls.R")
source("Normal plants/Bootstrap functions.R")
data <- read.csv("Normal plants/Derived variables.csv")
datadelta13 <- read.csv("Normal plants/delta13.csv")
colnames(datadelta13) <- c("MAP", "delta13")
MAP <- datadelta13$MAP

cica <- (-8.3-datadelta13$delta13-4.4)/22.6
nls1 <- nls(cica ~ a+b*MAP, start=list(a=0.411991, b=0.000619469), data=datadelta13)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,1))
Cols <- c("red","darkgreen","blue")
plot(MAP, cica, panel.first={
  addpoly(p$x, p$lwr, p$upr)
  lines(p$x, p$pred)},
  type="n", xlim=c(0, 3000), ylim=c(0, 1),# xaxt="n",
  xlab=NA, ylab=expression(italic(c[i]/c[a])), cex.lab=1.3)

plotBy(data$avercica ~ MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 1),
       xlab="MAP (mm)", ylab=expression(italic(c[i]/c[a])),
       cex.lab=1.3, col=Cols, add=TRUE)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
mtext("MAP (mm)", side=1, line=2.2, cex=1.3)
legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Prentice~italic(et~al.)~2011),
       col=c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd=c(2, 2, 2, 2))
