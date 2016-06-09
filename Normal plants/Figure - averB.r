
library(plotBy)
source("Normal plants/Bootstrap predict_nls.R")
source("Normal plants/Bootstrap functions.R")
data <- read.csv("Normal plants/Derived variables.csv")
dataNPP <- read.csv("Normal plants/NPP.csv")
colnames(dataNPP) <- c("MAP", "TNPP")
MAP <- dataNPP$MAP
TNPP <- dataNPP$TNPP

nls1 <- nls(TNPP ~ a*MAP^b/exp(c*MAP), start=list(a = 0.5346608789, b = 1.0787801564, c= 0.0003562646), data=dataNPP)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

# average B
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 1, 4), mfrow=c(1,1))
p1 <- data.frame(x = p$x, pred = p$pred, lwr = p$lwr, upr = p$upr)
p2 <- subset(p1, x<=2500)

plot(MAP, TNPP, panel.first = {
  addpoly(p2$x, p2$lwr, p2$upr)
  lines(p2$x, p2$pred)
},
type = "n",
xlim = c(0, 3000),ylim = c(0, 1420),
xlab = "", ylab = "",
axes = FALSE)

axis(4, ylim = c(0, 1200), pos = 3000, lwd = 2, at = c(0, 500, 1000))
mtext(expression(TNPP~(gC~m^-2~yr^-1)),side = 4,line = 2.8, cex = 1.3)
rect(0, 0, 100, 1420, col = "white", border = NA)

par(new=TRUE)

Cols <- c("black", "red", "blue")

plotBy(data$averB ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 15),
       xlab="MAP (mm)", ylab=NA,
       xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(bar(italic(B))~(mu*mol~m^-2~s^-1)), side=2, line=1.7, cex=1.3)
legend("topleft", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))
