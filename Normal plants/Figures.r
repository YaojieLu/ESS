
library(plotBy)
data <- read.csv("Normal plants/Derived variables.csv")

# Figures
# average B
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
Cols <- c("black", "red", "blue")
plotBy(data$averB ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 15),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(bar(italic(B))~(mu*mol~m^-2~s^-1)), side=2, line=1.7, cex=1.3)
legend("topleft", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))

# average w
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
Cols <- c("black", "red", "blue")
plotBy(data$averw ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 1),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(bar(italic(w))), side=2, line=1.7, cex=1.3)
legend("topleft", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))

# average fL
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
Cols <- c("black", "red", "blue")
plotBy(data$f0 ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 1),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(f[L])), side=2, line=1.7, cex=1.3)
legend("topright", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))

# average E
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
Cols <- c("black", "red", "blue")
plotBy(data$averE*500 ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 8),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 8), pos=0, lwd=2)
mtext(expression(italic(bar(E))~"(mm per day)"), side=2, line=1.7, cex=1.3)
legend("topleft", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))

# average E/MAP
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
Cols <- c("black", "red", "blue")
plotBy(data$averE*500*365/data$MAP ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 1),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(bar(E))/MAP), side=2, line=1.7, cex=1.3)
legend("bottomleft", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))
