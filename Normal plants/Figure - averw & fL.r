
library(plotBy)
data <- read.csv("Normal plants/Derived variables.csv")

# Figures
# average w
windows(8, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, oma=c(2.5, 0, 0, 0), mar=c(1, 3.5, 1, 1), mfrow=c(2, 1))
Cols <- c("red", "darkgreen", "blue")
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

axis(1, xlim=c(0, 3000), pos=0, lwd=2, labels=FALSE)
axis(2, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(bar(italic(w))), side=2, line=1.9, cex=1.3)
text(3000*0.03, 0.96, "a", cex=1.5)

# average fL
plotBy(data$fL*100 ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000),
       ylim=c(0, 100),
       xlab="MAP (mm)",
       ylab=NA,
       xaxt="n",
       yaxt="n",
       cex.lab=1.3,
       col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 100), pos=0, lwd=2)
mtext("MAP (mm)", side=1, line=2.2, cex=1.3)
mtext(expression(italic(f[L])~"(%)"), side=2, line=1.9, cex=1.3)
legend("topright", c("0.025", "0.05", "0.1"), lty=c(1, 1, 1), col=Cols, title=expression(italic(k)))
text(3000*0.03, 96, "b", cex=1.5)
