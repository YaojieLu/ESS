
source("Functions.r")

# Initializing
ca <- 400

Af1 <- function(gs)Af(gs, ca)

windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 1, 1), mfrow=c(1,1))

curve(Af1, 0, 1,
      xaxt="n", yaxt="n",
      xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 20), cex.lab=1.3, col="black")
axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=1,line=2.7, cex=1.3)
axis(2, ylim=c(0, 20), pos=0, lwd=2)
mtext(expression(italic(A)~(mu*mol~m^-2~s^-1)),side=2,line=2, cex=1.3)
