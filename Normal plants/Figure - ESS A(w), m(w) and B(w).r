
source("Normal plants/Functions.r")

# Initializing
ca <- 400
ESSAf1 <- Vectorize(function(w)ESSAf(w, ca))
ESSmf1 <- Vectorize(function(w)ESSmf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,1))
curve(ESSAf1, wL, 1,
      xlab=expression(italic(w)), yaxt="n", ylab=NA,
      xlim=c(0, 1), ylim=c(0, 15),
      cex.lab=1.3, col="red")

curve(ESSmf1, wL, 1, col="blue", add=T)
curve(ESSBf1, wL, 1, add=T)

axis(2, ylim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(A)*", "*italic(m~or~B)~(mu*mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)
text(0.2, 13.5, expression(A), col="red", cex=1.5)
text(0.28, 1.5, expression(m), col="blue", cex=1.5)
text(0.27, 12.1, expression(B), col="black", cex=1.5)
