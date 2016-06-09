
source("Normal plants/Functions.r")

# Initializing
ca1 <- 400
ESSf1 <- Vectorize(function(w)ESSf(w, ca1))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca1))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3, 3.9, 2, 4), mfrow=c(1,1))
curve(ESSf1, wL, 1,
      main=expression(ESS~italic(g[s])(italic(w))*" & "*italic(B)(italic(w))),
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 0.4),
      cex.lab=1.3, col="blue")

axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(w)),side=1,line=1.7, cex=1.3)
axis(2, ylim=c(0, 1.2), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

par(new=TRUE)

curve(ESSBf1, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 15),
      cex.lab=1.3, col="red")

axis(4, ylim=c(0, 15), pos=1, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=2.4, cex=1.3)
text(x=0.51, y=14.2, labels=expression(italic(B)), cex=1.3, col="red")
text(x=0.51, y=10.8, labels=expression(italic(g[s])), cex=1.3, col="blue")
