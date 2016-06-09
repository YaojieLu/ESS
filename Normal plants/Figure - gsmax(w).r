
source("Normal plants/Functions.r")

# Initializing
pkx <- 0.5
ca <- 400
ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
gsmaxf1 <- Vectorize(gsmaxf)

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(gsmaxf1, 0, 1, xlim=c(0, 1), ylim=c(0, 0.5), cex.lab=1.3,
      main=expression("Maximum"~italic(g[s])),
      xlab=expression(italic(w)),
      ylab=expression(italic(g[smax])~(mol~m^-2~s^-1)))
