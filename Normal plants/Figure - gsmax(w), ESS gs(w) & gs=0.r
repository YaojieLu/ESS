
source("Normal plants/Functions.r")

# Initializing
ca <- 400
gsmaxf1 <- Vectorize(gsmaxf)
ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
windows(8, 18)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 0.5, 1), mfrow=c(3, 1))
par(yaxs="r")
curve(gsmaxf1, 0, 1, xlim=c(0, 1), ylim=c(0, 0.5), cex.lab=1.3,
      xlab=expression(italic(w)),
      ylab=expression(italic(g[s])~(mol~m^-2~s^-1)), col="darkgreen")

curve(ESSf1, wL, 1, add=T, col="blue")
segments(wL, 0, wL, ESSf1(wL), col="blue", lty=2)

abline(h=0, col="red")

legend("topleft", title=expression(italic(g[s](w)*"=")), expression(italic(g[smax]), "0", ESS~italic(g[s])),
       col = c("darkgreen","red", "blue"), lty=c(1, 1, 1), lwd=c(2, 2, 2))
text(0.977, 0.505, "a", cex=1.5)
