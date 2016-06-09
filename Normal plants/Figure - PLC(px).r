
# PLC(w)
PLCf <- function(x, c=2.64, d=3.54)(1-exp(-(-x/d)^c))*100

f <- function(x)PLCf(x)-50
P50 <- uniroot(f, c(-8, 0), tol=.Machine$double.eps)$root

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 3, 1.5), mfrow=c(1,1))
curve(PLCf, -8, 0, xlim=c(-8, 0), ylim=c(0, 100),
      xlab=expression(psi[x]~(MPa)), ylab="PLC (%)",
      main=as.expression(bquote(P[50]==.(P50)~MPa)), cex.lab=1.3)

