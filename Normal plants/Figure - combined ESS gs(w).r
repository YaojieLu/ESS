
source("Normal plants/Functions.r")

# Initializing
ca <- 400
pe <- -1.58*10^-3

ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
ESSPLCf1 <- Vectorize(function(w)ESSPLCf(w, ca)*100)
PLCwf <- function(w)PLCf(psf(w))*100
ESSpsf1 <- Vectorize(function(ps)ESSpsf(ps, ca))
ESSBpsf1 <- Vectorize(function(ps)ESSBpsf(ps, ca))
ESSPLCpsf1 <- Vectorize(function(ps)ESSPLCpsf(ps, ca)*100)
PLCf1 <- function(px)PLCf(px)*100

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- psf(wL)

# Figures
windows(16, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 8), mfrow=c(2, 1))
# On w
# B
curve(ESSBf1, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 15),
      cex.lab=1.3, col="red")

segments(0, 0, wL, 0, col="red")

axis(4, ylim=c(0, 15), pos=1, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=2.6, cex=1.3)

# gs
par(new=TRUE)
curve(ESSf1, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 0.3),
      cex.lab=1.3, col="blue")

segments(0, 0, wL, 0, col="blue", lty=2)

axis(1, xlim=c(0, 1), pos=-0.3*0.04, lwd=2)
mtext(expression(italic(w)),side=1,line=1.9, cex=1.3)
axis(2, ylim=c(0, 0.3), pos=0, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

# PLC
par(new=TRUE)
curve(ESSPLCf1, wL, 1,
      xlab=NA, ylab=NA, xaxt="n", yaxt="n",
      xlim=c(0, 1), ylim=c(0, 100),
      cex.lab=1.3)

curve(PLCwf, 0, wL, add=T)

abline(v=wL, lwd=1, lty=2)

axis(4, ylim=c(0, 100), pos=1.082, lwd=2, at=c(0, 25, 50, 75, 100))
mtext("PLC (%)",side=4,line=6.3, cex=1.3)

# text
text(x=0.5, y=66, labels=expression(italic(g[s])), cex=1.3, col="blue")
text(x=0.5, y=89, labels=expression(italic(B)), cex=1.3, col="red")
text(x=0.5, y=12, labels="PLC", cex=1.3, col="black")
text(0.215, 93, expression(italic(w==w[L])), cex=1.3)
text(0.98, 100*1.04*0.95, "a", cex=1.5)
box()

# On ps
# B
curve(ESSBpsf1, psL, pe,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(-4, 0), ylim=c(0, 15),
      cex.lab=1.3, col="red")

segments(-4, 0, psL, 0, col="red")

axis(4, ylim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=2.6, cex=1.3)

# gs
par(new=TRUE)
curve(ESSpsf1, psL, pe,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(-4, 0), ylim=c(0, 0.3),
      cex.lab=1.3, col="blue")

segments(-4, 0, psL, 0, col="blue", lty=2)

axis(1, xlim=c(-4, 0), pos=-0.3*0.04, lwd=2)
mtext(expression(psi[s]*" (MPa)"),side=1,line=2.3, cex=1.3)
axis(2, ylim=c(0, 0.3), pos=-4, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

# PLC
par(new=TRUE)
curve(ESSPLCpsf1, psL, pe,
      xlab=NA, ylab=NA, xaxt="n", yaxt="n",
      xlim=c(-4, 0), ylim=c(0, 100),
      cex.lab=1.3)

curve(PLCf1, -4, psL, add=T)

abline(v=psf(wL), lwd=1, lty=2)

axis(4, ylim=c(0, 100), pos=0.33, lwd=2, at=c(0, 25, 50, 75, 100))
mtext("PLC (%)",side=4,line=6.3, cex=1.3)

# text
text(x=-2, y=19, labels=expression(italic(g[s])), cex=1.3, col="blue")
text(x=-2, y=49, labels=expression(italic(B)), cex=1.3, col="red")
text(x=-2, y=35, labels="PLC", cex=1.3, col="black")
text(-3, 93, expression(italic(w==w[L])), cex=1.3)
text(-4*0.02, 100*1.04*0.95, "b", cex=1.5)
box()

dev.copy2pdf(file = "Normal plants//Figures/Figure 4.pdf")
