
source("Normal plants/Functions.r")

# Parameters
pe <- -1.58*10^-3

# Initializing
ca <- 400
ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
ESSPLCpsf1 <- Vectorize(function(w)ESSPLCpsf(w, ca))
ESSpsf1 <- Vectorize(function(ps)ESSpsf(ps, ca))
#ESSBpsf1 <- Vectorize(function(ps)ESSBpsf(ps, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- with(data.frame(pe=c(-1.58*10^-3), b=c(4.38)), pe*wL^(-b))

# Figures
windows(8, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 3.5), mfrow=c(2, 1))
# ESS gs(w)
par(yaxs="r")
curve(ESSf1, wL, 1,
      #main=expression(ESS~italic(g[s])(italic(w))*" & "*italic(B)(italic(w))),
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 0.4),
      cex.lab=1.3, col="blue")

axis(1, xlim=c(0, 1), pos=-0.4*0.04, lwd=2)
mtext(expression(italic(w)),side=1,line=1.9, cex=1.3)
axis(2, ylim=c(0, 0.4), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

par(new=TRUE)

curve(ESSBf1, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 15),
      cex.lab=1.3, col="red")

segments(wL, 0, wL, ESSf1(wL)/0.4*15, col="blue", lty=2)
segments(0, 0, wL, 0, col="blue")
axis(4, ylim=c(0, 15), pos=1, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=2.6, cex=1.3)
text(x=0.5, y=13, labels=expression(italic(B)), cex=1.3, col="red")
text(x=0.5, y=9.55, labels=expression(italic(g[s])), cex=1.3, col="blue")
text(0.05, 15*1.04*0.95, "a", cex=1.5)

# ESS gs(ps)
curve(ESSpsf1, psL, pe,
      #main=expression(ESS~italic(g[s])(psi[soil])*" & "*italic(B)(psi[soil])),
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(-4, 0), ylim=c(0, 0.4),
      cex.lab=1.3, col="blue")

segments(psL, 0, psL, ESSpsf1(psL), col="blue", lty=2)
segments(-4, 0, psL, 0, col="blue")

axis(1, xlim=c(-4, 0), pos=-0.4*0.04, lwd=2)
mtext(expression(psi[soil]*" (MPa)"),side=1,line=2.2, cex=1.3)
axis(2, ylim=c(0, 0.4), pos=-4, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

par(new=TRUE)

# ESS PLC(w)
curve(ESSPLCpsf1, -4, pe,
      #main=expression(ESS~PLC(italic(w))),
      xlab=NA, ylab=NA, xaxt="n", yaxt="n",
      xlim=c(-4, 0), ylim=c(0, 100),
      cex.lab=1.3)

axis(4, ylim=c(0, 100), pos=0, lwd=2, at=c(0, 25, 50, 75, 100))
mtext("PLC (%)",side=4,line=2, cex=1.3)
text(x=-2, y=22, labels=expression(italic(g[s])), cex=1.3, col="blue")
text(x=-2, y=49, labels="PLC", cex=1.3)
text(4*0.95, 95, "b", cex=1.5)

## ESS gs(ps)
#curve(ESSpsf1, psL, pe,
#      #main=expression(ESS~italic(g[s])(psi[soil])*" & "*italic(B)(psi[soil])),
#      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
#      xlim=c(-4, 0), ylim=c(0, 0.4),
#      cex.lab=1.3, col="blue")
#
#axis(1, xlim=c(-4, 0), pos=-0.4*0.04, lwd=2)
#mtext(expression(psi[soil]*" (MPa)"),side=1,line=2.9, cex=1.3)
#axis(2, ylim=c(0, 0.4), pos=-4, lwd=2)
#mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)
#
#par(new=TRUE)
#
#curve(ESSBpsf1, psL, pe,
#      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
#      xlim=c(-4, 0), ylim=c(0, 15),
#      cex.lab=1.3, col="red")
#
#segments(psL, 0, psL, ESSpsf1(psL)/0.4*15, col="blue", lty=2)
#segments(-4, 0, psL, 0, col="blue")
#axis(4, ylim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
#mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=3.3, cex=1.3)
#text(x=-1.96, y=7.45, labels=expression(italic(B)), cex=1.3, col="red")
#text(x=-1.96, y=3.27, labels=expression(italic(g[s])), cex=1.3, col="blue")
#text(-4*0.95, 15*1.04*0.95, "b", cex=1.5)
#
## ESS PLC(w)
#par(yaxs="i")
#curve(ESSPLCf1, wL, 1,
#      #main=expression(ESS~PLC(italic(w))),
#      xlab=NA, ylab=NA,
#      xlim=c(0, 1), ylim=c(0, 100), yaxt="n",
#      cex.lab=1.3)
#
#axis(2, ylim=c(0, 100), pos=0, lwd=2, at=c(0, 25, 50, 75, 100))
#mtext("PLC (%)",side=2,line=2, cex=1.3)
#mtext(expression(italic(w)),side=1,line=1.9, cex=1.3)
#text(0.95, 95, "c", cex=1.5)
#