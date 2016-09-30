
source("Normal plants/Functions.r")

# Initializing
pe <- -1.58*10^-3

PLCf1 <- function(px)PLCf(px)*100

# Figures
Cols <- c("blue", "red", "darkgreen")
w <- c(0.17, 0.2, 1)
windows(8, 6)
par(mgp=c(2.2, 1, 0), lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(2, 2))
# Soil water retention curve
curve(psf, 0, 1, xlim=c(0, 1), ylim=c(-10, 0),
      xlab=expression(italic(w)), ylab=expression(psi[s]~(MPa)), cex.lab=1.3)

points(w[1], psf(w[1]), col=Cols[1], pch=4, cex=1.3)
points(w[2], psf(w[2]), col=Cols[2], pch=4, cex=1.3)
points(w[3], psf(w[3]), col=Cols[3], pch=4, cex=1.3)

text(-1*0.04+1*1.08*0.05, -10.4+10.8*0.95, "a", cex=1.5)
box()

# Vulnerability curve
curve(PLCf1, -20, pe, xlim=c(-10, 0), ylim=c(0, 100),
      xlab=expression(psi[x]~(MPa)), ylab="PLC (%)", cex.lab=1.3)

text(-10*1.04+10*1.08*0.95, -4+108*0.95, "b", cex=1.5)
box()

# Plant water use envelope
par(xaxs ="i", yaxs="i")
with(list(w=w[1], Col=Cols[1]),
     {
       gswpxf1 <- Vectorize(function(px)gswpxf(w, px))
       curve(gswpxf1, pxminf(w), pe,
             xlab=expression(psi[x]~(MPa)), ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
             xlim=c(-10, 0), ylim=c(0, 0.3), yaxt="n", cex.lab=1.3, col=Col)
       curve(gswpxf1, -10, pxminf(w), cex.lab=1.3, col=Col, lty=2, add=T)
     }
)

with(list(w=w[2], Col=Cols[2]),
     {
       gswpxf1 <- Vectorize(function(px)gswpxf(w, px))
       curve(gswpxf1, pxminf(w), pe, cex.lab=1.3, col=Col, add=T)
       curve(gswpxf1, -10, pxminf(w), cex.lab=1.3, col=Col, lty=2, add=T)
     }
)

with(list(w=w[3], Col=Cols[3]),
     {
       gswpxf1 <- Vectorize(function(px)gswpxf(w, px))
       curve(gswpxf1, pxminf(w), pe, cex.lab=1.3, col=Col, add=T)
       curve(gswpxf1, -10, pxminf(w), cex.lab=1.3, col=Col, lty=2, add=T)
     }
)

curve(gspxminf, pxminf(0.01), pxminf(1), lwd=1, add=T)

axis(2, ylim=c(0, 0.3), pos=-10, lwd=2, at=c(0, 0.1, 0.2, 0.3), cex.lab=1)
text(-2.2, 0.13, expression(italic(g[smax])), cex=1.3)
text(-0.5, 0.3*0.95, "c", cex=1.5)
legend("topleft", title=expression(italic(w)), as.character(w), col=Cols, lty=c(1, 1, 1))
box()

# PLC(w, gs)
with(list(w=w[1], Col=Cols[1]),
     {
       PLCwgsf1 <- Vectorize(function(gs)PLCwgsf(w, gs)*100)
       PLCwgsf21 <- Vectorize(function(gs)PLCwgsf2(w, gs)*100)
       
       curve(PLCwgsf1, 0, gsmaxf(w),
             xlab=NA, 
             ylab="PLC (%)",
             xlim=c(0, 0.3), ylim=c(0, 100), xaxt="n", cex.lab=1.3, col=Col)
       curve(PLCwgsf21, 0, gsmaxf(w), lty=2, col=Col, add=T)
       
       mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=1,line=2.6, cex=1)
     }
)
with(list(w=w[2], Col=Cols[2]),
     {
       PLCwgsf1 <- Vectorize(function(gs)PLCwgsf(w, gs)*100)
       PLCwgsf21 <- Vectorize(function(gs)PLCwgsf2(w, gs)*100)
       
       curve(PLCwgsf1, 0, gsmaxf(w), col=Col, add=T)
       curve(PLCwgsf21, 0, gsmaxf(w), lty=2, col=Col, add=T)
     }
)
with(list(w=w[3], Col=Cols[3]),
     {
       PLCwgsf1 <- Vectorize(function(gs)PLCwgsf(w, gs)*100)
       PLCwgsf21 <- Vectorize(function(gs)PLCwgsf2(w, gs)*100)
       
       curve(PLCwgsf1, 0, gsmaxf(w), col=Col, add=T)
       curve(PLCwgsf21, 0, gsmaxf(w), lty=2, col=Col, add=T)
     }
)

curve(PLCgsmaxf, gsmaxf(0.1), gsmaxf(1), lwd=1, add=T)

axis(1, xlim=c(0, 0.3), pos=0, lwd=2, at=c(0, 0.1, 0.2, 0.3), cex.lab=1)
text(0.16, 39, expression(italic(g[smax])), cex=1.3)
text(0.3*0.95, 95, "d", cex=1.5)
box()

dev.copy2pdf(file = "Normal plants//Figures/Figure 1.pdf")
