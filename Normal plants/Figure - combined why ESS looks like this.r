
source("Normal plants/Functions.r")

ca <- 400
gsmaxf1 <- Vectorize(gsmaxf)
ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
inversegsmaxf <- function(gs)optimize((function(w)(gsmaxf(w)-gs)^2), c(0, 1), tol=.Machine$double.eps)$minimum
Bgsmaxf <- Vectorize(function(gs)Bf(inversegsmaxf(gs), gs, ca))
inverseESSgsf <- function(gs)optimize((function(w)(ESSf(w, ca)-gs)^2), c(0, 1), tol=.Machine$double.eps)$minimum
BESSgsf <- Vectorize(function(gs)Bf(inverseESSgsf(gs), gs, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
windows(8, 12)
par(mgp=c(2.2, 1, 0), lwd=2, mar=c(3.5, 3.5, 1, 1), mfrow=c(2, 1))
# B(gs) at given w
par(xaxs="r", yaxs="r")
with(data.frame(w=c(1)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w),
             xlab=NA, ylab=NA, xaxt="n",
             xlim=c(0, 0.3), ylim=c(-10, 15), cex.lab=1.3, lwd=1)
       text(gsmaxf(w)+0.02, Bgsf(gsmaxf(w))-0.5, as.expression(bquote(italic(w)==.(w))), cex=1.3)
     }
)
with(data.frame(w=c(0.2)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), add=T, lwd=1)
       text(gsmaxf(w)+0.025, Bgsf(gsmaxf(w))-0.5, as.expression(bquote(italic(w)==.(w))), cex=1.3)
     }
)
with(data.frame(w=wL),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), add=T, lwd=1)
       text(gsmaxf(w)+0.025, Bgsf(gsmaxf(w))-0.8, expression(italic(w==w[L])), cex=1.3)
     }
)

# B(gs) at given w with gs is indicated by gs=0, gsmaxf and ESS
segments(0, Bf(0, 0, ca), 0, Bf(1, 0, ca), col="red")
curve(Bgsmaxf, gsmaxf(0.14), gsmaxf(1), col="orange", add=T)
curve(BESSgsf, ESSf(wL, ca), ESSf(1, ca), col="blue", add=T)
segments(-1, 0, ESSf(wL, ca), Bf(wL, ESSf(wL, ca), ca), lwd=1, lty=3)

axis(1, xlim=c(0, 0.3), pos=-10-25*0.04, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=1, line=2.6, cex=1.3)
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)), side=2, line=1.8, cex=1.3)
legend("bottomright", expression(italic(g[smax]), italic(g[s]==0), ESS~italic(g[s])),
       col = c("orange","red", "blue"), lty=c(1, 1, 1), lwd=c(2, 2, 2))
text(0.303, 15, "a", cex=1.5)
box()

# ESS gs and gsmax
par(xaxs="i", yaxs="r")
curve(gsmaxf1, 0, 1, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.3,
      xlab=NA, ylab=NA, yaxt="n", col="orange")

abline(h=0, col="red")

curve(ESSf1, wL, 1, add=T, col="blue")
segments(0, 0, wL, 0, col="blue", lty=2)

abline(v=wL, lwd=1, lty=2)

axis(2, ylim=c(0, 0.3), pos=0, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(w)), side=1, line=1.9, cex=1.3)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=2, line=1.8, cex=1.3)
text(0.11, 0.25, expression(italic(w==w[L])), cex=1.3)
text(0.9722, 0.3, "b", cex=1.5)
box()

dev.copy2pdf(file = "Normal plants//Figures/Figure 3.pdf")
