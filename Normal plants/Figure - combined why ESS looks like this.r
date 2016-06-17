
source("Normal plants/Functions.r")

ca <- 400
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
inversegsmaxf <- function(gs)optimize((function(w)(gsmaxf(w)-gs)^2), c(0, 1), tol=.Machine$double.eps)$minimum
Bgsmaxf <- Vectorize(function(gs)Bf(inversegsmaxf(gs), gs, ca))
inverseESSgsf <- function(gs)optimize((function(w)(ESSf(w, ca)-gs)^2), c(0, 1), tol=.Machine$double.eps)$minimum
BESSgsf <- Vectorize(function(gs)Bf(inverseESSgsf(gs), gs, ca))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
#windows(8, 6)
#par(mgp=c(2.2, 1, 0), lwd=2, yaxs="i", mar=c(4, 4, 1, 1), mfrow=c(1,1))#, xaxs="i"
par(xaxs="r", yaxs="i")
# B(gs) at given w
with(data.frame(w=c(1)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w),
             xlab=NA, ylab=expression(italic(B)~(mu*mol~m^-2~s^-1)),
             xlim=c(0, 0.5), ylim=c(-10, 15), xaxt="n", cex.lab=1.3, lwd=1)
       segments(-1, ESSBf1(w), gsmaxf(w), ESSBf1(w), lty=2, lwd=1)
     }
)
with(data.frame(w=c(0.2)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), add=T, lwd=1)
       segments(-1, ESSBf1(w), gsmaxf(w), ESSBf1(w), lty=2, lwd=1)
     }
)
with(data.frame(w=wL),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), add=T, lwd=1)
       segments(-1, ESSBf1(w), gsmaxf(w), ESSBf1(w), lty=2, lwd=1)
     }
)

# B(gs) at given w with gs is indicated by gs=0, gsmaxf and ESS
segments(0, Bf(wL, 0, ca), 0, Bf(1, 0, ca), col="red")
curve(Bgsmaxf, gsmaxf(wL), gsmaxf(1), col="darkgreen", add=T)
curve(BESSgsf, ESSf(wL, ca), ESSf(1, ca), col="blue", add=T)

# legend and labels
axis(1, xlim=c(0, 0.5), pos=-10, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=1, line=3)
text(0.435, 11, expression(italic(w==1)), cex=1.3)
text(0.183, 6.9, expression(italic(w==0.2)), cex=1.3)
text(0.06, -0.8, expression(italic(w==w[L])), cex=1.3)
legend("bottomright", title=expression(italic(g[s](w)*"=")), expression(italic(g[smax]), "0", ESS~italic(g[s])),
       col = c("darkgreen","red", "blue"), lty=c(1, 1, 1), lwd=c(2, 2, 2))
text(0.507, 14.3, "c", cex=1.5)
