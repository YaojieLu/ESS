
source("Normal plants/Functions.r")

ca <- 400
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,1))
with(data.frame(w=c(1)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w),
             xlab=expression(italic(g[s])~(mol~m^-2~s^-1)),
             ylab=expression(italic(B)~(mu*mol~m^-2~s^-1)),
             xlim=c(0, 0.5), ylim=c(-10, 20), cex.lab=1.3, col="black")
       ESSgs <- ESSf(w, ca)
       points(ESSgs, Bgsf(ESSgs), type='p', pch="|", col="black")
     }
)
with(data.frame(w=c(0.2)),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), col="red", add=T)
       ESSgs <- ESSf(w, ca)
       points(ESSgs, Bgsf(ESSgs), type='p', pch="|", col="red")
     }
)
with(data.frame(w=wL),
     {
       Bgsf <- Vectorize(function(gs)Bf(w, gs, ca))
       curve(Bgsf, 0, gsmaxf(w), col="blue", add=T)
       ESSgs <- ESSf(w, ca)
       points(ESSgs, Bgsf(ESSgs), type='p', pch="|", col="blue")
     }
)
abline(h=0, lty=2)
legend("topright", c("1", "0.2", expression(w[L])), lty=c(1, 1, 1), col=c("black", "red", "blue"), title=expression(italic(w)))
