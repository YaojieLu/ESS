
source("Functions.r")

# Initializing
ca <- 400
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
# mf(w, gs)
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
with(data.frame(w=c(1)),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       curve(mf1, 0, gsmaxf(w),
             main=expression(italic(m(w,g[s]))),
             xlab=expression(italic(g[s])), 
             ylab=expression(italic(m)~(mu*mol~m^-2~s^-1)),
             xlim=c(0, 1.6), ylim=c(0, 5), cex.lab=1.3, col="black")
     }
)
with(data.frame(w=c(0.2)),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       curve(mf1, 0, gsmaxf(w), col="red", add=T)
     }
)
with(data.frame(w=c(0.15)),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       curve(mf1, 0, gsmaxf(w), col="blue", add=T)
     }
)
legend("topleft", c("0.15", "0.2", "1"), lty=c(1, 1, 1), col=c("blue", "red", "black"), title=expression(italic(w)))

# ESSAf(w)
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
ESSAf <- Vectorize(function(w)Af(ESSf(w, ca), ca))
curve(ESSAf, wL, 1,
      main=expression(ESS~italic(A)(italic(w))),
      xlab=expression(italic(g[s])), 
      ylab=expression(italic(A)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1), ylim=c(0, 20), cex.lab=1.3, col="black")

# ESSmf(w)
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
ESSmf <- Vectorize(function(w)mf(w, ESSf(w, ca)))
curve(ESSmf, wL, 1,
      main=expression(ESS~italic(m)(italic(w))),
      xlab=expression(italic(g[s])), 
      ylab=expression(italic(m)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1), ylim=c(0, 20), cex.lab=1.3, col="black")
