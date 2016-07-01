
source("Normal plants/Functions.r")

# Initializing
ca <- 400

Af1 <- function(gs)Af(gs, ca)

# Figure
Cols <- c("blue", "red", "darkgreen")
w <- c(0.17, 0.2, 1)
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "r", lwd = 2, mar=c(3.5, 3.5, 0.5, 1), mfrow=c(1, 1))
# A(gs)
curve(Af1, 0, 0.5,
      xlab=NA, ylab=NA, xaxt="n", yaxt="n",
      xlim=c(0, 0.3), ylim=c(0, 15.1), cex.lab=1.3)

mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=1,line=2.6, cex=1.3)
axis(1, xlim=c(0, 0.3), pos=-15*0.04, lwd=2, at=c(0, 0.1, 0.2, 0.3))
axis(2, ylim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(italic(A[N])*" or "*italic(m)~(mu*mol~m^-2~s^-1)),side=2,line=1.7, cex=1.3)

# m(w, gs)
with(list(w=w[1], Col=Cols[1]),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       mf21 <- Vectorize(function(gs)mf2(w, gs))
       curve(mf1, 0, gsmaxf(w), col=Col, add=T)
       curve(mf21, 0, gsmaxf(w), col=Col, lty=2, add=T)
     }
)
with(list(w=w[2], Col=Cols[2]),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       mf21 <- Vectorize(function(gs)mf2(w, gs))
       curve(mf1, 0, gsmaxf(w), col=Col, add=T)
       curve(mf21, 0, gsmaxf(w), col=Col, lty=2, add=T)
     }
)
with(list(w=w[3], Col=Cols[3]),
     {
       mf1 <- Vectorize(function(gs)mf(w, gs))
       mf21 <- Vectorize(function(gs)mf2(w, gs))
       curve(mf1, 0, gsmaxf(w), col=Col, add=T)
       curve(mf21, 0, gsmaxf(w), col=Col, lty=2, add=T)
     }
)

curve(mgsmaxf, gsmaxf(0.1), gsmaxf(1), lwd=1, add=T)

text(0.17, 3.9, expression(italic(g[smax])), cex=1.3)
legend("topleft", expression(italic(A[N]), italic(m)*", at "*italic(w)*" = 0.17", italic(m)*", at "*italic(w)*" = 0.2", italic(m)*", at "*italic(w)*" = 1"), col=c("black", Cols), lty=c(1, 1, 1, 1))
box()

dev.copy2pdf(file = "Normal plants//Figures/Figure 2.pdf")
