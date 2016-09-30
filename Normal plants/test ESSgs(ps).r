
source("Normal plants/tf gs(ps).r")

# Initializing
ca <- 400
pe <- -1.58*10^-3

ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

# Figures
Cols <- c("black", "blue", "red")
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1), mfrow=c(1, 1))
plot(0, 0,
     type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
     xlim=c(0, 1), ylim=c(0, 1), cex.lab=1.3)

axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(psi[s]),side=1,line=2.3, cex=1.3)
axis(2, ylim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(g[s])),side=2,line=1.8, cex=1.3)

# Sensitivity analysis
SA <- c(10, 50, 100)
h3 <- SA[1]

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- psf(wL)
f1 <- Vectorize(function(x)ESSpsf((pe-psL)*x+psL, ca)/ESSpsf(pe, ca))
curve(f1, 0, 1, col=Cols[1], add=T)
h3 <- SA[2]
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- psf(wL)
f1 <- Vectorize(function(x)ESSpsf((pe-psL)*x+psL, ca)/ESSpsf(pe, ca))
curve(f1, 0, 1, col=Cols[2], add=T)
h3 <- SA[3]
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- psf(wL)
f1 <- Vectorize(function(x)ESSpsf((pe-psL)*x+psL, ca)/ESSpsf(pe, ca))
curve(f1, 0, 1, col=Cols[3], add=T)
