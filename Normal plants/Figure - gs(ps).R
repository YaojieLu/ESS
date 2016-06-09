
source("Normal plants/Functions.r")

# ESS gs along the soil water potential gradient
ESSpsf <- function(ps, ca, pe=-1.58*10^-3, b=4.38){
  w <- (ps/pe)^(-1/b)
  res <- ESSf(w, ca)
  return(res)
}

# Initializing
ca <- 400
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
ESSpsf1 <- Vectorize(function(ps)ESSpsf(ps, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- with(data.frame(pe=c(-1.58*10^-3), b=c(4.38)), pe*wL^(-b))

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
with(data.frame(pe=c(-1.58*10^-3)),
     curve(ESSpsf1, psL, pe,
           main=expression(ESS~italic(g[s])(psi[soil])),
           xlab=expression(psi[soil]*" (Mpa)"),
           ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
           xlim=c(-4, 0), ylim=c(0, 0.3),
           cex.lab=1.3, col="blue")
)
