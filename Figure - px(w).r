
source("Functions.r")

# ESS px(w)
ESSpxf <- function(w, ca,
                   a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  ps <- pe*w^(-b)
  gs <- ESSf(w, ca)
  f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
  res <- optimize(f1, c(pxminf(w), ps), tol=.Machine$double.eps)
  return(res$minimum)
}

# Initializing
ca <- 400
ESSpxf1 <- Vectorize(function(w)ESSpxf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(ESSpxf1, wL, 1,
      main=expression(ESS~psi[x](italic(w))),
      xlab=expression(italic(w)),
      ylab=expression(psi[x]("Mpa")),
      xlim=c(0, 1), ylim=c(-15, 0),
      cex.lab=1.3, col="blue")
