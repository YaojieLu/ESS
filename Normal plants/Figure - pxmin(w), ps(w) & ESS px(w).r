
source("Normal plants/Functions.r")

# pxmin(w)
pxminf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
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
  pxmin <- pxminf(w)
  return(pxmin)
}

# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# ESS px(w)
ESSpxf <- function(w,
                     a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                     pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  px <- pxf(w, ESSf(w, ca))
  return(px)
}

# Initializing
ca <- 400
pxminf1 <- Vectorize(pxminf)
ESSpxf1 <- Vectorize(ESSpxf)

ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

# Figure
#windows(8, 6)
#par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1, 1))
par(yaxs="i")
curve(pxminf1, wL, 1, xlim=c(0, 1), ylim=c(-5, 0), cex.lab=1.3,
      xlab=expression(italic(w)),
      ylab=expression(italic(psi[x])~(MPa)), col="darkgreen")
curve(psf, wL, 1, add=T, col="red")
curve(ESSpxf1, wL, 1, add=T, col="blue")
abline(v=wL, lwd=1, lty=2)
text(0.11, -4.8, expression(italic(w==w[L])), cex=1.5)
legend("bottomright", expression(italic(psi[xmin]), italic(psi[s]), ESS~italic(psi[x])),
       col = c("darkgreen","red", "blue"), lty=c(1, 1, 1), lwd=c(2, 2, 2))
text(0.025, -0.14, "b", cex=1.5)
