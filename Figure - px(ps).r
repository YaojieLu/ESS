
source("Functions.r")

# ESS px(ps)
ESSpxpsf <- function(ps,
                     a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                     pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
  w <- (ps/pe)^(-1/b)
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
ESSpxpsf1 <- Vectorize(ESSpxpsf)
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))

wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
psL <- with(data.frame(pe=c(-1.58*10^-3), b=c(4.38)), pe*wL^(-b))

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
with(data.frame(pe=c(-1.58*10^-3)),
     curve(ESSpxpsf1, psL, pe,
           main=expression(ESS~psi[xylem]*(psi[soil])),
           xlab=expression(psi[soil]*" (Mpa)"),
           ylab=expression(psi[xylem]*" (Mpa)"),
           xlim=c(-11, 0), ylim=c(-12, 0),
           cex.lab=1.3, col="blue")
)
