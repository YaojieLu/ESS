
# B(w, gs)
Bf <- function(w, gs, ca,
               Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
               a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
               pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10){
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
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
    return(res)
  }
  # m(w, gs)
  mf <- function(w, gs){
    px <- pxf(w, gs)
    kx <- kxf(px)
    PLC <- 1-kx/kxmax
    res <- h3*PLC
    return(res)
  }
  # A(gs)
  Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))
  
  res <- Af(gs)-mf(w, gs)
  return(res)
}
# ESS stomatal behaviour function
ESSf1 <- function(w, ca,
                  Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05, h3=10){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
    return(gsmax)
  }
  gsmax <- gsmaxf(w)
  
  f1 <- function(gs)Bf(w, gs, ca)
  res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}
# ESS B(w)
ESSBf1 <- function(w, ca)Bf(w, ESSf1(w, ca), ca)

# Initializing
ca1 <- 400
ESSf2 <- Vectorize(function(w)ESSf1(w, ca1))
ESSBf2 <- Vectorize(function(w)ESSBf1(w, ca1))

wL <- uniroot(ESSBf2, c(0.12, 1), tol=.Machine$double.eps)$root

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3, 3.9, 2, 4), mfrow=c(1,1))
curve(ESSf2, wL, 1,
      main=expression(ESS~italic(g[s])(italic(w))*" & "*italic(B)(italic(w))),
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 1.2),
      cex.lab=1.3, col="blue")

axis(1, xlim=c(0, 1), pos=0, lwd=2)
mtext(expression(italic(w)),side=1,line=1.7, cex=1.3)
axis(2, ylim=c(0, 1.2), pos=0, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3)

par(new=TRUE)

curve(ESSBf2, wL, 1,
      xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(0, 1), ylim=c(0, 20),
      cex.lab=1.3, col="red")

axis(4, ylim=c(0, 20), pos=1, lwd=2)
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=4,line=2.4, cex=1.3)
text(x=0.51, y=14.7, labels=expression(italic(B)), cex=1.3, col="red")
text(x=0.51, y=18.2, labels=expression(italic(g[s])), cex=1.3, col="blue")
