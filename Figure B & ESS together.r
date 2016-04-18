
options(digits=20)

# wL
f <- function(w,
              ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
              a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
              pe=-1.58*10^-3, b=4.38, h2=h/1000, xkmax=5, c=5.71, d=10.05, h3=10){
  
  # photosynthesis rate function
  Af <- function(gs){LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-sqrt(((ca-Km)*gs+Rd-Vcmax)^2+4*gs*((ca*gs+Rd)*Km+cp*Vcmax)))}
  # xylem conductance function
  xkf <- function(x)xkmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-(ps-x)*h2*xkf(x)
    res <- optimize(f1, c(-200,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*xkf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  # PLC cost function
  mf <- function(w, gs){
    px <- pxf(w, gs)
    xk <- xkf(px)
    PLC <- 1-xk/xkmax
    res <- h3*PLC
    return(res)
  }
  # maximum gs function for given w
  gsmaxf <- function(w){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    gsmax <- (ps-pxmin)*h2*xkf(pxmin)/(h*VPD)
    return(gsmax)
  }
  # net photosynthesis rate function
  Bf <- function(w, gs)Af(gs)-mf(w, gs)
  
  f1 <- function(gs)Bf(w, gs)
  gsmax <- gsmaxf(w)
  res <- optimize(f1, c(0, gsmax), tol=.Machine$double.eps, maximum=T)
  #browser()
  return(res)
}
ESSBf <- Vectorize(function(w)f(w)$objective)
ESSf <- Vectorize(function(w)f(w)$maximum)
wL <- uniroot(ESSBf, c(0.12, 1), tol=.Machine$double.eps)$root

# w0
w0f <- function(w0, gs=0,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w0^b*(-(pe/(w0^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w0^(2*b)*(-(pe/(w0^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w0^b*(-(pe/(w0^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)^2))/(w0^b*(-(pe/(w0^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c))
  return(abs(res))
}

w0 <- optimize(w0f, c(0.113, 1), tol=.Machine$double.eps)$minimum

# ESS A, m, and B
Af <- function(w, ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*ESSf(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*ESSf(w)+((ca+Km)*ESSf(w)+Rd)^2-2*Rd*Vcmax)^(1/2))

mf1 <- function(w,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  minimum <- optimize(f1, c(-20,0))$minimum
  f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*ESSf(w))^2
  px <- optimize(f2, c(minimum, ps))$minimum
  
  k <- kf(px)
  PLC <- 1-k/kmax
  res <- h3*PLC
  return(res)
}
mf2 <- Vectorize(mf1)

Bf <- function(w)Af(w)-mf2(w)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3, 3.9, 2, 4), mfrow=c(1,1))
curve(Bf, wL, 1,
      xaxt="n", yaxt="n",
      xlim=c(0, 1), ylim=c(0, 16),
      xlab=NA, ylab=NA, col="red",
      cex.lab=1.3)

axis(1, xlim=c(wL, 1), pos=0, lwd=2)
mtext(expression(italic(w)),side=1,line=1.7, cex=1.3)
axis(2, ylim=c(0, 12), pos=0, lwd=2)
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3, col="red")
#abline(h=0, lty=2, col="red")
#abline(v=w0)
#text(x=0.51, y=2, labels=expression(italic(w==w[0])), cex=1.3)
#abline(v=wL)
#text(x=0.73, y=2, labels=expression(italic(w==w[B0])), cex=1.3)

par(new=TRUE)
curve(ESSf, wL, 1,
      main="ESS stomatal behaviour",
      xaxt="n", yaxt="n",
      xlim=c(0, 1), ylim=c(0, 1.2),
      xlab=NA, ylab=NA, col="blue",
      cex.lab=1.3)

axis(4, ylim=c(0, 1.2), pos=1, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=4,line=2.4, cex=1.3, col="blue")
abline(h=ESSf(wL))
