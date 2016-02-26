
gsf <- function(w)ESS2(w)

Af <- function(w,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1){
  res <- LAI*1/2*(Vcmax+(Km+ca)*gsf(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gsf(w)+((ca+Km)*gsf(w)+Rd)^2-2*Rd*Vcmax)^(1/2))
  return(res)
}
  
mf1 <- function(w,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){

  ps <- pe*(w0+(1-w0)*w)^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  minimum <- optimize(f1, c(-20,0))$minimum
  f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gsf(w))^2
  px <- optimize(f2, c(minimum, ps))$minimum
  
  k <- kf(px)
  PLC <- 1-k/kmax
  res <- h3*PLC
  return(res)
}
mf2 <- Vectorize(mf1)
Bf1 <- function(w)Af(w)-mf2(w)
Bf2 <- Vectorize(function(w)Bf1((w-w0)/(1-w0)))

pxf1 <- function(w,
                 a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  minimum <- optimize(f1, c(-20,0))$minimum
  f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*ESS3(w))^2
  px <- optimize(f2, c(minimum, ps))$minimum
  return(px)
}

pxf2 <- Vectorize(pxf1)

# Figure
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(Bf2, w0, 1,
      main= expression(italic(B(w,g[sESS](w)))),
      xlab=expression(italic(w)), 
      ylab=expression(italic(B)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1), ylim=c(-12, 16),
      cex.lab=1.3)
abline(h=0, col="red")
abline(v=w0, col="blue")
text(x=0.51, y=2, labels=expression(italic(w==w[0])), col="blue", cex=1.3)
#curve(Af, 0, 1)
#curve(mf2, 0, 1, col="red", add=T)
#curve(Bf, 0, 1, col="blue", add=T)
curve(pxf2, w0, 1,
      xlim=c(0, 1), ylim=c(-20, -2),
      cex.lab=1.3)
