#options(digits=20)
w0f <- function(w0, gs=0,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w0^b*(-(pe/(w0^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w0^(2*b)*(-(pe/(w0^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w0^b*(-(pe/(w0^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)^2))/(w0^b*(-(pe/(w0^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c))
  return(abs(res))
}

w0 <- optimize(w0f, c(0.55, 0.9), tol=.Machine$double.eps)$minimum

ESS1 <- function(w,
                 ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                 a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  pxf1 <- function(gs){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol=.Machine$double.eps)$minimum
    return(px)
  }
  pxf2 <- Vectorize(pxf1)
  
  gsmax1 <- function(w){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    pxmin <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    gsmax <- (ps-pxmin)*h2*kf(pxmin)/(h*VPD)
    return(gsmax)
  }
  gsmax2 <- Vectorize(gsmax1)
  gsmax <- gsmax2(w)
  
  f <- function(gs){
    px <- pxf2(gs)
    test <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*px*(-(px/d))^c*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*(w+w0-w*w0)^(2*b)-c^2*h*h3*(-(px/d))^(2*c)*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*(w+w0-w*w0)^b*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b))-sqrt(c*h*h3*(cp+Km)*(-(px/d))^c*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*(w+w0-w*w0)^b*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b)))*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(2*h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b)))^2))/((-(px/d))^c*(w+w0-w*w0)^b)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*px*(w+w0-w*w0)^b+c*(-(px/d))^c*(h*h3*VPD*(w+w0-w*w0)^b+ca*h2*kmax*LAI*(pe-px*(w+w0-w*w0)^b)+h2*Km*kmax*LAI*(pe-px*(w+w0-w*w0)^b))))
    res <- (gs-test)^2
    return(res)
  }
  
  res <- optimize(f, c(0, gsmax), tol=.Machine$double.eps)$minimum
  return(res)
}
ESS2 <- Vectorize(ESS1)
ESS3 <- Vectorize(function(w)ifelse(w<w0, 0, ESS2((w-w0)/(1-w0))))


gsf <- function(w)ESS2(w)

Af <- function(w,
               ca=400, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1){
  res <- LAI*1/2*(Vcmax+(Km+ca)*gsf(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gsf(w)+((ca+Km)*gsf(w)+Rd)^2-2*Rd*Vcmax)^(1/2))
  return(res)
}

mf1 <- function(w,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
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
w0B1 <- uniroot(Bf1, c(0, 1))$root
w0B2 <- w0+(1-w0)*w0B1

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3, 3.9, 2, 4), mfrow=c(1,1))
curve(Bf2, w0, 1,
      xaxt="n", yaxt="n",
      xlim=c(0, 1), ylim=c(-12, 16), col="red",
      xlab=NA, ylab=NA,
      cex.lab=1.3)

axis(1, xlim=c(0, 1), pos=-12, lwd=2)
mtext(expression(italic(w)),side=1,line=1.7, cex=1.3)
axis(2, ylim=c(-12, 16), pos=0, lwd=2)
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)),side=2,line=1.8, cex=1.3, col="red")
abline(h=0, lty=2, col="red")
abline(v=w0)
text(x=0.51, y=2, labels=expression(italic(w==w[0])), cex=1.3)
abline(v=w0B2)
text(x=0.73, y=2, labels=expression(italic(w==w[B0])), cex=1.3)

par(new=TRUE)
curve(ESS3, 0, 1,
      main="ESS stomatal behaviour",
      xaxt="n", yaxt="n",
      xlim=c(0, 1), ylim=c(0, 1), col="blue",
      xlab=NA, ylab=NA,
      cex.lab=1.3)

axis(4, ylim=c(0, 1), pos=1, lwd=2)
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=4,line=2.4, cex=1.3, col="blue")
