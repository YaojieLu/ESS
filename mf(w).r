
#options(digits=20)
w0f <- function(w0, gs=0,
                ca=400, Vcmax=50, cp=30, Km=703, Rd=1,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  res <- ((-c)*h*h2*h3*(ca+Km)*kmax*LAI*pe*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD*w0^b*(-(pe/(w0^b*d)))^c-c^2*h^2*h3^2*(ca*(Rd-Vcmax)+2*cp*Vcmax+Km*(Rd+Vcmax))*VPD^2*w0^(2*b)*(-(pe/(w0^b*d)))^(2*c)-sqrt(c*h*h3*(cp+Km)*Vcmax*(Km*Rd+ca*(Rd-Vcmax)+cp*Vcmax)*VPD*w0^b*(-(pe/(w0^b*d)))^c*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)*(h2*(ca+Km)*kmax*LAI*pe+2*c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c)^2))/(w0^b*(-(pe/(w0^b*d)))^c)/(c*h*h3*(ca+Km)^2*VPD*(h2*(ca+Km)*kmax*LAI*pe+c*h*h3*VPD*w0^b*(-(pe/(w0^b*d)))^c))
  return(abs(res))
}

w0 <- optimize(w0f, c(0.55, 0.9), tol=.Machine$double.eps)$minimum

mf1 <- function(gs,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  ps <- pe*(w0+(1-w0)*w)^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  minimum <- optimize(f1, c(-20,0))$minimum
  f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
  px <- optimize(f2, c(minimum, ps))$minimum
  
  k <- kf(px)
  PLC <- 1-k/kmax
  res <- h3*PLC
  return(res)
}
mf2 <- Vectorize(mf1)

gsmax1 <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05){
  ps <- pe*(w0+(1-w0)*w)^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  pxmin <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
  f2 <- function(gs)((ps-pxmin)*h2*kf(pxmin)-h*VPD*gs)^2
  gsmax <- optimize(f2, c(0, 2), tol = .Machine$double.eps)$minimum
  return(gsmax)
}
gsmax2 <- Vectorize(gsmax1)

# Figure
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
w <- (0.9-w0)/(1-w0)
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax,
      main=expression(italic(m(w,g[s](w)))),
      xlab=expression(italic(g[s])), 
      ylab=expression(italic(m)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1.2), ylim=c(0, 6),
      cex.lab=1.3)
text(x=1.07, y=1.2, labels=expression(italic(w)==0.9), col="black", cex=1.3)
w <- (0.8-w0)/(1-w0)
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax, col="red", add=T)
text(x=0.74, y=1.5, labels=expression(italic(w)==0.8), col="red", cex=1.3)
w <- (0.7-w0)/(1-w0)
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax, col="blue", add=T)
text(x=0.28, y=3.8, labels=expression(italic(w)==0.7), col="blue", cex=1.3)
