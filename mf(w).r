
#options(digits=20)
mf1 <- function(gs,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  ps <- pe*w^(-b)
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
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05){
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  pxmin <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
  gsmax <- (ps-pxmin)*h2*kf(pxmin)/(h*VPD)
  return(gsmax)
}
gsmax2 <- Vectorize(gsmax1)

# Figure
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
w <- 0.15
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax,
      main=expression(italic(m(w,g[s]))),
      xlab=expression(italic(g[s])), 
      ylab=expression(italic(m)~(mu*mol~m^-2~s^-1)),
      xlim=c(0, 1.6), ylim=c(0, 5),
      cex.lab=1.3)
#text(x=1.07, y=1.2, labels=expression(italic(w)==0.6), col="black", cex=1.3)
w <- 0.2
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax, col="red", add=T)
#text(x=0.28, y=3.8, labels=expression(italic(w)==0.4), col="red", cex=1.3)
w <- 1
gsmax <- gsmax2(w)
curve(mf2, 0, gsmax, col="blue", add=T)
#text(x=1.07, y=1.2, labels=expression(italic(w)==0.2), col="blue", cex=1.3)
legend("topleft", c("0.15", "0.2", "1"), lty=c(1, 1, 1), col=c("black", "red", "blue"), title=expression(italic(w)))
