
AverBf <- function(k, MAP, ca=400,
                   nz=0.5, gamma=1/((MAP/365/k)/1000)*nZ,
                   Vcmax=50, cp=30, Km=703, Rd=1, LAI=1,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05, h3=10){
  
  Af <- function(gs)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))
  
  mf <- function(w, gs){
    ps <- pe*(w0+(1-w0)*w)^(-b)
    kf <- function(x)kmax*exp(-(-x/d)^c)
    f1 <- function(x)-((ps-x)*h2*kf(x))
    minimum <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
    f2 <- function(x)((ps-x)*h2*kf(x)-h*VPD*gs)^2
    px <- optimize(f2, c(minimum, ps), tol = .Machine$double.eps)$minimum
    
    k <- kf(px)
    PLC <- 1-k/kmax
    res <- h3*PLC
    return(res)
  }
  
  Bf1 <- function(w, gs)Af(gs)-mf(w, gs)
  Bf2 <- Vectorize(Bf1)
  Bf3 <- function(w)Bf2(w, ESS2(w))
  
  Ev <- function(w){h*VPD*ESS2(w)}
  rEv <- function(w){1/Ev(w)}
  integralEv <- Vectorize(function(w){integrate(rEv, 0.12, w)$value})
  fnoc <- function(w){1/Ev(w)*exp(-gamma*w+k*integralEv(w))}
  integralfnoc <- integrate(fnoc, 0.12, 1)$value
  cESS <- 1/integralfnoc
  f <- function(w){cESS*fnoc(w)}
  Bff <- function(w){f(w)*Bf3(w)}
  
  averB <- integrate(Bff, 0.12, 1)$value
  return(averB)
}

AverBf(k=0.025, MAP=365)
