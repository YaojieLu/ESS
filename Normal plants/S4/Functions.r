
############################### Scenario 4 ###############################
# psf(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# xylem conductance function
kxf <- function(px, kxmax=5, c=2.64, d=3.54)kxmax*exp(-(-px/d)^c)

# gsmaxf(w)
gsmaxf <- Vectorize(function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, h2=l*LAI/nZ*p/1000, VPD=0.02){
  ps <- psf(w)
  f1 <- function(px)(ps-px)*h2*kxf(px)/(h*VPD)
  res <- optimize(f1, c(-20, 0), tol=.Machine$double.eps, maximum=T)$objective
  return(res)
})

# Af(gs, ca)
Af <- function(gs, ca, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# integralfnoc of PDF
integralfnocf <- function(ca, k, MAP,
                          LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                          gamma=1/((MAP/365/k)/1000)*nZ){
  
  Ef <- function(w)h*VPD*gsmaxf(w)
  rEf <- function(w)1/Ef(w)
  integralrEf <- Vectorize(function(w)integrate(rEf, wmin, w, rel.tol=.Machine$double.eps^0.25)$value)
  fnoc <- function(w)1/Ef(w)*exp(-gamma*w+k*integralrEf(w))
  
  res <- integrate(fnoc, wmin, 1, rel.tol=.Machine$double.eps^0.25)
  return(res)
}

# PDF of w
PDFf <- function(w, cPDF, ca, k, MAP,
                 LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 gamma=1/((MAP/365/k)/1000)*nZ){
  
  Ef <- function(w)h*VPD*gsmaxf(w)
  rEf <- function(w)1/Ef(w)
  integralrEf <- Vectorize(function(w)integrate(rEf, wmin, w, rel.tol=.Machine$double.eps^0.25)$value)
  res <- cPDF/Ef(w)*exp(-gamma*w+k*integralrEf(w))
  return(res)
}

# Average A
averAf <- function(cPDF, ca, k, MAP){
  Awf <- function(w)Af(gsmaxf(w), ca)
  f1 <- function(w)Awf(w)*PDFf(w, cPDF, ca, k, MAP)
  res <- integrate(f1, wmin, 1, rel.tol=.Machine$double.eps^0.25)
  return(res)
}

# Average E
averEf <- function(cPDF, ca, k, MAP,
                   LAI=1, a=1.6, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02){
  Ef <- function(w)h*VPD*gsmaxf(w)
  f1 <- function(w)Ef(w)*PDFf(w, cPDF, ca, k, MAP)
  res <- integrate(f1, wmin, 1, rel.tol=.Machine$double.eps^0.25)
  return(res)
}

# Average w
averwf <- function(cPDF, ca, k, MAP){
  f1 <- function(w)w*PDFf(w, cPDF, ca, k, MAP)
  res <- integrate(f1, wmin, 1, rel.tol=.Machine$double.eps^0.25)
  return(res)
}
