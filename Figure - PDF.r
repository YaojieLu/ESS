
options(digits=20)
source("Functions.r")

# Initializing
ca1 <- 400
k1 <- 0.1
MAP1 <- 365*2

# Results
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca1))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root

integralfnoc <- integralfnocf(ca1, k1, MAP1, wL)
cPDF <- 1/(integralfnoc$value+1/k1)
f0 <- cPDF/k1

PDFf1 <- Vectorize(function(w)PDFf(w, ca1, k1, MAP1, wL, cPDF))

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 4), mfrow=c(1,1))
curve(PDFf1, wL, 1, xlim=c(0, 1), cex.lab=1.3)
