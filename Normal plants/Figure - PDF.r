
source("Normal plants/Functions.r")
# Initializing
ca <- 400
k <- 0.05
MAP <- 2100

ESSf1 <- Vectorize(function(w)ESSf(w, ca))
ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
integralfnoc <- integralfnocf(ca, k, MAP, wL)
cPDF <- 1/(integralfnoc$value+1/k)
PDFf1 <- Vectorize(function(w)PDFf(w, ca, k, MAP, wL, cPDF))

# Result
x <- seq(wL, 1, by=0.008)
y <- PDFf1(x)

# Figures
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 1.5, 3, 1.5), mfrow=c(1,1))
plot(x, y, type="l",
     main=as.expression(bquote(MAP==.(MAP))),
     xlab=expression(italic(w)), ylab=NA, yaxt="n",
     xlim=c(0, 1), cex.lab=1.3)
