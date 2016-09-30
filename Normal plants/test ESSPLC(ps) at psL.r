
source("Normal plants/tf gs(ps).r")

# Initializing
ca <- 400
SA <- seq(1, 100, by=1)
data <- data.frame(wL=numeric(length(SA)), psL=numeric(length(SA)), PLCwL=numeric(length(SA)))

# Sensitivity analysis
for(i in 1:length(SA)){
  h3 <- SA[i]
  data[i, 1] <- uniroot(ESSBf, c(0.12, 1), tol=.Machine$double.eps, ca=ca)$root
  data[i, 2] <- psf(data[i, 1])
  data[i, 3] <- ESSPLCf(data[i, 1], ca)*100
}

#Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1), mfrow=c(1, 1))
plot(data[, 2], data[, 3],
      type="l", xaxt="n", yaxt="n", xlab=NA, ylab=NA,
      xlim=c(-4, 0), ylim=c(0, 100), cex.lab=1.3)

axis(1, xlim=c(-4, 0), pos=0, lwd=2)
mtext(expression(psi[s]~(MPa)),side=1,line=2.3, cex=1.3)
axis(2, ylim=c(0, 100), pos=-4, lwd=2)
mtext("PLC (%)",side=2,line=1.8, cex=1.3)