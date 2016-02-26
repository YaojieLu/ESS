
library(rgl)
gs1 <- function(w, w0,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05){
  h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)}
gs2 <- Vectorize(function(w, w0)ifelse(w>w0, gs1(w, w0), 0))

# data
w <- seq(0.6, 1, by=(1-0.6)/100)
data <- matrix(data=numeric(), nrow=101, ncol=101)
for(i in 1:101){
  data[, i] <- sapply(w, gs2, w0=w[i])
}

# Optimize
gs3 <- function(w0, w)-gs2(w, w0)
gsmax <- Vectorize(function(x)optimize(gs3, c(0, x), w=x)$minimum)

# Figure
terrain3d(x=w, y=w, z=data, color="blue", xlab="w", ylab="w0")
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(gsmax, 0.5, 1, ylim=c(0.5, 1))
abline(a=0, b=1, col="red")

#gs4 <- function(w0)-gs2(w=1, w0)
#curve(gs4, 0, 1)
