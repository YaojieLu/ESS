
library(inline)
source("Two - functions.r")
options(digits=22)

hk <- 10
L <- 8
ca <- 400
k <- 0.05
MAP <- 1825

# stationary point
netA <- function(g,
                 ca=eval(ca), hk=eval(hk),
                 Vcmax=50, Km=703, cp=30, Rd=1)-(0.5*(Vcmax+(Km+ca)*g-Rd-(Vcmax^2+2*Vcmax*(Km-ca+2*cp)*g+((ca+Km)*g+Rd)^2-2*Rd*Vcmax)^0.5)-hk*g)
gX <- optimize(netA, c(0, 1))
gX

# run models
case1 <- source("Two - SANN.r")
case1

curve()