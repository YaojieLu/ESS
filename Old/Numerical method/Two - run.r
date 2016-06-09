
library(inline)
source("Two - functions.r")
options(digits=22)
# envirnonmental conditions
hk <- 10
ca <- 400
k <- 0.05
MAP <- 1095
# stationary point
netA <- function(g,
                 ca1=ca, hk1=hk,
                 Vcmax=50, Km=703, cp=30, Rd=1)-(0.5*(Vcmax+(Km+ca1)*g-Rd-(Vcmax^2+2*Vcmax*(Km-ca1+2*cp)*g+((ca1+Km)*g+Rd)^2-2*Rd*Vcmax)^0.5)-hk1*g)
gX <- optimize(netA, c(0, 1))
gX
# run models
L <- 5 
case1 <- source("Two - SANN.r")
case1
