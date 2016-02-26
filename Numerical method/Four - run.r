
library(inline)
source("Four - functions.r")
options(digits=22)
# envirnonmental conditions
hk <- 25
ca <- 400
k <- 0.025
MAP <- 365
# stationary point
netA <- function(g,
                 ca1=ca, hk1=hk,
                 Vcmax=50, Km=703, cp=30, Rd=1)-(0.5*(Vcmax+(Km+ca1)*g-Rd-(Vcmax^2+2*Vcmax*(Km-ca1+2*cp)*g+((ca1+Km)*g+Rd)^2-2*Rd*Vcmax)^0.5)-hk1*g)
gX <- optimize(netA, c(0, 1))
gX
# run models
L <- 7
case7 <- source("Four - SANN.r")
case7
L <- 6
case6 <- source("Four - SANN.r")
case6
L <- 5
case5 <- source("Four - SANN.r")
case5
L <- 4
case4 <- source("Four - SANN.r")
case4
L <- 3
case3 <- source("Four - SANN.r")
case3
L <- 2
case2 <- source("Four - SANN.r")
case2
L <- 1
case1 <- source("Four - SANN.r")
case1
L <- 0
case0 <- source("Four - SANN.r")
case0
