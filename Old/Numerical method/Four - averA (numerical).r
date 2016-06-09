
library(optimx)
options(digits=20)
f <- function(par,
              hk,
              ca, k, MAP,
              Vcmax=50, cp=30, Km=703, Rd=1, VPD=0.02, p=43200, LAI=1, nZ=0.5,
              a=1.6, l=1.8e-5, MDP=MAP/365, gamma=1/((MDP/k)/1000)*nZ, h=l*a*LAI/nZ*p){
  
  set.seed(1)
  
  #fAD: generating AD
  fAD <- function(n){
    # LD[i] - the length of the interval between (i-1)th and ith rainfall events
    LD <- rexp(n, k)
    # DR[i] - the day when the ith rainfall event occurs
    DR <- floor(cumsum(LD))
    #AR[i] <- the amount of ith rainfall event
    AR <- rexp(n, gamma)
    #AD[i] - the amount of rainfall on the ith day
    AD <- rep(0, length=tail(DR, 1)) # the length of simulation run
    for (i in 1:n){AD[DR[i]] <- AD[DR[i]]+AR[i]}
    return(AD)
  }
  
  AD <- fAD(100) # the number of rainfall events
  LS <- length(AD)
  
  gsw <- function(w)par[1]*w^par[2]+par[3]*w^par[4]
  
  #g[i] - average stomatal conductance on the ith day
  g <- vector("numeric")
  g[1] <- 0
  # A[i] - photosynthesis rate on the ith day; μmol per second
  A <- vector("numeric")
  A[1] <- 0
  # E[i] - the amount of soil water transpired on ith day  
  E <- vector("numeric")
  E[1] <- 0
  # I[i] - the amount (on relative scale) of rain water infiltrating the soil on ith day
  I <- vector("numeric")
  I[1] <- 0
  # dRWdt[i] - the change of relative soil water content on ith day
  dRWdt <- vector("numeric")
  dRWdt[1] <- 0
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric")
  RW[1] <- 0.5
  
  for(i in 2:LS){
    g[i] <- ifelse(RW[i-1]==0, 0, gsw(RW[i-1]))
    A[i] <- LAI*(0.5*(Vcmax+(Km+ca)*g[i]-Rd-(Vcmax^2+2*Vcmax*(Km-ca+2*cp)*g[i]+((ca+Km)*g[i]+Rd)^2-2*Rd*Vcmax)^0.5)-hk*g[i])
    E[i] <- h*VPD*g[i]
    #sum of infiltration and current RW must be smaller than 1
    I[i] <- min(AD[i], 1-RW[i-1])
    dRWdt[i] <- I[i]-E[i]
    x <- RW[i-1]+dRWdt[i]
    RW[i] <- ifelse(x<0, 0, x)
  }
  res <- mean(A)
  return(res)
}

# environmental conditions
ca1 <- 400
k1 <- 0.1
MAP1 <- 3650
hk1 <- 0
# Optimal
int1 <- c(1, 1, 1, 1)
opt <- optimx(int1, f, itnmax=5000, method="BFGS", control=list(maximize=T), ca=ca1, k=k1, MAP=MAP1, hk=hk1)
opt
# ESS
netA <- function(g,
                 ca=ca1, hk=hk1,
                 Vcmax=50, Km=703, cp=30, Rd=1)-(0.5*(Vcmax+(Km+ca)*g-Rd-(Vcmax^2+2*Vcmax*(Km-ca+2*cp)*g+((ca+Km)*g+Rd)^2-2*Rd*Vcmax)^0.5)-hk*g)
gX <- optimize(netA, c(0, 1))$minimum
gX
averAESS <- f(par=c(gX, 0, 0, 0), ca=ca1, k=k1, MAP=MAP1, hk=hk1)
averAESS
