
# daily rainfall (AD) function
fAD <- function(n, MDP=MAP/365, nZ=0.5, gamma=1/(MDP/k/1000)*nZ, rseed=1){
  
  set.seed(rseed)
  # DR[i] - the day when the ith rainfall event occurs
  DR <- floor(cumsum(rexp(n, k)))
  #AR[i] <- the amount of ith rainfall event
  AR <- rexp(n, gamma)
  #AD[i] - the amount of total rainfall on the ith day
  AD <- rep(0, length=tail(DR, n=1))
  for (i in 1:n){AD[DR[i]] <- AD[DR[i]]+AR[i]}
  return(AD)
}

# growth with shared water resource 
make_f <- function(AD, LS,
                   par1, par2,
                   Vcmax=50, Km=703, cp=30, Rd=1,VPD=0.02,p=43200,
                   LAI=1, nZ=0.5, a=1.6, l=1.8e-5, h=l*a*LAI/nZ*p){
  
  # to be safe
  AD <- as.numeric(AD)
  LS <- as.integer(LS)
  par1 <- as.numeric(par1)
  par2 <- as.numeric(par2)
  Vcmax <- as.numeric(Vcmax)
  Km <- as.numeric(Km)
  cp <- as.numeric(cp)
  Rd <- as.numeric(Rd)
  VPD <- as.numeric(VPD)
  p <- as.numeric(p)
  LAI <- as.numeric(LAI)
  nZ <- as.numeric(nZ)
  a <- as.numeric(a)
  l <- as.numeric(l)
  h <- as.numeric(h)
  
  #g[i] - average stomatal conductance on the ith day
  g1 <- vector("numeric", length=LS)
  g1[1] <- as.numeric(0)
  g2 <- vector("numeric", length=LS)
  g2[1] <- as.numeric(0)
  # A[i] - photosynthate production on the ith day; Î¼mol per leaf area per second
  A1 <- vector("numeric", length=LS)
  A1[1] <- as.numeric(0)
  A2 <- vector("numeric", length=LS)
  A2[1] <- as.numeric(0)
  # E[i] - the amount of soil water transpired on ith day; per day
  E1 <- vector("numeric", length=LS)
  E1[1] <- as.numeric(0)
  E2 <- vector("numeric", length=LS)
  E2[1] <- as.numeric(0)
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric", length=LS)
  RW[1] <- 0.5
  
  code <- "
  integer i
  
  do i=2, LS(1)
  
  g1(i)=par1(1)*RW(i-1)**par1(2)+par1(3)*RW(i-1)**par1(4)
  A1(i)=0.5*(Vcmax(1)+(Km(1)+ca(1))*g1(i)-Rd(1)-(Vcmax(1)**2.+  &
  2.*Vcmax(1)*(Km(1)-ca(1)+2.*cp(1))*g1(i)+((ca(1)+Km(1))*g1(i)+  &
  Rd(1))**2.-2.*Rd(1)*Vcmax(1))**0.5)-hk(1)*g1(i)
  E1(i)=h(1)*VPD(1)*g1(i)
  
  g2(i)=par2(1)*RW(i-1)**par2(2)+par2(3)*RW(i-1)**par2(4)
  A2(i)=0.5*(Vcmax(1)+(Km(1)+ca(1))*g2(i)-Rd(1)-(Vcmax(1)**2.+  &
  2.*Vcmax(1)*(Km(1)-ca(1)+2.*cp(1))*g2(i)+((ca(1)+Km(1))*g2(i)+  &
  Rd(1))**2.-2.*Rd(1)*Vcmax(1))**0.5)-hk(1)*g2(i)
  E2(i)=h(1)*VPD(1)*g2(i)
  
  RW(i)=RW(i-1)+AD(i)-E1(i)
  if(RW(i).gt.1.0)RW(i)=1.0
  if(RW(i).lt.0.0)RW(i)=0.0
  
  enddo
  "
  daymodel <- cfunction(signature(AD="numeric",
                                  LS="integer", 
                                  par1="numeric",
                                  par2="numeric",
                                  hk="numeric",
                                  ca="numeric",
                                  Vcmax ="numeric",
                                  Km="numeric",
                                  cp="numeric",
                                  Rd="numeric",
                                  VPD="numeric",
                                  h="numeric",
                                  RW="numeric",
                                  g1="numeric",
                                  A1="numeric",
                                  E1="numeric",
                                  g2="numeric",
                                  A2="numeric",
                                  E2="numeric"),
                        code, convention=".Fortran", language="F95")
  return(daymodel)
}

run_f <- function(mod,
                  AD, LS,
                  par1, par2,
                  Vcmax=50, Km=703, cp=30, Rd=1,VPD=0.02,p=43200,
                  LAI=1,nZ=0.5, a=1.6, l=1.8e-5, h=l*a*LAI/nZ*p){
  
  AD <- as.numeric(AD)
  LS <- as.integer(LS)
  par1 <- as.numeric(par1)
  par2 <- as.numeric(par2)
  Vcmax <- as.numeric(Vcmax)
  Km <- as.numeric(Km)
  cp <- as.numeric(cp)
  Rd <- as.numeric(Rd)
  VPD <- as.numeric(VPD)
  p <- as.numeric(p)
  LAI <- as.numeric(LAI)
  nZ <- as.numeric(nZ)
  a <- as.numeric(a)
  l <- as.numeric(l)
  h <- as.numeric(h)
  
  #g[i] - average stomatal conductance on the ith day
  g1 <- vector("numeric", length=LS)
  g1[1] <- as.numeric(0)
  g2 <- vector("numeric", length=LS)
  g2[1] <- as.numeric(0)
  # A[i] - photosynthate production on the ith day; Î¼mol per leaf area per second
  A1 <- vector("numeric", length=LS)
  A1[1] <- as.numeric(0)
  A2 <- vector("numeric", length=LS)
  A2[1] <- as.numeric(0)
  # E[i] - the amount of soil water transpired on ith day; per day
  E1 <- vector("numeric", length=LS)
  E1[1] <- as.numeric(0)
  E2 <- vector("numeric", length=LS)
  E2[1] <- as.numeric(0)
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric", length=LS)
  RW[1] <- 0.5
  
  z <- mod(AD, LS, par1, par2, hk, ca, Vcmax, Km, cp, Rd, VPD, h, RW, g1, A1, E1, g2, A2, E2)
  avernetA2 <- mean(z$A2)
  return(-avernetA2)
}

run_f_wrapper <- function(par, ...){run_f(par2=par, ...)}

# growth in monoculture 
make_fmono <- function(AD, LS,
                       par,
                       Vcmax=50, Km=703, cp=30, Rd=1,VPD=0.02,p=43200,
                       LAI=1, nZ=0.5, a=1.6, l=1.8e-5, h=l*a*LAI/nZ*p){
  
  # to be safe
  AD <- as.numeric(AD)
  LS <- as.integer(LS)
  par <- as.numeric(par)
  Vcmax <- as.numeric(Vcmax)
  Km <- as.numeric(Km)
  cp <- as.numeric(cp)
  Rd <- as.numeric(Rd)
  VPD <- as.numeric(VPD)
  p <- as.numeric(p)
  LAI <- as.numeric(LAI)
  nZ <- as.numeric(nZ)
  a <- as.numeric(a)
  l <- as.numeric(l)
  h <- as.numeric(h)
  
  #g[i] - average stomatal conductance on the ith day
  g <- vector("numeric", length=LS)
  g[1] <- as.numeric(0)
  # A[i] - photosynthate production on the ith day; Î¼mol per leaf area per second
  A <- vector("numeric", length=LS)
  A[1] <- as.numeric(0)
  # E[i] - the amount of soil water transpired on ith day; per day
  E <- vector("numeric", length=LS)
  E[1] <- as.numeric(0)
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric", length=LS)
  RW[1] <- 0.5
  
  code <- "
  integer i
  
  do i=2, LS(1)
  
  g(i)=par(1)*RW(i-1)**par(2)+par(3)*RW(i-1)**par(4)
  A(i)=0.5*(Vcmax(1)+(Km(1)+ca(1))*g(i)-Rd(1)-(Vcmax(1)**2.+  &
  2.*Vcmax(1)*(Km(1)-ca(1)+2.*cp(1))*g(i)+((ca(1)+Km(1))*g(i)+  &
  Rd(1))**2.-2.*Rd(1)*Vcmax(1))**0.5)-hk(1)*g(i)
  E(i)=h(1)*VPD(1)*g(i)
  
  RW(i)=RW(i-1)+AD(i)-E(i)
  if(RW(i).gt.1.0)RW(i)=1.0
  if(RW(i).lt.0.0)RW(i)=0.0
  
  enddo
  "
  daymodel <- cfunction(signature(AD="numeric",
                                  LS="integer", 
                                  par="numeric",
                                  hk="numeric",
                                  ca="numeric",
                                  Vcmax="numeric",
                                  Km="numeric",
                                  cp="numeric",
                                  Rd="numeric",
                                  VPD="numeric",
                                  h="numeric",
                                  RW="numeric",
                                  g="numeric",
                                  A="numeric",
                                  E="numeric"),
                        code, convention=".Fortran", language="F95")
  return(daymodel)
}

run_fmono <- function(mod,
                      AD, LS,
                      par,
                      Vcmax=50, Km=703, cp=30, Rd=1,VPD=0.02,p=43200,
                      LAI=1, nZ=0.5, a=1.6, l=1.8e-5, h=l*a*LAI/nZ*p){
  
  AD <- as.numeric(AD)
  LS <- as.integer(LS)
  par <- as.numeric(par)
  Vcmax <- as.numeric(Vcmax)
  Km <- as.numeric(Km)
  cp <- as.numeric(cp)
  Rd <- as.numeric(Rd)
  VPD <- as.numeric(VPD)
  p <- as.numeric(p)
  LAI <- as.numeric(LAI)
  nZ <- as.numeric(nZ)
  a <- as.numeric(a)
  l <- as.numeric(l)
  h <- as.numeric(h)
  
  #g[i] - average stomatal conductance on the ith day
  g <- vector("numeric", length=LS)
  g[1] <- as.numeric(0)
  # A[i] - photosynthate production on the ith day; Î¼mol per leaf area per second
  A <- vector("numeric", length=LS)
  A[1] <- as.numeric(0)
  # E[i] - the amount of soil water transpired on ith day; per day
  E <- vector("numeric", length=LS)
  E[1] <- as.numeric(0)
  # RW[i] - relative soil water content on ith day
  RW <- vector("numeric", length=LS)
  RW[1] <- 0.5
  
  z <- mod(AD, LS, par, hk, ca, Vcmax, Km, cp, Rd, VPD, h, RW, g, A, E)
  avernetA <- mean(z$A)
  return(avernetA)
}

obj <- function(par1, par2, modmono, modmult, mult, Amin, ...){
  
  avernetA1 <- run_fmono(par=par1, mod=modmono, ...)
  avernetA2 <- run_fmono(par=par2, mod=modmono, ...)
  advantage <- run_f(par1=par1, par2=par2, mod=modmult, ...)
  res <- advantage+100*((1-1/(10^(mult*(Amin-avernetA1))+1))+(1-1/(exp(mult*(Amin-avernetA2))+1)))
  return(res)
}

obj_wrapper <- function(par, ...){obj(par2=par, ...)}