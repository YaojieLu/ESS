
# daily rainfall
AD1 <- fAD(n=500) # the number of simulated rainfall events
LS1 <- length(AD1)

# compile Frotran code
par0 <- c(20, 20)
daymodel <- make_f(AD1, LS1, par0, par0)
daymodelmono <- make_fmono(AD1, LS1, par0)

wp <- txtProgressBar(title = "Running simulation", label = "", min=0, max=44, initial=0, width=50, style=3)
p <- 0

par1 <- c(1, 1)
#result <- data.frame(a1= numeric(0), b1= numeric(0), avernetA1=numeric(0))
result <- data.frame(a1= numeric(0), b1= numeric(0), avernetA1=numeric(0), a20= numeric(0), b20= numeric(0), avernetA20=numeric(0))
result[1, 1:3] <- c(par1, run_fmono(mod=daymodelmono, par=par1, AD=AD1, LS=LS1))

for(i in 1:22){
  result[i+1, 1:2] <- optim(result[i, 1:2], obj_wrapper, par1=result[i, 1:2], modmono=daymodelmono, modmult=daymodel, AD=AD1, LS=LS1, mult=1e10, Amin=L, method = "SANN")$par
  result[i+1, 3] <- run_fmono(mod=daymodelmono, par=result[i+1, 1:2], AD=AD1, LS=LS1)
  p <- p + 1
  setTxtProgressBar(wp, p)
}

par20 <- c(20, 20)
result[1, 4:6] <- c(par20, run_fmono(mod=daymodelmono, par=par20, AD=AD1, LS=LS1))

for(i in 1:22){
  result[i+1, 4:5] <- optim(result[i, 4:5], obj_wrapper, par1=result[i, 4:5], modmono=daymodelmono, modmult=daymodel, AD=AD1, LS=LS1, mult=1e10, Amin=L, method = "SANN")$par
  result[i+1, 6] <- run_fmono(mod=daymodelmono, par=result[i+1, 4:5], AD=AD1, LS=LS1)
  p <- p + 1
  setTxtProgressBar(wp, p)
}

close(wp)
result
