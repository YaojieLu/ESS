
# daily rainfall
AD1 <- fAD(n=500) # the number of simulated rainfall events
LS1 <- length(AD1)

# compile Frotran code
par0 <- c(20, 20, 20, 20)
daymodel <- make_f(AD1, LS1, par0, par0)
daymodelmono <- make_fmono(AD1, LS1, par0)

wp <- txtProgressBar(title ="Running simulation", label ="", min=0, max=100, initial=0, width=50, style=3)
p <- 0

par1 <- c(0.028027107,0.116094747,0.016460075,0.119243726)
result <- data.frame(a1=numeric(0), b1=numeric(0), c1=numeric(0), d1=numeric(0), avernetA1=numeric(0))
#result <- data.frame(a1=numeric(0), b1=numeric(0), c1=numeric(0), d1=numeric(0), avernetA1=numeric(0), a20=numeric(0), b20=numeric(0), c20=numeric(0), d20=numeric(0), avernetA20=numeric(0))
result[1, 1:5] <- c(par1, run_fmono(mod=daymodelmono, par=par1, AD=AD1, LS=LS1))

for(i in 1:100){
  result[i+1, 1:4] <- optim(result[i, 1:4], obj_wrapper, par1=result[i, 1:4], modmono=daymodelmono, modmult=daymodel, AD=AD1, LS=LS1, mult=1e10, Amin=L, method ="SANN")$par
  result[i+1, 5] <- run_fmono(mod=daymodelmono, par=result[i+1, 1:4], AD=AD1, LS=LS1)
  p <- p + 1
  setTxtProgressBar(wp, p)
}

close(wp)
result
