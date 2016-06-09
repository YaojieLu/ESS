
options(digits=20)
source("Normal plants/Functions.r")

# function for all the derived variables
dvsf <- function(ca, k, MAP){
  
  ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
  wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
  integralfnoc <- integralfnocf(ca, k, MAP, wL)
  cPDF <- 1/(integralfnoc$value+1/k)
  #browser()
  f0 <- cPDF/k
  averB <- averBf(ca, k, MAP, wL, cPDF)
  averE <- averEf(ca, k, MAP, wL, cPDF)
  averwp1 <- averwp1f(ca, k, MAP, wL, cPDF)
  averw <- averwp1$value+f0*wL
  
  res <- c(integralfnoc$abs.error/integralfnoc$value*100, averwp1$abs.error/averwp1$value*100, averE$abs.error/averE$value*100, averB$abs.error/averB$value*100, f0, averw, averE$value, averB$value)
  return(res)
}

#environmental conditions
ca <- c(400)  # Atmospheric CO2 concentration (ppm)
k <- c(0.025, 0.05, 0.1) # Rainfall frequency (per day)
MAP <- seq(100, 3000, by=100) # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=8)

# Run every parameter combination
for(i in 1:nrow(env)){
  begin <- proc.time()
  dvs[i,] <- dvsf(env[i, 1], env[i, 2], env[i, 3])
  end <- proc.time()
  message(sprintf("%s/%s completed in %.2f min",i, nrow(env), (end[3]-begin[3])/60))
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "integralfnoc rel error (%)", "averwp1 rel error (%)", "averE rel error (%)", "averB rel error (%)", "f0", "averw", "averE", "averB") 

write.csv(res, "Normal plants/Derived variables.csv", row.names = FALSE)
