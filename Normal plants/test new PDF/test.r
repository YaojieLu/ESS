
options(digits=20)
source("Normal plants/test new PDF/tf.r")

# function for all the derived variables
dvsf <- function(ca, k, MAP){
  
  ESSBf1 <- Vectorize(function(w)ESSBf(w, ca))
  wL <- uniroot(ESSBf1, c(0.12, 1), tol=.Machine$double.eps)$root
  integralfnoc <- integralfnocf(ca, k, MAP, wL)
  gamma <- 1/((MAP/365/k)/1000)*0.5
  cPDF <- 1/(integralfnoc$value+1/k*exp(-gamma*wL))
  #browser()
  fL <- cPDF/k*exp(-gamma*wL)
  averA <- averAf(ca, k, MAP, wL, cPDF)
  averm <- avermf(ca, k, MAP, wL, cPDF)
  averE <- averEf(ca, k, MAP, wL, cPDF)
  averwp1 <- averwp1f(ca, k, MAP, wL, cPDF)
  averw <- averwp1$value+fL*wL
  avercica <- avercicaf(ca, k, MAP, wL, cPDF)
  EMAP <- averE$value*500*365/MAP
  res <- c(integralfnoc$abs.error/integralfnoc$value*100, averwp1$abs.error/averwp1$value*100, averE$abs.error/averE$value*100, averA$abs.error/averA$value*100, averm$abs.error/averm$value*100, avercica$abs.error/avercica$value*100, fL, averw, averE$value, averA$value, averm$value, avercica$value, EMAP)
  return(res)
}

#environmental conditions
ca <- c(400)  # Atmospheric CO2 concentration (ppm)
k <- c(0.025, 0.05, 0.1) # Rainfall frequency (per day)
MAP <- seq(100, 3500, by=100) # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=13)

# Run every parameter combination
for(i in 1:nrow(env)){
  begin <- proc.time()
  dvs[i,] <- dvsf(env[i, 1], env[i, 2], env[i, 3])
  end <- proc.time()
  message(sprintf("%s/%s completed in %.2f min",i, nrow(env), (end[3]-begin[3])/60))
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "integralfnoc rel error (%)", "averwp1 rel error (%)", "averE rel error (%)", "averA rel error (%)", "averm rel error (%)", "avercica rel error (%)", "fL", "averw", "averE", "averA", "averm", "avercica", "E/MAP") 

write.csv(res, "Normal plants/test new PDF/test.csv", row.names = FALSE)
