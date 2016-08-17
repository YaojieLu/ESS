
############################### Scenario 4 ###############################
options(digits=20)
source("Normal plants/S4/Functions.r")

wmin <- 0.145

# function for all the derived variables
dvsf <- function(ca, k, MAP){
  integralfnoc <- integralfnocf(ca, k, MAP)
  cPDF <- 1/(integralfnoc$value)

  averA <- averAf(cPDF, ca, k, MAP)
  averE <- averEf(cPDF, ca, k, MAP)
  averw <- averwf(cPDF, ca, k, MAP)
  res <- c(integralfnoc$abs.error/integralfnoc$value*100, averw$abs.error/averw$value*100, averE$abs.error/averE$value*100, averA$abs.error/averA$value*100, averw$value, averE$value, averA$value)
  return(res)
}

#environmental conditions
ca <- c(400)  # Atmospheric CO2 concentration (ppm)
k <- c(0.025, 0.05, 0.1) # Rainfall frequency (per day)
MAP <- seq(100, 3000, by=100) # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=7)

# Run every parameter combination
for(i in 1:nrow(env)){
  begin <- proc.time()
  dvs[i,] <- dvsf(env[i, 1], env[i, 2], env[i, 3])
  end <- proc.time()
  message(sprintf("%s/%s completed in %.2f seconds",i, nrow(env), end[3]-begin[3]))
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "integralfnoc rel error (%)", "averw rel error (%)", "averE rel error (%)", "averA rel error (%)", "averw", "averE", "averA") 

write.csv(res, "Normal plants/S4/Derived variables.csv", row.names=FALSE)
