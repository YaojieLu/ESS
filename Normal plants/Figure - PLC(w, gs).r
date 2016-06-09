
source("Normal plants/Functions.r")

# Figure
# mf(w, gs)
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
with(data.frame(w=c(1)),
     {
       f <- Vectorize(function(gs)mf(w, gs)*10)
       curve(f, 0, gsmaxf(w),
             main=expression(PLC(italic(w,g[s]))),
             xlab=expression(italic(g[s])~(mol~m^-2~s^-1)), 
             ylab="PLC (%)",
             xlim=c(0, 0.5), ylim=c(0, 100), cex.lab=1.3, col="black")
     }
)
with(data.frame(w=c(0.22)),
     {
       f <- Vectorize(function(gs)mf(w, gs)*10)
       curve(f, 0, gsmaxf(w), col="red", add=T)
     }
)
with(data.frame(w=0.18),
     {
       f <- Vectorize(function(gs)mf(w, gs)*10)
       curve(f, 0, gsmaxf(w), col="blue", add=T)
     }
)
legend("topright", c("1", "0.22", "0.18"), lty=c(1, 1, 1), col=c("black", "red", "blue"), title=expression(italic(w)))
