
library(plotBy)
source("Normal plants/Bootstrap predict_nls.R")
source("Normal plants/Bootstrap functions.R")
data <- read.csv("Normal plants/Derived variables.csv")

#data from Zhang et al. 2001
MAP <- c(
  1260, 1819, 645, 1639, 2978, 995, 1469, 1494, 1272, 1451, 1602, 596, 596, 701, 721, 788, 788, 836, 959, 1128, 1212, 1017, 2049, 2080, 1023, 1153, 1113, 1768, 2015, 1309, 704, 2648, 2425, 2851, 1983, 1800, 1750, 1854, 1611, 1193, 1207, 1333, 813, 813, 639, 457, 457, 1153, 1113, 2098, 2081, 1460, 1100, 579, 597
)
E <- c(
  1115, 1363, 622, 1190, 1483, 898, 1204, 1200, 1127, 1199, 1295, 552, 515, 683, 696, 701, 706, 714, 863, 879, 939, 956, 1328, 1291, 938, 859, 823, 1059, 1154, 1033, 631, 1311, 1440, 1481, 1265, 1145, 1019, 1016, 1063, 1168, 1122, 1180, 727, 726, 568, 437, 439, 860, 823, 1394, 1295, 990, 910, 467, 442
)
MAPE <- data.frame(MAP=MAP, E=E)
# the functional form follows Eq. 6 Zhang et al. 2001
nls1 <- nls(E ~ MAP*((1+w*E0/MAP)/(1+w*E0/MAP+MAP/E0)), start=list(w=3.822, E0=1128.347), data=MAPE)
p <- predict_nls(nls1, from=min(MAP), to=max(MAP), interval="confidence")

# Figures
windows(8, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, oma=c(2.5, 0, 0, 0), mar=c(1, 3.5, 1, 1), mfrow=c(2, 1))
Cols <- c("red", "darkgreen", "blue")

# average E
plot(MAP, E/365, panel.first={
  addpoly(p$x, p$lwr/365, p$upr/365)
  lines(p$x, p$pred/365)
},
type="n", cex.lab=1.3,
xlim=c(0, 3000), ylim=c(0, 8),
xlab=NA, ylab=NA, xaxt="n", yaxt="n")

plotBy(data$averE*500 ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 8),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols, add=T)

axis(1, xlim=c(0, 3000), pos=0, lwd=2, labels=FALSE)
axis(2, xlim=c(0, 8), pos=0, lwd=2)
mtext(expression(italic(bar(E))~"(mm per day)"), side=2, line=1.9, cex=1.3)
text(3000*0.03, 8*0.96, "a", cex=1.5)

# average E/MAP
plot(MAP, E/MAP, panel.first={
  addpoly(p$x, p$lwr/p$x, p$upr/p$x)
  lines(p$x, p$pred/p$x)
},
type="n", cex.lab=1.3,
xlim=c(0, 3000), ylim=c(0, 1),
xlab=NA, ylab=NA, xaxt="n", yaxt="n")

plotBy(data$averE*500*365/data$MAP ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 1),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols, add=T)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 1), pos=0, lwd=2)
mtext("MAP (mm)", side=1, line=2.2, cex=1.3)
mtext(expression(italic(bar(E))/MAP), side=2, line=1.9, cex=1.3)
legend("bottomleft", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Zhang~italic(et~al.)~2001),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1))
text(3000*0.03, 0.96, "b", cex=1.5)
