
# Initialize
options(digits=2)
library(plotBy)
source("Normal plants/Bootstrap predict_nls.R")
source("Normal plants/Bootstrap functions.R")

# Data
# Model
data <- read.csv("Normal plants/Derived variables.csv")
# TNPP
dataTNPP <- read.csv("Normal plants/NPP.csv")
colnames(dataTNPP) <- c("MAP", "TNPP")
MAPTNPP <- dataTNPP$MAP
TNPP <- dataTNPP$TNPP
# E (Zhang et al. 2001)
MAPE <- c(
  1260, 1819, 645, 1639, 2978, 995, 1469, 1494, 1272, 1451, 1602, 596, 596, 701, 721, 788, 788, 836, 959, 1128, 1212, 1017, 2049, 2080, 1023, 1153, 1113, 1768, 2015, 1309, 704, 2648, 2425, 2851, 1983, 1800, 1750, 1854, 1611, 1193, 1207, 1333, 813, 813, 639, 457, 457, 1153, 1113, 2098, 2081, 1460, 1100, 579, 597
)
E <- c(
  1115, 1363, 622, 1190, 1483, 898, 1204, 1200, 1127, 1199, 1295, 552, 515, 683, 696, 701, 706, 714, 863, 879, 939, 956, 1328, 1291, 938, 859, 823, 1059, 1154, 1033, 631, 1311, 1440, 1481, 1265, 1145, 1019, 1016, 1063, 1168, 1122, 1180, 727, 726, 568, 437, 439, 860, 823, 1394, 1295, 990, 910, 467, 442
)
dataE <- data.frame(MAP=MAPE, E=E)

# Curve fit
# TNPP
nls1 <- nls(TNPP ~ a*MAPTNPP^b/exp(c*MAPTNPP), start=list(a=0.5346608789, b=1.0787801564, c= 0.0003562646), data=dataTNPP)
pTNPP <- predict_nls(nls1, from=min(MAPTNPP), to=max(MAPTNPP), interval="confidence")
pTNPP1 <- data.frame(x=pTNPP$x, pred=pTNPP$pred, lwr=pTNPP$lwr, upr=pTNPP$upr)
pTNPP2 <- subset(pTNPP1, x<=2500)
# E (The functional form below follows Eq. 6 in Zhang et al. 2001)
nls1 <- nls(E ~ MAPE*((1+w*E0/MAPE)/(1+w*E0/MAPE+MAPE/E0)), start=list(w=3.822, E0=1128.347), data=dataE)
pE <- predict_nls(nls1, from=min(MAPE), to=max(MAPE), interval="confidence")

# Figure
Cols <- c("red","darkgreen","blue")
windows(16, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, oma=c(0, 0, 0, 2.5), mar=c(3.5, 3.5, 1, 1), mfrow=c(1, 2))
# average E
plot(MAPE, E, panel.first={
  addpoly(pE$x, pE$lwr, pE$upr)
  lines(pE$x, pE$pred)},
  type="n", cex.lab=1.3,
  xlim=c(0, 3000), ylim=c(0, 3000),
  xlab="MAP (mm)", ylab=NA)

plotBy(data$averE*500*365 ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 3000),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols, add=T)

abline(a=0, b=1, lwd=1, lty=2)

mtext(expression(italic(bar(E))~"(mm per year)"), side=2, line=1.9, cex=1.3)
legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Zhang~italic(et~al.)~2001),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1))
text(3000*0.05, 3000*0.95, "a", cex=1.5)
box()

# TNPP
plot(MAPTNPP, TNPP, panel.first={
  addpoly(pTNPP2$x, pTNPP2$lwr, pTNPP2$upr)
  lines(pTNPP2$x, pTNPP2$pred)},
  type="n", xlab=NA, ylab=NA, axes=FALSE,
  xlim=c(0, 3000), ylim=c(0, 1420))

axis(4, ylim=c(0, 1200), pos=3000, lwd=2, at=c(0, 500, 1000))
mtext(expression(TNPP~(gC~m^-2~yr^-1)),side=4,line=2.5, cex=1.3)
rect(0, 0, 100, 1420, col="white", border=NA)

# average A
par(new=TRUE)
plotBy(data$averA ~ data$MAP | k, data=data,
       type='l', legend=FALSE, legendwhere="topleft",
       xlim=c(0, 3000), ylim=c(0, 15),
       xlab="MAP (mm)", ylab=NA,
       xaxt="n", yaxt="n",
       cex.lab=1.3, col=Cols)

axis(1, xlim=c(0, 3000), pos=0, lwd=2)
axis(2, xlim=c(0, 15), pos=0, lwd=2, at=c(0, 5, 10, 15))
mtext(expression(bar(italic(A[N]))~(mu*mol~m^-2~s^-1)), side=2, line=1.8, cex=1.3)#*" or "*bar(italic(m))
legend("bottomright", expression(italic(k==0.025), italic(k==0.05), italic(k==0.1), Del~Grosso~italic(et~al.)~2008),
       col = c("red","darkgreen","blue", "black"), lty=c(1, 1, 1, 1), lwd=c(2, 2, 2, 2))
text(3000*0.05, 15*0.95, "b", cex=1.5)
box()

dev.copy2pdf(file = "Normal plants//Figures/Figure 5.pdf")
