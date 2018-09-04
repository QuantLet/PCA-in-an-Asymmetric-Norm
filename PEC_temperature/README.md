[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **PEC_temperature** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿Name of Quantlet: "PEC_temperature"

Published in: 'Principal Component Analysis in an Asymmetric Norm'

Description: 'This is the application of TopDown, BottomUp and PrincipalExpectile algorithms to the Chinese weather data'

Keywords: expectile, pca, principal-components, quantile

See also: 'PEC_algorithm_princdir.R , PEC_algorithm_topdown.R'

Author: Ngoc M. Tran, Maria Osipenko, Petra Burdejova

Submitted:  Fri, July 29 2016 by Lining Yu

Input: geopos.csv, temperature.txt

Output:  Plots of expectile components and scores

Example: 'TopDown (TD), BottomUp (BUP), and PrincipalExpetile (PEC) components for the expectile levels 0.05, 0.5, and 0.95.'

```

### R Code
```r

rm(list = ls())
graphics.off()

#### enter main directory containning the functional libraries fpca.R and
#### PrincipalDirection.R below mainDir= 'H:/Asymmetric norm/Rcodes/application'

#### optional directory for figures figDir = 'H:/Asymmetric norm/talk/figures'

library("fda")
source(file.path(mainDir, "fpca.R"))
source(file.path(mainDir, "principalDirection.R"))

# -------------------download data
temp = read.table(file.path(mainDir, "temperature.txt"), header = T)

# -------------------substitute missing observations with a mean of neighbour
# values function to substitute some periods of NAs by the mean of neighbours
subna = function(x) {
  indna = which(is.na(x) == 1)
  if (length(indna) > 0) {
    x[indna] = (x[min(indna) - 1] + x[max(indna) + 1])/2
  }
  return(x)
}
temp = apply(temp, 2, subna)

# ------------------ function to remove the season Lambda_it
remseas = function(x) {
  time  = seq(1, length(x))
  reg = cbind(time, sin(2 * pi * time/365), cos(2 * pi * time/365), sin(pi * time/365), 
    cos(pi * time/365))
  fit = lm(x ~ reg)
  return(x - fit$fitted.values)
}
tempdes = apply(temp[, -1], 2, remseas)
# ------------------function to remove the AR(p) process
arcust = function(x) {
  mod = ar(x, aic = F, order.max = 10)
  return(mod$resid)
}
tempres = apply(tempdes, MARGIN = 2, arcust)
# -------------------function to average the residuals daywise X temp data
# vector, over - timespan
tempav = function(X, over) {
  Xmat = matrix(X, over, round(length(X)/365, 0))
  Xav = apply(Xmat, 1, mean, na.rm = T)
  return(Xav)
}
tempresan = apply(tempres, 2, tempav, over = 365)
# -------------------demean the residuals
man = apply(tempresan, 1, mean)
tempresan = tempresan - man
# --------------------smooth with Fourier series
daybasis = create.fourier.basis(c(0, 365), nbasis = 23, period = 365)  #nbasis=9
tempmfd = smooth.basis(day.5, tempresan, daybasis)$fd
tempsm = eval.fd(day.5, tempmfd)
# --------------------Compute TD, BUP PCs, and PEC for expectiles
set.seed(1234)
tau = list(0.05, 0.95, 0.5)
alpha = list(-0.45, 0.45, 0)

# PEC
print("PEC")
output.pec0 = pec.k(tempsm, nk = 2, alpha = -0.43)  #set as the starting value to speed up the computation
output.pec = lapply(alpha, pec.k, Y = tempsm, nk = 2, reset.tol = 50, lab.ini = -output.pec0[[4]])
# BUP
print("BUP")
output.bup = lapply(tau, bup, Y = tempsm, k = 2, basis.start = as.matrix(output.pec[[1]][[2]][, 
  1]))
# TD
print("TD")
output.td = lapply(X = tau, FUN = laws.main, Y = tempsm, nc = 2, max.iter = 30, 
  max.reset = 50, tol = 1e-10, mu = NA, preB = NA, mu.estimate = TRUE, basis.start = output.bup[[2]][[3]])

# ---------------------Project the bases to the span of the PEC(0.05) for visual
# comparison (utility functions listBasis and CompareBasis are in the functional
# library PrincipalDirection.R)
basis.pec = lapply(output.pec, listBasis, basis0 = output.pec[[1]][[2]], num = 2)
basis.bup = lapply(output.bup, listBasis, basis0 = output.pec[[1]][[2]], num = 3)
basis.td  = lapply(output.td, listBasis, basis0 = output.pec[[1]][[2]], num = 3)

# -------------------------- Table 5: explained
# variance----------------------------------
# ----------------------------------------------------------------------------------------
ess = matrix(NA, length(tau), 3)
for (k in 1:length(tau)) {
  tau0 = tau[[k]]
  sstot = sum(tempsm^2 * (sign(tempsm) * (tau0 - 0.5) + 0.5))
  ess[k, 1] = 1 - resid(tempsm, output.bup[[k]], tau = tau0, rss = T)/sstot
  ess[k, 2] = 1 - resid(tempsm, output.td[[k]], tau = tau0, rss = T)/sstot
  ess[k, 3] = 1 - resid.pec(tempsm, output.pec[[k]], tau = tau0, rss = T)/sstot
}
ess = data.frame(round(ess, 2))
rownames(ess) = c("tau=0.05", "tau=0.95", "tau=0.5")
colnames(ess) = c("BUP", "TD", "PEC")
# write.table(ess, 'rss.txt')
ess

# -------------------------------------Figures 3 and
# 4-------------------------------------
# -----------------------------------------------------------------------------------------
months = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")  # labels for the months
colh = (colh = c("darkgreen", 2, 4)  # color definition
)
line = 2  #line width
figDir = getwd()

# mean

plot.vec(tempsm, type = "l", lwd = 2, ylab = "", xlab = "", ylim = c(-1, 1), color = expgr)
lines(output.bup[[1]][[2]], lwd = line, col = colh[1])
lines(output.bup[[2]][[2]], lwd = line, col = colh[1], lty = 2)
lines(output.td[[1]][[2]], lwd = line, col = colh[2])
lines(output.td[[2]][[2]], lwd = line, col = colh[2], lty = 2)
lines(apply(resid.pec(tempsm, output.pec[[1]], fitted = T, tau = 0.05), 1, mean) + 
  apply(resid.pec(tempsm, output.pec[[1]], tau = 0.05), 1, expectile, alpha = 0.05 - 
    0.5), lwd = line, col = colh[3])
lines(apply(resid.pec(tempsm, output.pec[[2]], fitted = T, tau = 0.95), 1, mean) + 
  apply(resid.pec(tempsm, output.pec[[2]], tau = 0.95), 1, expectile, alpha = 0.95 - 
    0.5), lwd = line, col = colh[3], lty = 2)
axis(1, lwd = 1, cex.axis = 1.3, at = seq(1, 365, length = 12), lab = months)
axis(2, lwd = 1, cex.axis = 1.3, at = seq(-1, 1, by = 0.5))
legend(2, -0.8, c("BUP", "TD", "PEC"), lty = 1, lwd = line, horiz = TRUE, col = colh, 
  bty = "n", cex = 1)
legend(2, -0.6, c(expression(paste(tau, "=0.05")), expression(paste(tau, "=0.95"))), 
  lty = c(1, 2), lwd = line, horiz = TRUE, bty = "n", cex = 1)

# PCs

plot.new()
plot(basis.bup[[1]][, 1], type = "l", lwd = line, axes = FALSE, ylab = "", xlab = "", 
  ylim = c(-0.2, 0.2), col = colh[1])
lines(-basis.bup[[2]][, 1], lwd = line, col = colh[1], lty = 2)
lines(basis.td[[1]][, 1], lwd = line, col = colh[2])
lines(-basis.td[[2]][, 1], lwd = line, col = colh[2], lty = 2)
lines(basis.pec[[1]][, 1], lwd = line, col = colh[3])
lines(-basis.pec[[2]][, 1], lwd = line, col = colh[3], lty = 2)
lines(seq(1, 365), rep(0, 365), lwd = 1, lty = 2)
axis(1, lwd = 1, cex.axis = 1.3, at = seq(1, 365, length = 12), lab = months)
axis(2, lwd = 1, cex.axis = 1.3, at = seq(-0.2, 0.2, by = 0.1))
legend(2, -0.15, c("BUP", "TD", "PEC"), lty = 1, lwd = line, horiz = TRUE, col = colh, 
  bty = "n", cex = 1)
legend(2, -0.12, c(expression(paste(tau, "=0.05")), expression(paste(tau, "=0.95"))), 
  lty = c(1, 2), lwd = line, horiz = TRUE, bty = "n", cex = 1)



plot.new()
plot(basis.bup[[1]][, 2], type = "l", lwd = line, axes = FALSE, ylab = "", xlab = "", 
  ylim = c(-0.2, 0.2), col = colh[1])
lines(-basis.bup[[2]][, 2], lwd = line, col = colh[1], lty = 2)
lines(basis.td[[1]][, 2], lwd = line, col = colh[2])
lines(-basis.td[[2]][, 2], lwd = line, col = colh[2], lty = 2)
lines(basis.pec[[1]][, 2], lwd = line, col = colh[3])
lines(-basis.pec[[2]][, 2], lwd = line, col = colh[3], lty = 2)
lines(seq(1, 365), rep(0, 365), lwd = 1, lty = 2)
axis(1, lwd = 1, cex.axis = 1.3, at = seq(1, 365, length = 12), lab = months)
axis(2, lwd = 1, cex.axis = 1.3, at = seq(-0.2, 0.2, by = 0.1))
legend(2, -0.15, c("BUP", "TD", "PEC"), lty = 1, lwd = line, horiz = TRUE, col = colh, 
  bty = "n", cex = 1)
legend(2, -0.12, c(expression(paste(tau, "=0.05")), expression(paste(tau, "=0.95"))), 
  lty = c(1, 2), lwd = line, horiz = TRUE, bty = "n", cex = 1)


# ------------------------------------Scores on the map
libraries = c("maps","mapdata", "maptools", "classInt") 
lapply(libraries, function(x) if (!(x %in% installed.packages())) {   
install.packages(x) }) 
lapply(libraries, library, quietly = TRUE, character.only = TRUE)


scores_pec = rbind(basis.pec[[1]][, 1] %*% tempresan, basis.pec[[1]][, 2] %*% (tempresan - 
  apply(basis.pec[[1]][, 1] %*% tempresan, 2, "*", basis.pec[[1]][, 1])))
geo = read.table("geopos.csv", header = T, sep = ",")  #coordinates
geo[, 2] = substr(geo[, 2], 1, 4)

substrRight = function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}
geo[, 2] = as.numeric(substrRight(geo[, 2], 2))/60 + floor(as.numeric(geo[, 2])/100)  # transfer to GPS coordinates
geo[, 3] = as.numeric(substrRight(geo[, 3], 2))/60 + floor(as.numeric(geo[, 3])/100)
stn = as.numeric(substr(colnames(tempresan), 2, 6))

getcoor = function(stn, geo) {
  # asigns coordinates to stations
  coords = matrix((geo[geo[, 1] == stn, 2:3]), 1, 2)  #geo - table with station numbers and corresponding coords; stn - station number
  return(coords)
}
stn = as.matrix(stn)
coords = apply(stn, 1, getcoor, geo = geo)
coords = do.call(rbind, coords)

# set breaks for the 9 colors
brks = cut(scores_pec[1, ], breaks = c(-3, -1, 0, 1, 4))  #create levels for scores PEC1 for 0.05 exp
brks2 = cut(scores_pec[2, ], breaks = c(-3, -1, 0, 1, 4))  #create levels for scores PEC2 for 0.05 exp

# define colors
colors = c("darkgreen", "gold", "darkorange", "red")

plot.new()
map("china")
points(coords[10, 2], coords[10, 1])
points(coords[, 2], coords[, 1], col = colors[brks], pch = 16)
# add a legend
legend(x = 122, y = 35, levels(brks), fill = colors, bty = "n")

plot.new()
map("china")
points(coords[10, 2], coords[10, 1])
points(coords[, 2], coords[, 1], col = colors[brks2], pch = 16)
# add a legend
legend(x = 122, y = 35, levels(brks), fill = colors, bty = "n")
 

```

automatically created on 2018-09-04