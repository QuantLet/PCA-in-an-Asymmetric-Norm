[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **PEC_fMRI** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

﻿Name of Quantlet: "PEC_fMRI"

Published in: Principal Component Analysis in an Asymmetric Norm

Description: 'This is the application of  PrincipalExpectile algorithms to the fMRI data'

Keywords: expectile, pca, principal-components, quantile

See also: "PEC_algorithm_princdir.R , PEC_algorithm_topdown.R"

Author: Ngoc M. Tran, Maria Osipenko, Petra Burdejova

Submitted: "Sept 07 2016 by Petra Burdejova"

Input: "per_voxel_data.RData, Risk_Att.mat"

Output:  Plots of expectile components and scores

Example:  "PrincipalExpetile (PEC) components for the expectile levels 0.05, …, 0.95."

```

### R Code
```r

# load libraries
library(R.matlab)

man=c(1,2,3,4,5,6,7,8,9,11,13,14,16,17,18,19,21,23,24)

# vector  of risk attitudes (0...1 scale)
riskatt=readMat(file.path(dataDir,"Risk_Att.mat"))
risk_att <- as.vector(riskatt$Risk.Att[,2])

#plot risk attitudes
risk_att2 <- data.frame(cbind(man,risk_att))
risk_att2 <- risk_att2[order(-risk_att),]
#pdf(file = file.path(figDir,"risk_attitudes.pdf"),paper="special",height=5,width=9)
plot(risk_att2[,2], axes=FALSE, xlab="Individual ID", ylab="Risk attitude", ylim=c(-0.1,1.2), col="blue", pch=16, lwd=5, cex=1.5)
segments(1:19,rep(-0.25,19), x1 = 1:19, y1 = risk_att2[,2]-0.02, lty=2, col="grey", lwd=2)
axis(1,lwd=1,cex.axis=0.8,at=1:19,lab=risk_att2[,1])
axis(1,lwd=1,cex.axis=0.8,at=1:19,lab=risk_att2[,1]) # repeat to make bold
axis(2, ylim=c(-0.1,1.2))
box()
#dev.off()

# load data
load("/Users/PetraB/Github/PEC/output_voxel_data.RData")
# 2st coordinate = tau level
# 2st coordinate = type of output (1=1st scores, 2=2nd scores, 3= 1st component, 4 = 2nd component)
# 3rd coordinate = individual 
# 4th coordinate = brain area

# stimuluses
stim_list <- matrix (nrow=256, ncol=19)
for (j in 1:19){
  file <- paste0(100+man[j],"_output.txt")
  temp_data <- read.table(file.path(dataDir,file))	
  stim_list[,j] <- temp_data[,4]
} 
stim_list = round(stim_list/1000)
for (i in c(1:256)) {
  for (j in c(1:19)) {
    if (stim_list[i,j] == 0){stim_list[i,j] = max(stim_list[i,j-1],stim_list[i,j+1])}
  }}


stim_mtx = matrix ( nrow= 19, ncol = 6)
colnames(stim_mtx) <-  c ("PC1-aINS-L", "PC1-aINS-R", "PC1-DMPFC", "PC2-aINS-L", "PC2-aINS-R", "PC2-DMPFC")
rownames(stim_mtx) <- paste0("Ind_",man)
stimuluses = rbind(stim_list, stim_list+1,stim_list+2)

results = rep(0,19)

fitall <- vector(mode='list', length=19)

for (tau in c(1:19)){
    for (j in 1:19){ 
        for (brainpart_id  in 1:3) {
          print(tau,j, brainpart_id)
          stim_mtx[j,brainpart_id] = mean(output[[tau]][[1]][[j]][[brainpart_id]][stimuluses[,j]])
          stim_mtx[j,brainpart_id+3] = mean(output[[tau]][[2]][[j]][[brainpart_id]][stimuluses[,j]])
          stim_mtx[j,brainpart_id] = mean(output[[tau]][[1]][[j]][[brainpart_id]])
          stim_mtx[j,brainpart_id+3] = mean(output[[tau]][[2]][[j]][[brainpart_id]])
          }
          }
  results[tau] <- (summary(lm (risk_att ~ stim_mtx[,1:6]))$r.squared)
  fitall[[tau]] <- lm (risk_att ~ stim_mtx[,1:6])$fitted.values
} # of tau


#linear regresion for specific! tau=0.6
#pdf(file = file.path(figDir,"pec_fmri_fitted_riskatt.pdf"),paper="special",height=6,width=6)
fit <- unlist(fitall[12])
fin_model <-lm( risk_att  ~ fit ) 
plot( fit, risk_att,cex=2.5, xlab="Fitted value", ylab="Risk attitude", cex.lab=1.3, main="", ylim=c(-0.2,1.2),xlim=c(-0.25,0.9))
abline(fin_model, lwd=2, col="blue")
newx <- seq(min(fit), 0.9, 0.01)
lm_int <- predict(fin_model, newdata=data.frame(fit=newx), interval="confidence", level=0.99) 
lines(newx, lm_int[,2], lty=2, lwd=3, col="blue")
lines(newx, lm_int[,3], lty=2,lwd=3, col="blue") 
for (j in 1:19){text(fit[j], risk_att[j],man[j])}
#dev.off()


#plot loadings  no1 and no16(ie19)
par(mfrow=c(1,1))
pdf(file = file.path(figDir,"loadings_ind1_area1.pdf"),paper="special",height=3,width=9)
plot(output2[[0.6*20]][[1]][[1]][[1]], type="l", xlab="time", ylab="", ylim=c(-1500,1500))
dev.off()
pdf(file = file.path(figDir,"loadings_ind19_area1.pdf"),paper="special",height=3,width=9)
plot(output2[[0.6*20]][[1]][[16]][[1]], type="l", xlab="time", ylab="", ylim=c(-1500,1500))
dev.off()


# plot all R^2
plot(results, ylim=c(0,0.7), xlab="tau", ylab="expression(R^2)", type="l", xaxt="n")
points(results, ylim=c(0,0.7))
axis(1, at=c(1,10,19), lab=c(0.05,0.5,0.95))

```

automatically created on 2018-09-04