####################################################################
# Fused Sparse Group Lasso Simulation Study
# Make Supplementary Figure S1: Bias-Variance Decomposition of MSE
####################################################################
# This script is used to report simulation results
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# 9 csv files that are simulation study results output from the script
# s02_FSGLassoSimulation_TestDataPredictionError.R, named as follows: 
# compaggbv.csv, partaggbv.csv, compdistbv.csv,
# spcompaggbv.csv, sppartaggbv.csv, spcompdistbv.csv, 
# exspbv.csv, misspbv.csv, misspspbv.csv
### OUTPUTS:
# pdf figure 'biasvar.pdf
####################################################################

########################################
# LOAD DATASETS
########################################
# set the working and data directories
# setwd()
datadir <- paste0(getwd(), '/BiasVarDecomp/')
compaggbv <- read.csv(paste0(datadir, 'compaggbv.csv'))
spcompaggbv <- read.csv(paste0(datadir, 'spcompaggbv.csv'))
partaggbv <- read.csv(paste0(datadir, 'partaggbv.csv'))
sppartaggbv <- read.csv(paste0(datadir, 'sppartaggbv.csv'))
compdistbv <- read.csv(paste0(datadir, 'compdistbv.csv'))
spcompdistbv <- read.csv(paste0(datadir, 'spcompdistbv.csv'))
exspbv <- read.csv(paste0(datadir, 'exspbv.csv'))
misspbv <- read.csv(paste0(datadir, 'misspbv.csv'))
misspspbv <- read.csv(paste0(datadir, 'misspspbv.csv'))

########################################
# function to plot means over all 100 obs 
# xorder = 1 means alpha changes fastest, gamma slowest
# xorder = 2 means alpha changes slowest, gamma fastest
########################################
bvplot <- function(bvdata, title, xorder=1, ylimits=NULL){
  means <- colMeans(bvdata[,3:86])
  if (xorder==1){
    vars <- means[(1:21*4) - 3][c(1,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,5,9,13,17,21)]
    bias2 <- means[(1:21*4) - 1][c(1,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,5,9,13,17,21)]
    mse <- means[(1:21*4)][c(1,2,6,10,14,18,3,7,11,15,19,4,8,12,16,20,5,9,13,17,21)]
  } else if (xorder==2){
    vars <- means[(1:21*4) - 3]
    bias2 <- means[(1:21*4) - 1]
    mse <- means[(1:21*4)]
  }
  xval <- 1:21
  par(mar=c(5, 3, 3, 1), xpd=TRUE)
  if (!is.null(ylimits)){
    plot(xval, mse, ylim=ylimits, type='l', main=title, ylab='', xlab='', xaxt='n', cex.main=0.9, adj=0, cex.axis=0.9)
  } else {
    plot(xval, mse, ylim=c(0, max(mse) + 1), type='l', main=title, ylab='', xlab='', xaxt='n', cex.main=0.9, adj=0, cex.axis=0.9)
  }
  if (xorder==1){
    axis(side=1, at=c(1,2,6,7,11,12,16,17,21), labels=c('(0, 0)', '(0, 0.2)', '(1, 0.2)', '(0, 0.5)', '(1, 0.5)', '(0, 0.8)', '(1, 0.8)', '(0, 1)', '(1, 1)'), las=2, cex.axis=0.5)
  } else if (xorder==2){
    axis(side=1, at=c(1,5,6,9,10,13,14,17,18,21), labels=c('(0, 0)', '(0, 1)', '(0.2, 0.2)', '(0.2, 1)', '(0.5, 0.2)', '(0.5, 1)', '(0.8, 0.2)', '(0.8, 1)', '(1, 0.2)', '(1, 1)'), las=2, cex.axis=0.5)
  }
  title(xlab=expression(bold(paste("(", alpha, ", ", gamma,") values"))), line=3, cex.lab=0.9)
  lines(xval, vars, type='l', col='blue')
  lines(xval, bias2, type='l', col='red')
  legend('topright', legend=c('MSE', 'Minimum MSE', 'Squared Bias', 'Variance'), lty=c(1,2,1,1), col=c('black', 'black', 'red', 'blue'), bty='n', inset=c(-0.10,-0.27), cex=0.5)
  par(xpd=FALSE)
  if (xorder==1){
    polygon(c(1, 2, 2, 1), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(6, 7, 7, 6), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(11, 12, 12, 11), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(16, 17, 17, 16), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    points(1, mse[1], col='black', pch=15, cex=0.7)
    points(1, vars[1], col='blue', pch=15, cex=0.7)
    points(1, bias2[1], col='red', pch=15, cex=0.7)
  } else if (xorder==2){
    polygon(c(5, 6, 6, 5), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(9, 10, 10, 9), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(13, 14, 14, 13), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
    polygon(c(17, 18, 18, 17), c(-10, -10, max(mse)+10, max(mse)+10), col='white', border='white')
  }
  abline(v=which.min(mse), lty=2)
  box()
}

########################################
# make plot
########################################
pdf('biasvar.pdf', width=6, height=7.5)
par(mfrow=c(3,3), pty='s', las=1)
bvplot(compaggbv, 'True Coefficients 1A', ylimits=c(0, 255))
bvplot(partaggbv, 'True Coefficients 2B', ylimits=c(0, 255))
bvplot(compdistbv, 'True Coefficients 3C', ylimits=c(0, 255))
bvplot(spcompaggbv, 'True Coefficients 4A', ylimits=c(0, 155))
bvplot(sppartaggbv, 'True Coefficients 5B', ylimits=c(0, 155))
bvplot(spcompdistbv, 'True Coefficients 6C', ylimits=c(0, 155))
bvplot(exspbv, 'True Coefficients 7B')
bvplot(misspbv, 'True Coefficients 8B', ylimits=c(0, 225))
bvplot(misspspbv, 'True Coefficients 9B', ylimits=c(0, 225))
dev.off()
