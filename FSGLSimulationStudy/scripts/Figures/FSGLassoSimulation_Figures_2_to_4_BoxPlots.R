####################################################################
# Fused Sparse Group Lasso Simulation Study
# Make Figures 2 to 4: Box Plots Mean CVE, MSE beta, Prediction MSE
####################################################################
# This script is used to report simulation results
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# 9 csv files that are simulation study results output from the script
# s02_FSGLassoSimulation_TestDataPredictionError.R, named as follows: 
# compaggplot.csv, partaggplot.csv, compdistplot.csv,
# spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, 
# exspplot.csv, misspplot.csv, misspspplot.csv
### OUTPUTS:
# 3 pdf figures, 'boxplots_1to3.pdf', 'boxplots_4to6.pdf', 'boxplots_7to9.pdf'
####################################################################

########################################
# LOAD DATASETS
########################################
# set the working and data directories
# setwd()
datadir1 <- paste0(getwd(), '/DataForPlots/')
compagg <- read.csv(paste0(datadir1, 'compaggplot.csv'))
partagg <- read.csv(paste0(datadir1, 'partaggplot.csv'))
compdist <- read.csv(paste0(datadir1, 'compdistplot.csv'))
spcompagg <- read.csv(paste0(datadir1, 'spcompaggplot.csv'))
sppartagg <- read.csv(paste0(datadir1, 'sppartaggplot.csv'))
spcompdist <- read.csv(paste0(datadir1, 'spcompdistplot.csv'))
exsp <- read.csv(paste0(datadir1, 'exspplot.csv'))
missp <- read.csv(paste0(datadir1, 'misspplot.csv'))
misspsp <- read.csv(paste0(datadir1, 'misspspplot.csv'))

########################################
# function to make boxplot over 100 sim 
# for each of 21 (alpha, gamma) combinations
# xorder = 1 means alpha changes fastest, gamma slowest
# xorder = 2 means alpha changes slowest, gamma fastest
########################################
bxplot <- function(simdata, variable, title, xorder=1, xaxslab=FALSE, color='grey', margins=c(1, 4, 3, 1), ylabel, ylimits=NULL){
  if (xorder==1){
    simdata$gp <- as.numeric(interaction(simdata$alpha, simdata$gamma, lex.order=FALSE))
  } else if (xorder==2){
    simdata$gp <- as.numeric(interaction(simdata$alpha, simdata$gamma, lex.order=TRUE))
  }
  simdata <- simdata[order(simdata$gp),]
  par(mar=margins, xpd=TRUE)
  if (!is.null(ylimits)){
    boxplot(eval(parse(text=variable)) ~ gp, data=simdata, main=title, ylab='', xlab='', xaxt='n', cex.main=1, adj=0, cex.axis=0.9, col=color, medlwd=1, ylim=ylimits)
  } else {
    boxplot(eval(parse(text=variable)) ~ gp, data=simdata, main=title, ylab='', xlab='', xaxt='n', cex.main=1, adj=0, cex.axis=0.9, col=color, medlwd=1)
  }
  if (xaxslab==TRUE){
    if (xorder==1){
      axis(side=1, at=c(1,2,6,7,11,12,16,17,21), labels=c('(0, 0)', '(0, 0.2)', '(1, 0.2)', '(0, 0.5)', '(1, 0.5)', '(0, 0.8)', '(1, 0.8)', '(0, 1)', '(1, 1)'), las=2, cex.axis=0.5)
    } else if (xorder==2){
      axis(side=1, at=c(1,5,6,9,10,13,14,17,18,21), labels=c('', '', '', '', '', '', '', '', '', ''), las=2, cex.axis=0.5)
    }
  } else if (xaxslab==FALSE){
    if(xorder==1){
      axis(side=1, at=c(1,2,6,7,11,12,16,17,21), labels=c('', '', '', '', '', '', '', '', ''), las=2, cex.axis=0.5)
    }
  }
  title(ylab=ylabel, cex.lab=1, font.lab=2, line=2.25)
  par(xpd=FALSE)
  # plot vertical line(s) at alpha gamma with most frequent low error
  lowgp <- c()
  for (i in 1:100){
    simdataseed <- simdata[simdata$seed==i,]
    lowgp[i] <- simdataseed$gp[which.min(simdataseed[,variable])]
  }
  lowest <- as.numeric(names(which(table(lowgp) == max(table(lowgp)))))
  abline(v=sapply(lowest, function(x) which(unique(simdata$gp) == x)), lty=1, col='red')
}

########################################
# True Coefficients 1 to 3
########################################
pdf('boxplots_1to3.pdf', width=6.5, height=6)
par(mfrow=c(3,3), pty='s', oma=c(4, 3, 3, 0) + 0.1, las=1)
# ROW 1
bxplot(compagg, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,310))
mtext("Cross-Validation MSE", side=2, line=3, font=2, cex=0.75, las=0)
mtext('True Coefficients 1A\nCompletely Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(partagg, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,310))
mtext('True Coefficients 2B\nPartially Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(compdist, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,310))
mtext('True Coefficients 3C\nCompletely Distributed', side=3, line=1, font=2, cex=0.75)
# ROW 2
bxplot(compagg, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.9))
mtext(expression(bold(paste("MSE(", widehat(beta), ")"))), side=2, line=3, font=2, cex=0.75, las=0)
bxplot(partagg, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.9))
bxplot(compdist, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.9))
# ROW 3
bxplot(compagg, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits = c(0,460))
mtext(expression(bold(paste("MSE(", widehat(y)[test], ")"))), side=2, line=3, font=2, cex=0.75, las=0)
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(partagg, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits = c(0,460))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(compdist, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits = c(0,460))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
dev.off()

########################################
# True Coefficients 4 to 6
########################################
pdf('boxplots_4to6.pdf', width=6.5, height=6)
par(mfrow=c(3,3), pty='s', oma=c(4, 3, 3, 0) + 0.1, las=1)
# ROW 1
bxplot(spcompagg, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,190))
mtext("Cross-Validation MSE", side=2, line=3, font=2, cex=0.75, las=0)
mtext('True Coefficients 4A\nSparse Completely Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(sppartagg, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,190))
mtext('True Coefficients 5B\nSparse Partially Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(spcompdist, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,190))
mtext('True Coefficients 6C\nSparse Completely Distributed', side=3, line=1, font=2, cex=0.75)
# ROW 2
bxplot(spcompagg, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 0.46))
mtext(expression(bold(paste("MSE(", widehat(beta), ")"))), side=2, line=3, font=2, cex=0.75, las=0)
bxplot(sppartagg, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 0.46))
bxplot(spcompdist, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 0.46))
# ROW 3
bxplot(spcompagg, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 255))
mtext(expression(bold(paste("MSE(", widehat(y)[test], ")"))), side=2, line=3, font=2, cex=0.75, las=0)
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(sppartagg, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 255))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(spcompdist, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 255))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
dev.off()

########################################
# True Coefficients 7 to 9
########################################
pdf('boxplots_7to9.pdf', width=6.5, height=6.25)
par(mfrow=c(3,3), pty='s', oma=c(4, 3, 4, 0) + 0.1, las=1)
# ROW 1
bxplot(exsp, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,40))
mtext("Cross-Validation MSE", side=2, line=3, font=2, cex=0.75, las=0)
mtext('True Coefficients 7B\nExtra Sparse\nPartially Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(missp, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 325))
mtext('True Coefficients 8B\nMisspecified\nPartially Aggregated', side=3, line=1, font=2, cex=0.75)
bxplot(misspsp, variable='mean.cve', title='', xorder=1, color='lightblue', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0, 325))
mtext('True Coefficients 9B\nMisspecified Sparse\nPartially Aggregated', side=3, line=1, font=2, cex=0.75)
# ROW 2
bxplot(exsp, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.15))
mtext(expression(bold(paste("MSE(", widehat(beta), ")"))), side=2, line=3, font=2, cex=0.75, las=0)
bxplot(missp, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.68))
bxplot(misspsp, variable='mse.betas', title='', xorder=1, color='pink', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,0.68))
# ROW 3
bxplot(exsp, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,80))
mtext(expression(bold(paste("MSE(", widehat(y)[test], ")"))), side=2, line=3, font=2, cex=0.75, las=0)
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(missp, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,360))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
bxplot(misspsp, variable='mse.y.test', title='', xorder=1, xaxslab=TRUE, color='lightgreen', margins=c(0.5, 1, 0.5, 0), ylabel='', ylimits=c(0,360))
mtext(expression(bold(paste("(", alpha, ", ", gamma,") values"))), side=1, line=3, font=2, cex=0.75)
dev.off()



