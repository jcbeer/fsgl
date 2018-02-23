####################################################################
# Fused Sparse Group Lasso Simulation Study
# Make Supplemental Tables S2 to S10
# for summarizing simulation results in detail
####################################################################
# This script is used to analyze the simulation results
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
## original simulation statistics plus test set MSE
# compaggplot.csv, partaggplot.csv, compdistplot.csv,
# spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, 
# exspplot.csv, misspplot.csv, misspspplot.csv
## bias - variance results
# compaggcv.csv, partaggcv.csv, compdistcv.csv,
# spcompaggcv.csv, sppartaggcv.csv, spcompdistcv.csv, 
# exspcv.csv, misspcv.csv, misspspcv.csv
### OUTPUTS:
# text printed to console that can be cut and pasted into LaTeX document
####################################################################

library(xtable) # creates LaTeX tables

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

# setwd()
datadir2 <- paste0(getwd(), '/BiasVarDecomp/')
compaggbv <- read.csv(paste0(datadir2, 'compaggbv.csv'))
spcompaggbv <- read.csv(paste0(datadir2, 'spcompaggbv.csv'))
partaggbv <- read.csv(paste0(datadir2, 'partaggbv.csv'))
sppartaggbv <- read.csv(paste0(datadir2, 'sppartaggbv.csv'))
compdistbv <- read.csv(paste0(datadir2, 'compdistbv.csv'))
spcompdistbv <- read.csv(paste0(datadir2, 'spcompdistbv.csv'))
exspbv <- read.csv(paste0(datadir2, 'exspbv.csv'))
misspbv <- read.csv(paste0(datadir2, 'misspbv.csv'))
misspspbv <- read.csv(paste0(datadir2, 'misspspbv.csv'))

########################################
# function to make table
########################################
maketable <- function(simdata, simdatabv){
  simdata$alphagamma <- as.numeric(interaction(simdata$alpha, simdata$gamma))
  simdata <- simdata[order(simdata$seed, simdata$gamma, simdata$alpha),]
  # frequency of minimum CVE, MSE Beta, MSE Y test out of 100 simulations
  freq.min <- matrix(nrow=100, ncol=3)
  for(seed in 1:100){
    simdataseed <- simdata[simdata$seed==seed,]
    freq.min[seed,] <- apply(simdataseed[,c('mean.cve', 'mse.betas', 'mse.y.test')], 2, function(x) which.min(x))
  }
  # add alpha, gamma columns to frequency columns
  columns1to5 <- cbind(simdata$alpha[1:21], simdata$gamma[1:21], table(factor(freq.min[,1], levels=1:21)), table(factor(freq.min[,2], levels=1:21)), table(factor(freq.min[,3], levels=1:21)))
  # mean (sd) for mean CVE, cd CVE, MSE Beta, MSE Y test, Optimal lambda
  means <- matrix(nrow=21, ncol=5)
  sds <- matrix(nrow=21, ncol=5)
  for(alphagamma in simdata$alphagamma[1:21]){
    index <- which(simdata$alphagamma[1:21]==alphagamma)
    simdataalphagamma <- simdata[simdata$alphagamma==alphagamma, c('mean.cve', 'sd.cve', 'mse.betas', 'mse.y.test', 'opt.lambda')]
    means[index,] <- colMeans(simdataalphagamma)
    sds[index,] <- apply(simdataalphagamma, 2, sd)
  }
  # mean (sd) of squared bias and variance for 100 test observations
  bias2 <- simdatabv[,c(5,9,25,41,57,73,13,29,45,61,77,17,33,49,65,81,21,37,53,69,85)]
  variance <- simdatabv[,c(5,9,25,41,57,73,13,29,45,61,77,17,33,49,65,81,21,37,53,69,85)-2]
  biasvarmeans <- cbind(colMeans(bias2), colMeans(variance))
  biasvarsds <- cbind(apply(bias2, 2, sd), apply(variance, 2, sd))
  # function to create character vector with mean (sd)
  meansd <- function(meanvalues, sdvalues, meandigits, sddigits){
    meanschar <- formatC(round(meanvalues, meandigits), format='f', digits=meandigits)
    sdschar <- formatC(round(sdvalues, sddigits), format='f', digits=sddigits)
    meansdchar <- cbind(meanschar, sdschar)
    output <- apply(meansdchar, 1, function(x) paste0(x[1], ' (', x[2], ')'))
    return(output)
  }
  meancve <- meansd(meanvalues=means[,1], sdvalues=sds[,1], meandigits=2, sddigits=2)
  sdcve <- meansd(meanvalues=means[,2], sdvalues=sds[,2], meandigits=2, sddigits=2)
  msebeta <- meansd(meanvalues=means[,3], sdvalues=sds[,3], meandigits=3, sddigits=3)
  mseytest <- meansd(meanvalues=means[,4], sdvalues=sds[,4], meandigits=2, sddigits=2)
  optlambda <- meansd(meanvalues=means[,5], sdvalues=sds[,5], meandigits=2, sddigits=2)
  sqbias <- meansd(meanvalues=biasvarmeans[,1], sdvalues=biasvarsds[,1], meandigits=2, sddigits=2)
  varnc <- meansd(meanvalues=biasvarmeans[,2], sdvalues=biasvarsds[,2], meandigits=2, sddigits=2)
  tablevalues <- cbind(formatC(columns1to5[,1], format='f', digits=1), formatC(columns1to5[,2], format='f', digits=1), columns1to5[,3:5], meancve, sdcve, msebeta, mseytest, optlambda, sqbias, varnc)
  tableheader <- c('alpha', 'gamma', 'Min CVE', 'MSE Beta', 'MSE Y Test', 'Mean CVE', 'SD CVE', 'MSE Beta', 'MSE Y Test', 'Optimal Lambda', 'Squared Bias', 'Variance')
  colnames(tablevalues) <- tableheader
  print(xtable(tablevalues, caption='EDIT CAPTION', label='EDIT LABEL', align=c('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c')), include.rownames=FALSE)
}

########################################
# make tables
########################################
maketable(compagg, compaggbv)
maketable(partagg, partaggbv)
maketable(compdist, compdistbv)
maketable(spcompagg, spcompaggbv)
maketable(sppartagg, sppartaggbv)
maketable(spcompdist, spcompdistbv)
maketable(exsp, exspbv)
maketable(missp, misspbv)
maketable(misspsp, misspspbv)
