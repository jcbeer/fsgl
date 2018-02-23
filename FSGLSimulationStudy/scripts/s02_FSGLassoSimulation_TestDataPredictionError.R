####################################################################
# Fused Sparse Group Lasso Simulation Study
# Calculate Test Set Prediction Error
# and Bias - Variance decomposition
####################################################################
# This script is used to analyze the simulation results
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# 9 csv files that are simulation study results output from the script
# s01_FSGLassoSimulation.R, named as follows: 
# compagg.csv, partagg.csv, compdist.csv,
# spcompagg.csv, sppartagg.csv, spcompdist.csv, 
# exsp.csv, missp.csv, misspsp.csv
### OUTPUTS:
## augmented files for plotting and tables, named as follows:
## original simulation statistics plus test set MSE
# compaggplot.csv, partaggplot.csv, compdistplot.csv,
# spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, 
# exspplot.csv, misspplot.csv, misspspplot.csv
## bias - variance results
# compaggcv.csv, partaggcv.csv, compdistcv.csv,
# spcompaggcv.csv, sppartaggcv.csv, spcompdistcv.csv, 
# exspcv.csv, misspcv.csv, misspspcv.csv
####################################################################

########################################
# 9 scenarios:
# 1A. completely aggregated
# 2B. partially aggregated
# 3C. completely distributed
# 4A. sparse group completely aggregated
# 5B. sparse group parially aggregated
# 6C. sparse group completely distributed
# 7B. partially aggregated extra sparse
# 8B. partially aggregated misspecified
# 9B. partially aggregated misspecified sparse
########################################

########################################
# LOAD DATASETS
########################################
# set the working and data directories
# setwd()
# datadir <- paste0(getwd(), '/MainSimulationResults/')
compagg <- read.csv(paste0(datadir, 'compagg.csv'))
partagg <- read.csv(paste0(datadir, 'partagg.csv'))
compdist <- read.csv(paste0(datadir, 'compdist.csv'))
spcompagg <- read.csv(paste0(datadir, 'spcompagg.csv'))
sppartagg <- read.csv(paste0(datadir, 'sppartagg.csv'))
spcompdist <- read.csv(paste0(datadir, 'spcompdist.csv'))
exsp <- read.csv(paste0(datadir, 'exsp.csv'))
missp <- read.csv(paste0(datadir, 'missp.csv'))
misspsp <- read.csv(paste0(datadir, 'misspsp.csv'))

# remove the first column (row names)
compagg <- compagg[,2:409]
partagg <- partagg[,2:409] 
compdist <- compdist[,2:409]
spcompagg <- spcompagg[,2:409]
sppartagg <- sppartagg[,2:409]
spcompdist <- spcompdist[,2:409]
exsp <- exsp[,2:409]
missp <- missp[,2:409]
misspsp <- misspsp[,2:409]

####################################################################
# CACLULATE ERROR ON TEST DATA
####################################################################
########################################
# SET SOME SIMULATION PARAMETERS
########################################
# number of subjects
n <- 50
# 20*20 grid = 400 pixels
# dim 1 is number of rows of image
dim1 <- 20
# dim 2 is number of columns of image
dim2 <- 20
p <- dim1*dim2

######################################################
# DEFINE GROUP STRUCTURES
# 16 groups of 25 pixels
######################################################
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
GroupsA <- c(rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5),
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 4,
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 8,
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 12)
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
numbers <- cbind(1:16, c(4:16, 1:3), c(7:16, 1:6), c(10:16, 1:9))
numbers.list <- apply(numbers, 1, list)
# function to create subblock (5*5 square)
subblock <- function(x){
  x <- unlist(x)
  data <- c(rep(c(rep(x[1], 3), rep(x[2], 2)), 2), rep(x[1], 3), rep(x[4], 2), rep(c(rep(x[2], 2), x[4], rep(x[3], 2)), 2))
  return(matrix(data, nrow=5, byrow=TRUE))
}
# create subblocks
subblocks <- lapply(numbers.list, subblock)
# combine subblocks in rows
row1 <- cbind(subblocks[[1]], subblocks[[2]], subblocks[[3]], subblocks[[4]])
row2 <- cbind(subblocks[[5]], subblocks[[6]], subblocks[[7]], subblocks[[8]])
row3 <- cbind(subblocks[[9]], subblocks[[10]], subblocks[[11]], subblocks[[12]])
row4 <- cbind(subblocks[[13]], subblocks[[14]], subblocks[[15]], subblocks[[16]])
GroupsB <- as.vector(rbind(row1, row2, row3, row4))
# remove stuff we don't need
rm('row1', 'row2', 'row3', 'row4', 'numbers', 'numbers.list', 'subblocks')
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
GroupsC <- c(rep(1:16, 25))
########################################
# OPTIONAL STEP: Visualize the groups
########################################
# library(fields)
# image.plot(matrix(GroupsA, nrow=20, byrow=TRUE), main='Group Structure')
# image.plot(matrix(GroupsB, nrow=20, byrow=TRUE), main='Group Structure')
# image.plot(matrix(GroupsC, nrow=20, byrow=TRUE), main='Group Structure')

######################################################
# DEFINE TRUE COEFFICIENTS
######################################################
# set coefficient magnitude
beta.value <- 3  
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
### TRUE COEFFICIENTS 1A: Complete Group
trueBetas1A <- as.numeric(GroupsA == 7)*beta.value
### TRUE COEFFICIENTS 4A: Sparse Group 
trueBetas4A <- as.numeric(GroupsA == 7)*beta.value
trueBetas4A[c(113, 114, 115, 133, 134, 135, 154, 155, 175, 195)] <- 0
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
### TRUE COEFFICIENTS 2B: Complete Group
trueBetas2B <- as.numeric(GroupsB == 10)*beta.value
### TRUE COEFFICIENTS 5B: Sparse Group
trueBetas5B <- as.numeric(GroupsB == 10)*beta.value
trueBetas5B[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 7B: Extra Sparse Group 
trueBetas7B <- as.numeric(GroupsB == 10)*beta.value
trueBetas7B[c(44, 45, 111, 112, 113, 131, 132, 133, 151, 152, 153, 209, 210, 229, 230, 267, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 8B: Misspecified Group 
trueBetas8B <- as.numeric(GroupsB == 10)*beta.value
rotate <- function(x) t(apply(x, 2, rev))
trueBetas8B <- as.vector(rotate(matrix(trueBetas8B, nrow=20)))
### TRUE COEFFICIENTS 9B: Misspecified Sparse Group 
trueBetas9B <- as.numeric(GroupsB == 10)*beta.value
trueBetas9B[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
rotate <- function(x) t(apply(x, 2, rev))
trueBetas9B <- as.vector(rotate(matrix(trueBetas9B, nrow=20)))
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
### TRUE COEFFICIENTS 3C: Complete Group
trueBetas3C <- as.numeric(GroupsC == 7)*beta.value
### TRUE COEFFICIENTS 6C: Sparse Group
trueBetas6C <- as.numeric(GroupsC == 7)*beta.value
set.seed(1)
trueBetas6C[which(trueBetas6C == 3)[sample(1:25, 10)]] <- 0
########################################
# OPTIONAL STEP: Visualize the true coefficients
########################################
# image.plot(matrix(trueBetas1A, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas2B, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas3C, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas4A, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas5B, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas6C, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas7B, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas8B, nrow=20, byrow=TRUE), main='True Coefficients')
# image.plot(matrix(trueBetas9B, nrow=20, byrow=TRUE), main='True Coefficients')

######################################################
# CACLULATE ERROR ON TEST DATA
######################################################
########################################
# function to calculate prediction error 
# on test set for each seed
########################################
prederr <- function(simdata, trueBetas){
  n <- 50
  p <- 400
  simdata$mse.y.test <- rep(NA, 2100)
  for (seed in 101:200){
    set.seed(seed)
    X <- scale(matrix(rnorm(n*p), nrow=n, ncol=p), center=TRUE, scale=FALSE)
    epsilon <- rnorm(n, mean=0, sd=2)
    y.test <- X %*% trueBetas + epsilon
    # get the subset of 21 rows corresponding to a given seed
    betahat.set <- t(simdata[simdata$seed==(seed - 100), 9:408])
    yhat.test <- X %*% betahat.set
    # calculate msey.test
    simdata$mse.y.test[simdata$seed==(seed - 100)] <- colMeans(sweep(yhat.test, 1, y.test)^2)
  }
  return(simdata)
}

########################################
# CALCULATE TEST MEAN SQUARED ERROR
# for each set of true coefficients
########################################
compagg2 <- prederr(compagg, trueBetas1A)
spcompagg2 <- prederr(spcompagg, trueBetas4A)
partagg2 <- prederr(partagg, trueBetas2B)
sppartagg2 <- prederr(sppartagg, trueBetas5B)
compdist2 <- prederr(compdist, trueBetas3C)
spcompdist2 <- prederr(spcompdist, trueBetas6C)
exsp2 <- prederr(exsp, trueBetas7B)
missp2 <- prederr(missp, trueBetas8B)
misspsp2 <- prederr(misspsp, trueBetas9B)

########################################
# save stats output data for making tables and plots
########################################
write.csv(compagg2[,c(1:8,409)], 'compaggplot.csv', row.names=FALSE)
write.csv(spcompagg2[,c(1:8,409)], 'spcompaggplot.csv', row.names=FALSE)
write.csv(partagg2[,c(1:8,409)], 'partaggplot.csv', row.names=FALSE)
write.csv(sppartagg2[,c(1:8,409)], 'sppartaggplot.csv', row.names=FALSE)
write.csv(compdist2[,c(1:8,409)], 'compdistplot.csv', row.names=FALSE)
write.csv(spcompdist2[,c(1:8,409)], 'spcompdistplot.csv', row.names=FALSE)
write.csv(exsp2[,c(1:8,409)], 'exspplot.csv', row.names=FALSE)
write.csv(missp2[,c(1:8,409)], 'misspplot.csv', row.names=FALSE)
write.csv(misspsp2[,c(1:8,409)], 'misspspplot.csv', row.names=FALSE)

####################################################################
# ESTIMATE BIAS-VARIANCE DECOMPOSITION OF MSE
# for 100 new observations
####################################################################
########################################
# function to estimate bias-variance
# decomposition of mse
# for 100 new observations
########################################
biasvar <- function(simdata, trueBetas){
  set.seed(201)
  n <- 100
  p <- 400
  X <- scale(matrix(rnorm(n*p), nrow=n, ncol=p), center=TRUE, scale=FALSE)
  testdata <- data.frame(id=1:n)
  testdata$X.beta <- X %*% trueBetas
  # get alpha / gamma values and estimated betas
  betahat <- simdata[,c(2:3,9:408)]
  # sort betahat by the alpha / gamma values
  betahat.sorted <- betahat[order(betahat$alpha, betahat$gamma),]
  yhat.test <- X %*% t(betahat.sorted[,3:402])
  # for each alpha / gamma combination
  # and each observation
  # calculate variance of the predicted value
  # also calculate estimated bias and bias squared
  alphagamma <- betahat[1:21,1:2]
  for (i in 1:21){
    predvals <- yhat.test[,(100*i - 99):(100*i)]
    testdata[[paste0('yhat.var.alpha', alphagamma[i,1] ,'gamma', alphagamma[i,2])]] <- apply(predvals, 1, var)
    testdata[[paste0('yhat.bias.alpha', alphagamma[i,1] ,'gamma', alphagamma[i,2])]] <- (rowMeans(predvals) - testdata$X.beta)
    testdata[[paste0('yhat.bias2.alpha', alphagamma[i,1] ,'gamma', alphagamma[i,2])]] <- (rowMeans(predvals) - testdata$X.beta)^2
    testdata[[paste0('yhat.mse.alpha', alphagamma[i,1] ,'gamma', alphagamma[i,2])]] <- apply(cbind(testdata$X.beta, predvals), 1, function(x) mean((x[2:101] - (rnorm(100, mean=0, sd=2) + x[1]))^2))
  }
  return(testdata)
}

compaggbv <- biasvar(compagg, trueBetas1A)
spcompaggbv <- biasvar(spcompagg, trueBetas4A)
partaggbv <- biasvar(partagg, trueBetas2B)
sppartaggbv <- biasvar(sppartagg, trueBetas5B)
compdistbv <- biasvar(compdist, trueBetas3C)
spcompdistbv <- biasvar(spcompdist, trueBetas6C)
exspbv <- biasvar(exsp, trueBetas7B)
misspbv <- biasvar(missp, trueBetas8B)
misspspbv <- biasvar(misspsp, trueBetas9B)

########################################
# save output data for making tables and plots
########################################
write.csv(compaggbv, 'compaggbv.csv', row.names=FALSE)
write.csv(spcompaggbv, 'spcompaggbv.csv', row.names=FALSE)
write.csv(partaggbv, 'partaggbv.csv', row.names=FALSE)
write.csv(sppartaggbv, 'sppartaggbv.csv', row.names=FALSE)
write.csv(compdistbv, 'compdistbv.csv', row.names=FALSE)
write.csv(spcompdistbv, 'spcompdistbv.csv', row.names=FALSE)
write.csv(exspbv, 'exspbv.csv', row.names=FALSE)
write.csv(misspbv, 'misspbv.csv', row.names=FALSE)
write.csv(misspspbv, 'misspspbv.csv', row.names=FALSE)
