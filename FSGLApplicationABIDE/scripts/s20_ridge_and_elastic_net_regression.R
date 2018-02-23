####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Do ridge and elastic net regression on ABIDE training data 
# Save ridge coefficients for making adaptive penalty weights
# Calculate Sum of Squares Total (SST) for training and test sets
# Results reported in Table 3
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# folds10.csv
# srs_train.txt
# trainXstd_5476.txt
# srs_test.txt
# testXstd_5476.txt
### OUTPUTS:
# betaridge.csv
# results reported in Table 3
####################################################################

######################################################
# LOAD DATA
######################################################
library(glmnet)
setwd('./data')
folds <- read.csv('folds10.csv', header=FALSE)
# TRAINING SET
Ytraindata <- read.table('srs_train.txt', header=TRUE)
Ytrain <- Ytraindata[,3]
rm('Ytraindata')
Xtrain <- read.table('trainXstd_5476.txt', sep=',')
# TEST SET
Ytestdata <- read.table('srs_test.txt', header=TRUE)
Ytest <- Ytestdata[,3]
rm('Ytestdata')
Xtest <- read.table('testXstd_5476.txt', sep=',')

######################################################
# ELASTIC NET
######################################################
# do cross-validation over grid of alpha and lambda
# alpha = 0 : ridge regression
# alpha = 1 : lasso
alpha.grid <- seq(0, 1, length=100) 
# set lambda grid
exponent <- seq(-2, 9.5, length=100)
lambda.grid <- rev(exp(exponent))
# do cross-validations for each alpha
# NOTE: will do 5 rounds of CV, using 5 different 
# CV Fold assignments and average the cross-validation error
elasticcv1 <- lapply(alpha.grid, function(curAlpha){
  cv.glmnet(as.matrix(Xtrain), Ytrain, alpha=curAlpha, lambda=lambda.grid, foldid=folds[,1])
})
elasticcv2 <- lapply(alpha.grid, function(curAlpha){
  cv.glmnet(as.matrix(Xtrain), Ytrain, alpha=curAlpha, lambda=lambda.grid, foldid=folds[,2])
})
elasticcv3 <- lapply(alpha.grid, function(curAlpha){
  cv.glmnet(as.matrix(Xtrain), Ytrain, alpha=curAlpha, lambda=lambda.grid, foldid=folds[,3])
})
elasticcv4 <- lapply(alpha.grid, function(curAlpha){
  cv.glmnet(as.matrix(Xtrain), Ytrain, alpha=curAlpha, lambda=lambda.grid, foldid=folds[,4])
})
elasticcv5 <- lapply(alpha.grid, function(curAlpha){
  cv.glmnet(as.matrix(Xtrain), Ytrain, alpha=curAlpha, lambda=lambda.grid, foldid=folds[,5])
})
# collect the optimum lambda for each alpha
# lam: value of lambda that gives minimum mean cross validated error
# alph: current alpha value
# cvup: upper curve
# cvm: mean cross validation error
optimumPerAlpha1 <- sapply(seq_along(alpha.grid), function(curi){
  curcvs <- elasticcv1[[curi]]
  curAlpha <- alpha.grid[curi]
  indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
  c(lam=curcvs$lambda.min, alph=curAlpha, cvup=curcvs$cvup[indOfMin], cvm=curcvs$cvm[indOfMin])
})
optimumPerAlpha2 <- sapply(seq_along(alpha.grid), function(curi){
  curcvs <- elasticcv2[[curi]]
  curAlpha <- alpha.grid[curi]
  indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
  c(lam=curcvs$lambda.min, alph=curAlpha, cvup=curcvs$cvup[indOfMin], cvm=curcvs$cvm[indOfMin])
})
optimumPerAlpha3 <- sapply(seq_along(alpha.grid), function(curi){
  curcvs <- elasticcv3[[curi]]
  curAlpha <- alpha.grid[curi]
  indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
  c(lam=curcvs$lambda.min, alph=curAlpha, cvup=curcvs$cvup[indOfMin], cvm=curcvs$cvm[indOfMin])
})
optimumPerAlpha4 <- sapply(seq_along(alpha.grid), function(curi){
  curcvs <- elasticcv4[[curi]]
  curAlpha <- alpha.grid[curi]
  indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
  c(lam=curcvs$lambda.min, alph=curAlpha, cvup=curcvs$cvup[indOfMin], cvm=curcvs$cvm[indOfMin])
})
optimumPerAlpha5 <- sapply(seq_along(alpha.grid), function(curi){
  curcvs <- elasticcv5[[curi]]
  curAlpha <- alpha.grid[curi]
  indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
  c(lam=curcvs$lambda.min, alph=curAlpha, cvup=curcvs$cvup[indOfMin], cvm=curcvs$cvm[indOfMin])
})
# collect all cvm values and take mean
cvm1 <- matrix(nrow=100, ncol=100)
for(i in 1:100){
  cvm1[i,] <- elasticcv1[[i]]$cvm
}
cvm2 <- matrix(nrow=100, ncol=100)
for(i in 1:100){
  cvm2[i,] <- elasticcv2[[i]]$cvm
}
cvm3 <- matrix(nrow=100, ncol=100)
for(i in 1:100){
  cvm3[i,] <- elasticcv3[[i]]$cvm
}
cvm4 <- matrix(nrow=100, ncol=100)
for(i in 1:100){
  cvm4[i,] <- elasticcv4[[i]]$cvm
}
cvm5 <- matrix(nrow=100, ncol=100)
for(i in 1:100){
  cvm5[i,] <- elasticcv5[[i]]$cvm
}
# first row is alpha = 0 (ridge), last row alpha = 1 (lasso)
# first column is largest lambda, last column is smallest
cvm_mean <- (cvm1 + cvm2 + cvm3 + cvm4 + cvm5)/5
# heat map of cross validation error
image(t(cvm_mean), xlab='Lambda', ylab='Alpha') # plot is flipped around
# calculate minimum for each row
alpha_mincvm <- apply(cvm_mean, 1, min)
# find index of minimum for each row (which lambda for each alpha)
lambda_index <- apply(cvm_mean, 1, function(x) which(x==min(x)))
# find minimum row (which alpha)
min(alpha_mincvm)
which.min(alpha_mincvm)
# alpha = 1 (ridge regression) is minimum
# will use next largest alpha for elastic net
alpha_index <- which(alpha_mincvm==min(alpha_mincvm[2:100]))
alpha_index
lambda_index[alpha_index]
# row 2, column 34
cvm_mean[2,34]
alpha.grid[2]
lambda.grid[34]
# plot this point
cvm_mean[2,34] <- 1500
image(t(-cvm_mean), xlab='Lambda', ylab='Alpha', col=terrain.colors(500))
# save optimal values
optAlpha <- alpha.grid[2]
optLambda <- lambda.grid[34]
# # make plots
# par(mfcol=c(5, 2), mar=c(2,2,2,2))
# # plot optimal lambda for each alpha
# plot(alpha.grid, log(optimumPerAlpha1[1,]), type='l', ylab='Optimum Lambda', xlab='Alpha', main='Alpha vs log(Lambda)')
# plot(alpha.grid, log(optimumPerAlpha2[1,]), type='l', ylab='Optimum Lambda', xlab='Alpha', main='Alpha vs log(Lambda)')
# plot(alpha.grid, log(optimumPerAlpha3[1,]), type='l', ylab='Optimum Lambda', xlab='Alpha', main='Alpha vs log(Lambda)')
# plot(alpha.grid, log(optimumPerAlpha4[1,]), type='l', ylab='Optimum Lambda', xlab='Alpha', main='Alpha vs log(Lambda)')
# plot(alpha.grid, log(optimumPerAlpha5[1,]), type='l', ylab='Optimum Lambda', xlab='Alpha', main='Alpha vs log(Lambda)')
# # plot alpha versus the mean cv error for optimum lambda 
# plot(alpha.grid, optimumPerAlpha1[4,], type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs CVE', ylim=c(1500, 1750))
# plot(alpha.grid, optimumPerAlpha2[4,], type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs CVE', ylim=c(1500, 1750))
# plot(alpha.grid, optimumPerAlpha3[4,], type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs CVE', ylim=c(1500, 1750))
# plot(alpha.grid, optimumPerAlpha4[4,], type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs CVE', ylim=c(1500, 1750))
# plot(alpha.grid, optimumPerAlpha5[4,], type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs CVE', ylim=c(1500, 1750))
# # calculate the mean CV error across alpha
# meancve_alpha <- colMeans(rbind(optimumPerAlpha1[4,], optimumPerAlpha2[4,], optimumPerAlpha3[4,], optimumPerAlpha4[4,], optimumPerAlpha5[4,]))
# # plot
# dev.off()
# plot(alpha.grid, meancve_alpha, type='l', ylab='Mean CV error', xlab='Alpha', main='Alpha vs Mean CVE', ylim=c(1600, 1700))
# points(alpha.grid[which(meancve_alpha==min(meancve_alpha[2:100]))], min(meancve_alpha[2:100]))
# fit to the training data
elastic.mod <- glmnet(as.matrix(Xtrain), Ytrain, alpha=optAlpha, lambda=lambda.grid)
elastic.pred.train <- predict(elastic.mod, s=optLambda, newx=as.matrix(Xtrain))
elastic.pred.test <- predict(elastic.mod, s=optLambda, newx=as.matrix(Xtest))
# estimated test error
mean((elastic.pred.train - Ytrain)^2)
cor.test(elastic.pred.train, Ytrain)
mean((elastic.pred.test - Ytest)^2)
cor.test(elastic.pred.test, Ytest)

######################################################
# RIDGE REGRESSION
######################################################
# do cross-validation over grid of lambda
# alpha = 0 : ridge regression
# set lambda grid
exponent <- seq(-2, 9.5, length=100)
lambda.grid <- rev(exp(exponent))
# NOTE: will do 5 rounds of CV, using 5 different 
# CV Fold assignments and average the cross-validation error
ridgecv1 <- cv.glmnet(as.matrix(Xtrain), Ytrain, foldid=folds[,1], lambda=lambda.grid, alpha=0)
ridgecv2 <- cv.glmnet(as.matrix(Xtrain), Ytrain, foldid=folds[,2], lambda=lambda.grid, alpha=0)
ridgecv3 <- cv.glmnet(as.matrix(Xtrain), Ytrain, foldid=folds[,3], lambda=lambda.grid, alpha=0)
ridgecv4 <- cv.glmnet(as.matrix(Xtrain), Ytrain, foldid=folds[,4], lambda=lambda.grid, alpha=0)
ridgecv5 <- cv.glmnet(as.matrix(Xtrain), Ytrain, foldid=folds[,5], lambda=lambda.grid, alpha=0)
# calculate the mean CV error
meancve <- colMeans(rbind(ridgecv1$cvm, ridgecv2$cvm, ridgecv3$cvm, ridgecv4$cvm, ridgecv5$cvm))
min(meancve)
which(meancve==min(meancve))
overall.lambda.min <- lambda.grid[which(meancve==min(meancve))]
# plot results
# par(mar=c(5, 5, 6, 2), mfrow=c(2,3))
# plot(ridgecv1, main='Ridge penalty (alpha = 0)')
# plot(ridgecv2, main='Ridge penalty (alpha = 0)')
# plot(ridgecv3, main='Ridge penalty (alpha = 0)')
# plot(ridgecv4, main='Ridge penalty (alpha = 0)')
# plot(ridgecv5, main='Ridge penalty (alpha = 0)')
# plot(log(lambda.grid), meancve, main='Ridge penalty mean CVE', type='l')
# abline(v=log(lambda.grid[which(meancve==min(meancve))]), lty=2)
# fit to the training data
ridge.mod <- glmnet(as.matrix(Xtrain), Ytrain, alpha=0, lambda=lambda.grid)
ridge.pred.train <- predict(ridge.mod, s=overall.lambda.min, newx=as.matrix(Xtrain))
ridge.pred.test <- predict(ridge.mod, s=overall.lambda.min, newx=as.matrix(Xtest))
# estimated test error
mean((ridge.pred.train - Ytrain)^2)
cor.test(ridge.pred.train, Ytrain)
mean((ridge.pred.test - Ytest)^2)
cor.test(ridge.pred.test, Ytest)
# save ridge coefficients to csv
# NOTE: for FSGL adaptive weights, 
# lambda = 2403 (approx) was used
# to calculate the ridge coefficient estimates
# ridge.mod2403 <- glmnet(as.matrix(Xtrain), Ytrain, alpha=0, lambda=2403)
# beta.ridge <- ridge.mod2403$beta
# write.table(beta.ridge, 'betaridge.csv', sep=',', row.names = FALSE, col.names = FALSE)

#################################################
# Calculate SST for training and test sets
#################################################
sum((Ytrain - mean(Ytrain))^2)/175
sum((Ytest - mean(Ytest))^2)/44