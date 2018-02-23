####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Calculate descriptive statistics (Table S11)
# Adjust the SRS scores using linear regression (Table S12)
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# phenosrsall.csv
# fileids_training.txt
# fileids_test.txt
### OUTPUTS:
# srs_train.txt: original and adjusted SRS for train data
# srs_test.txt: original and adjusted SRS for test data
####################################################################
# Calculate descriptive statistics and 
# adjust Social Responsiveness Scale scores (SRS_RAW_TOTAL) 
# for
# age (AGE_AT_SCAN)
# full-scale IQ (FIQ)
# site of acquisition (SITE_ID)
# eye status at scan (EYE_STATUS_AT_SCAN, 1=open, 2=closed)
# mean framewise displacement (func_mean_fd)
# using linear regression model fit to the training data
# will also tabulate diagnosis group (DX_GROUP, 1=Autism, 2=Control)
####################################################################

######################################################
# LOAD DATA
######################################################
setwd('./data')
pheno <- read.csv('phenosrsall.csv')
trainids <- read.table('fileids_training.txt')
testids <- read.table('fileids_test.txt')

# remove the 's' from each ID
trainids_num <- sapply(trainids[,1], function(x) as.numeric(strsplit(as.character(x), 's')[[1]][2]))
testids_num <- sapply(testids[,1], function(x) as.numeric(strsplit(as.character(x), 's')[[1]][2]))

# divide phenotype dataset into train and test
phenotrain <- pheno[pheno$SUB_ID %in% trainids_num,]
phenotest <- pheno[pheno$SUB_ID %in% testids_num,]

######################################################
# CALCULATE DESCRIPTIVE STATISTICS (Table S11)
# overall
# training set
# test set
######################################################
meansd <- function(datasets, variable, digits){
  overall <- datasets[[1]][,variable]
  train <- datasets[[2]][,variable]
  test <- datasets[[3]][,variable]
  overall.mean <- round(mean(overall, na.rm=TRUE), digits)
  overall.sd <- round(sd(overall, na.rm=TRUE), digits)
  train.mean <- round(mean(train, na.rm=TRUE), digits)
  train.sd <- round(sd(train, na.rm=TRUE), digits)
  test.mean <- round(mean(test, na.rm=TRUE), digits)
  test.sd <- round(sd(test, na.rm=TRUE), digits)
  return(data.frame(overall.mean, overall.sd, train.mean, train.sd, test.mean, test.sd))
}

npercent <- function(datasets, variable, digits){
  overall <- datasets[[1]][,variable]
  train <- datasets[[2]][,variable]
  test <- datasets[[3]][,variable]
  overall.n <- table(overall, useNA='always')
  overall.pct <- round(table(overall, useNA='always')/sum(table(overall, useNA='always'))*100, digits)
  train.n <- table(train, useNA='always')
  train.pct <- round(table(train, useNA='always')/sum(table(train, useNA='always'))*100, digits)
  test.n <- table(test, useNA='always')
  test.pct <- round(table(test, useNA='always')/sum(table(test, useNA='always'))*100, digits)
  return(list(overall.n, overall.pct, train.n, train.pct, test.n, test.pct))
}

meansd(datasets=list(pheno, phenotrain, phenotest), variable='AGE_AT_SCAN', digits=1)
meansd(datasets=list(pheno, phenotrain, phenotest), variable='FIQ', digits=1)
meansd(datasets=list(pheno, phenotrain, phenotest), variable='func_mean_fd', digits=2)
meansd(datasets=list(pheno, phenotrain, phenotest), variable='SRS_RAW_TOTAL', digits=1)

npercent(datasets=list(pheno, phenotrain, phenotest), variable='DX_GROUP', digits=1)
npercent(datasets=list(pheno, phenotrain, phenotest), variable='SITE_ID', digits=1)
npercent(datasets=list(pheno, phenotrain, phenotest), variable='EYE_STATUS_AT_SCAN', digits=1)

######################################################
# test for significant differences between train and test data
######################################################
t.test(phenotrain$AGE_AT_SCAN, phenotest$AGE_AT_SCAN)
t.test(phenotrain$FIQ, phenotest$FIQ)
t.test(phenotrain$func_mean_fd, phenotest$func_mean_fd)
t.test(phenotrain$SRS_RAW_TOTAL, phenotest$SRS_RAW_TOTAL)

fisher.test(table(c(phenotrain$DX_GROUP, phenotest$DX_GROUP), c(rep('train', 175), rep('test', 44))))
fisher.test(table(c(phenotrain$SITE_ID, phenotest$SITE_ID), c(rep('train', 175), rep('test', 44))))
fisher.test(table(c(phenotrain$EYE_STATUS_AT_SCAN, phenotest$EYE_STATUS_AT_SCAN), c(rep('train', 175), rep('test', 44))))

######################################################
# fit linear regression to impute a missing framewise 
# displacement in the training set
######################################################
fd_fit <- lm(func_mean_fd ~ AGE_AT_SCAN + FIQ + SRS_RAW_TOTAL + as.factor(SITE_ID) + as.factor(EYE_STATUS_AT_SCAN), data=phenotrain)
missing_data_subject <- phenotrain[is.na(phenotrain$func_mean_fd),c('AGE_AT_SCAN', 'FIQ', 'SRS_RAW_TOTAL', 'SITE_ID', 'EYE_STATUS_AT_SCAN')]
missing_data_subject_predicted_fd <- predict(fd_fit, newdata=missing_data_subject)
# replace this value in  the training data
phenotrain$func_mean_fd_imputed <- phenotrain$func_mean_fd
phenotrain$func_mean_fd_imputed[is.na(phenotrain$func_mean_fd_imputed)] <- missing_data_subject_predicted_fd

######################################################
# FIT LINEAR REGRESSION TO ADJUST SRS USING TRAINING DATA
# (Table S12)
######################################################
# fit the model
trainfit <- lm(SRS_RAW_TOTAL ~ AGE_AT_SCAN + FIQ + func_mean_fd_imputed + as.factor(SITE_ID) + as.factor(EYE_STATUS_AT_SCAN), data=phenotrain)
summary(trainfit)
# plot(phenotrain$SRS_RAW_TOTAL, trainfit$fitted.values)
# plot(phenotrain$SRS_RAW_TOTAL, trainfit$residuals)

# predict SRS for the test set and calculate residuals
phenotest$func_mean_fd_imputed <- phenotest$func_mean_fd
test.residuals <- phenotest$SRS_RAW_TOTAL - predict(trainfit, newdata=phenotest)
plot(phenotest$SRS_RAW_TOTAL, test.residuals)

######################################################
# save SRS raw and adjusted 
######################################################
srs_train <- data.frame(SUB_ID=phenotrain$SUB_ID, SRS_RAW_TOTAL=phenotrain$SRS_RAW_TOTAL, SRS_ADJUSTED=trainfit$residuals)
srs_test <- data.frame(SUB_ID=phenotest$SUB_ID, SRS_RAW_TOTAL=phenotest$SRS_RAW_TOTAL, SRS_ADJUSTED=test.residuals)
write.table(srs_train, 'srs_train.txt', row.names=FALSE)
write.table(srs_test, 'srs_test.txt', row.names=FALSE)
