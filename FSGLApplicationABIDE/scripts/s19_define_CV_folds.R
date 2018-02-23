####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Make 10 sets of 5-fold cross validation assignments
# stratified by adjusted SRS
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# srs_train.txt
### OUTPUTS:
# folds10.csv
####################################################################

######################################################
# LOAD DATA
######################################################
setwd('./data')
Ytraindata <- read.table('srs_train.txt', header=TRUE)
Ytrain <- Ytraindata[,3]

######################################################
# MAKE CV FOLD ASSIGNMENTS
######################################################
# get index for ordering of adjusted SRS
index <- sort.int(Ytrain, index.return=TRUE, decreasing=TRUE)$ix

# create 10 random CV fold assignments
# the numbers 1:5 sorted randomly, repeated 35 times
folds_unsorted <- matrix(nrow=175, ncol=10)
set.seed(13)
for(i in 1:10){
  folds_unsorted[,i] <- as.vector(replicate(35, sample(1:5)))
}

# reorder the rows
folds <- matrix(nrow=175, ncol=10)
for(j in 1:175){
  folds[index[j],] <- folds_unsorted[j,]
}

# check distributions
boxplot(Ytrain ~ folds[,1])
boxplot(Ytrain ~ folds[,2])
boxplot(Ytrain ~ folds[,3])
boxplot(Ytrain ~ folds[,4])
boxplot(Ytrain ~ folds[,5])
boxplot(Ytrain ~ folds[,6])
boxplot(Ytrain ~ folds[,7])
boxplot(Ytrain ~ folds[,8])
boxplot(Ytrain ~ folds[,9])
boxplot(Ytrain ~ folds[,10])

# save as csv
write.table(folds, 'folds10.csv', col.names=FALSE, row.names=FALSE, sep=',')
