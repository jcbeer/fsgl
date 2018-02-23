####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Randomly assign subjects to training and test set
# total N = 219
# training set N = round(0.8*219, 0) = 175
# test set N = round(0.2*219, 0) = 44
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# fileids.csv
### OUTPUTS:
# fileids_training.txt
# fileids_test.txt
####################################################################

######################################################
# LOAD DATA
######################################################
setwd('./data')
fileids <- read.csv('fileids.csv', stringsAsFactors=FALSE)
# reformat file IDs to match directory naming: "s"+ < last 5 digits of ID >
fileids_reformat <- sapply(fileids[,1], function(x) paste0('s', substr(x, nchar(x)-5+1, nchar(x))))

######################################################
# RANDOMLY ASSIGN TO TRAINING AND TEST SETS
######################################################
set.seed(13)
fileids_training <- fileids_reformat[sample(1:219, 175, replace=FALSE)]
fileids_test <- fileids_reformat[!(fileids_reformat %in% fileids_training)]

# reorder and unname
fileids_training <- unname(fileids_training[order(fileids_training)])
fileids_test <- unname(fileids_test[order(fileids_test)])

# save files of IDs
write.table(fileids_training, 'fileids_training.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(fileids_test, 'fileids_test.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
