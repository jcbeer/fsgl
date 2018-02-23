####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Save diagnosis info for test subjects for plotting
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# phenosrsall.csv
# srs_test.txt
### OUTPUTS:
# phenotest.txt
####################################################################

######################################################
# LOAD DATA
######################################################
setwd('./data')
abideall <- read.csv('phenosrsall.csv')
test <- read.table('srs_test.txt', header=TRUE)

######################################################
# REFORMAT AND SAVE PHENOTEST DATA
######################################################
phenotest <- merge(test, abideall, by='SUB_ID', all.x=TRUE, all.y=FALSE)
phenotest2 <- phenotest[,c('SUB_ID', "SRS_ADJUSTED", "DX_GROUP", "AGE_AT_SCAN")]
# DX_GROUP 1 = autism, 2 = typically developing
write.table(phenotest2, 'phenotest.txt', col.names = FALSE, row.names = FALSE, quote=FALSE)