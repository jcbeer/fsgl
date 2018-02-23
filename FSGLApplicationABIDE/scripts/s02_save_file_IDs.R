####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Save phenotype data and subject file IDs 
# for the subset with nonmissing SRS scores
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# phenomat.csv
# phenolab.csv
# Phenotypic_V1_0b_preprocessed1.csv
# downloaded from http://preprocessed-connectomes-project.org/abide/download.html
### OUTPUTS:
# phenosrs.csv: phenotype data for 219 subjects with nonmissing SRS score
# phenosrsall.csv: more phenotype data for 219 subjects with nonmissing SRS score
# fileids.csv: file IDs for 219 subjects with nonmissing SRS score
####################################################################

######################################################
# LOAD DATA
######################################################
# set working directory
setwd('./data/')
pheno <- read.csv('phenomat.csv', header=FALSE)
phenolab <- read.csv('phenolab.csv')
colnames(pheno) <- names(phenolab)
rm(list='phenolab')
# load entire ABIDE phenotype data
abide <- read.csv('Phenotypic_V1_0b_preprocessed1.csv')

######################################################
# get phenotype data for only those who have 
# nonmissing Social Responsiveness Scale (SRS) score
######################################################
srs <- pheno[!is.na(pheno$SRS_RAW_TOTAL),]
# save the srs data
write.csv(srs, 'phenosrs.csv', row.names=FALSE)
# get complete phenotype data for srs subset
srsids <- srs$SUB_ID
abidesrs <- abide[abide$SUB_ID %in% srsids,]
# save this data
write.csv(abidesrs, 'phenosrsall.csv', row.names=FALSE)

######################################################
# create the filenames
######################################################
parts <- cbind(as.character(abidesrs$SITE_ID), abidesrs$SUB_ID)
fileids <- apply(parts, 1, function(x) paste0(x[1], '_00', x[2]))
write.csv(filenames, 'fileids.csv', row.names=FALSE)
