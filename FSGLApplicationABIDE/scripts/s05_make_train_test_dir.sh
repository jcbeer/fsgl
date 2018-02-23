#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Move training and test data into different parent directories
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# fileids_training.txt
# fileids_test.txt
### OUTPUTS:
# moves subject fMRI data into train and test parent directories 
####################################################################

#######################################################
# NAVIGATE TO THE PROPER DIRECTORY
#######################################################
cd ABIDE/data/subjectdata

#######################################################
# MAKE NEW DIRECTORIES FOR TRAIN AND TEST DATA
#######################################################
mkdir ABIDE/data/subjectdata/train
mkdir ABIDE/data/subjectdata/test

#######################################################
# READ IN SUBJECT IDs FOR TRAIN AND TEST DATA
#######################################################
set train=`cat /ABIDE/data/fileids_training.txt`
set test=`cat /ABIDE/data/fileids_test.txt`

#######################################################
# MOVE EACH SUBJECT SUBDIRECTORY TO THE PROPER PLACE
#######################################################
foreach subj ($train)
mv $subj train
end

foreach subj ($test)
mv $subj test
end