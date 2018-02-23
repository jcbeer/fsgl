#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Make MNI152 brain mask 
# using 3*3*3 mm^3 voxel template from subject resting state dataset
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# MNI152_T1_2009c+tlrc, from AFNI files
# s50433_func_preproc.nii.gz, one subject’s fMRI data
### OUTPUTS:
# MNI152mask_3mm+tlrc: 3mm voxel brain mask, AFNI format
# MNI152maskbin_3mm+tlrc: 3mm voxel binary brain mask, AFNI format
# MNI152maskbin_3mm.nii: 3mm voxel binary brain mask, NIFTI format
####################################################################

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

# move to data directory
cd '/ABIDE/data'

# copy a RS dataset to the current directory 
# to use as template for resampling
cp subjectdata/train/s50433/s50433_func_preproc.nii.gz .

######################################################
# CREATE BRAIN MASK
######################################################
# resample MNI152 anatomical to 3mm voxels
3dresample -master s50433_func_preproc.nii.gz -prefix MNI152mask_3mm -input MNI152_T1_2009c+tlrc

# create a 3mm binary mask
3dcalc -a 'MNI152mask_3mm+tlrc' -expr 'ispositive(a)' -prefix MNI152maskbin_3mm+tlrc 

# save 3mm binary mask as NIFTI
3dAFNItoNIFTI -prefix MNI152maskbin_3mm MNI152maskbin_3mm+tlrc

# delete dataset
rm s50433_func_preproc.nii.gz