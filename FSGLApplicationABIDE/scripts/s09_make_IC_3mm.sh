#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Resample the independent components (4mm voxels) to 3mm voxels 
# to match the resting state data
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# MNI152mask_3mm+tlrc
# melodic_IC.nii (4mm voxels)
### OUTPUTS:
# melodic_IC_3mm+tlrc
# melodic_IC_3mm.nii
####################################################################

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

######################################################
# resample IC components (4mm) to 3mm voxels
######################################################
cd '/ABIDE/data'

# resample to 3mm
3dresample -master MNI152mask_3mm+tlrc -prefix melodic_IC_3mm -input melodic_IC.nii

# convert 3mm file to NIFTI
3dAFNItoNIFTI -prefix melodic_IC_3mm melodic_IC_3mm+tlrc
