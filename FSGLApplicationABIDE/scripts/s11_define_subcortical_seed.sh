#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Make a seed consisting of peak voxels of 
# IC17 from Cerliani et al 2015 paper (Figure 5A)
# This is coded as subbrick 16 in the melodic_IC_3mm file
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# melodic_IC_3mm+tlrc.HEAD
### OUTPUTS:
# subcortical_seed_mask+tlrc
# subcortical_seed_mask_binary+tlrc
####################################################################

# navigate to working directory
cd '/ABIDE/data'

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

# save mask of peak voxels -- 2 clusters
3dclust -1Dformat -nosum -1dindex 16 -1tindex 16 -2thresh -20.27 20.27 -dxyz=1 -savemask subcortical_seed_mask 1.01 2 melodic_IC_3mm+tlrc.HEAD

# make a binary mask of these two clusters
3dcalc -a subcortical_seed_mask+tlrc -expr 'step(a)' -prefix subcortical_seed_mask_binary