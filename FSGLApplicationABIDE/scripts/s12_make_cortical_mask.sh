#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Make cortical brain mask 
# of the 2 sensorimotor cortical independent components
# volume 4: IC5 -- dorsal sensorimotor
# volume 18: IC29 -- ventral sensorimotor
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# ICgroups.nii
### OUTPUTS:
# cortical_mask+tlrc
####################################################################

# navigate to working directory
cd '/ABIDE/data'

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

3dcalc -a ICgroups.nii -expr 'ispositive(-(a-4)^2+0.5) + ispositive(-(a-18)^2+0.5)' -prefix cortical_mask
