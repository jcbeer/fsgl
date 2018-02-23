!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Calculate pairwise correlation between all voxel time series 
# in the cortical mask for each subject in the training set. 
# Will be used to define voxel groups.
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# cortical_mask+tlrc
# ${subj}_func_preproc.nii.gz: individual subject fMRI data
### OUTPUTS: For each training subject:
# ${subj}_pairwise_corr.1D
####################################################################

# navigate to top directory
cd '/ABIDE/data'

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

# set cortical mask file
set corticalmask = cortical_mask

#######################################################
# READ IN SUBJECT IDs FOR TRAIN DATA
#######################################################
set train=`cat fileids_training.txt`

#######################################################
# Make Correlation Matrices for TRAINING DATA
#######################################################
# navigate to directory 
cd subjectdata/train

# begin loop over subjects
foreach subj ($train)
  
# move to the subject subdirectory 
cd $subj

# make pairwise correlation matrix 
3dAutoTcorrelate -pearson -mask ../../../$corticalmask+tlrc -mask_source ../../../$corticalmask+tlrc -out1D ${subj}_pairwise_corr.1D ${subj}_func_preproc.nii.gz

# delete the BRIK and HEAD files that were created
rm ATcorr+tlrc.BRIK ATcorr+tlrc.HEAD

# go back up a level
cd ..

end