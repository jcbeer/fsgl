#!/bin/tcsh

####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Make connectivity maps for voxels in cortical mask
# based on correlation of the 1st PC of the subcortical seed 
# voxel time series with all other voxels in the mask
# (Figure 5B)
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# subcortical_seed_mask_binary+tlrc
# cortical_mask+tlrc
# ${subj}_func_preproc.nii.gz: individual subject fMRI data
### OUTPUTS: For each subject:
# ${subj}_subcortical_pc1_vec.1D: first PC time series
# ${subj}_corr_pc1+tlrc: Pearson correlation map, AFNI format
# ${subj}_corr_pc1_z+tlrc: correlation z-map, AFNI format
# ${subj}_corr_pc1_z.nii: correlation z-map, NIFTI format
####################################################################

# navigate to top directory
cd '/ABIDE/data'

# set some variables to run afni
setenv DYLD_LIBRARY_PATH /opt/X11/lib/flat_namespace
setenv DYLD_FALLBACK_LIBRARY_PATH $HOME/abin

# set seed mask file
set seedmask = subcortical_seed_mask_binary

# set cortical mask file
set corticalmask = cortical_mask

#######################################################
# READ IN SUBJECT IDs FOR TRAIN AND TEST DATA
#######################################################
set train=`cat fileids_training.txt`
set test=`cat fileids_test.txt`

#######################################################
# TRAINING DATA
#######################################################
# navigate to directory 
cd subjectdata/train

# begin loop over subjects
foreach subj ($train)
  
# move to the subject subdirectory 
cd $subj

# extract first PC of the time series in the seed mask
3dpc -nscale -pcsave 1 -prefix ${subj}_subcortical_pc1 -mask ../../../$seedmask+tlrc ${subj}_func_preproc.nii.gz

# make correlation map with the first PC
3dTcorr1D -pearson -prefix ${subj}_corr_pc1 -mask ../../../$corticalmask+tlrc ${subj}_func_preproc.nii.gz ${subj}_subcortical_pc1_vec.1D

# do Fisher r-to-z transformation for first PC
3dcalc -a ${subj}_corr_pc1+tlrc -expr '0.5*(log((1+a)/(1-a)))' -prefix ${subj}_corr_pc1_z

# convert the 3mm PC1 z-map to a nii file 
# (will use this for the predictor matrix)
3dAFNItoNIFTI -prefix ${subj}_corr_pc1_z ${subj}_corr_pc1_z+tlrc

# go back up a level
cd ..

end


#######################################################
# TEST DATA
#######################################################
# navigate to directory 
cd ../test

# begin loop over subjects
foreach subj ($test)
  
# move to the subject subdirectory 
cd $subj

# extract first PC of the time series in the mask
3dpc -nscale -pcsave 1 -prefix ${subj}_subcortical_pc1 -mask ../../../$seedmask+tlrc ${subj}_func_preproc.nii.gz

# make correlation map with the first PC
3dTcorr1D -pearson -prefix ${subj}_corr_pc1 -mask ../../../$corticalmask+tlrc ${subj}_func_preproc.nii.gz ${subj}_subcortical_pc1_vec.1D

# do Fisher r-to-z transformation for first PC
3dcalc -a ${subj}_corr_pc1+tlrc -expr '0.5*(log((1+a)/(1-a)))' -prefix ${subj}_corr_pc1_z

# convert the PC1 z-map to a nii file 
# (will use this for the predictor matrix)
3dAFNItoNIFTI -prefix ${subj}_corr_pc1_z ${subj}_corr_pc1_z+tlrc

# go back up a level
cd ..

end