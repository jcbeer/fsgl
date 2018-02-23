# Fused Sparse Group Lasso ABIDE Application

As reported in the manuscript *Incorporating Prior Information with Fused Sparse Group Lasso: Application to Prediction of Clinical Measures from Neuroimages*. This analysis uses a combination of R, Matlab (including SPM functions), and AFNI (tcsh) software.

## Workflow
### 1. s01_save_pheno_mat_as_csv.m
Save mat phenotype data as csv for ABIDE subjects (n = 359) used in Cerliani et al (2015) paper.

**Input:** pheno_359.mat, downloaded from https://github.com/sblnin/rsfnc \
**Output:** phenomat.csv: matrix of phenotype data for 359 subjects \
phenolabels.csv: header labels for phenomat.csv


### 2. s02_save_file_IDs.R
Save phenotype data and subject file IDs for the subset (n = 219) with nonmissing SRS scores.

**Input:** phenomat.csv, phenolab.csv, Phenotypic_V1_0b_preprocessed1.csv downloaded from http://preprocessed-connectomes-project.org/abide/download.html \
**Output:** phenosrs.csv: phenotype data for 219 subjects with nonmissing SRS score \
phenosrsall.csv: more phenotype data for 219 subjects with nonmissing SRS score \
fileids.csv: file IDs for 219 subjects with nonmissing SRS score


### 3. s03_download_data.R
Download fMRI data from ABIDE website.

**Input:** fileids.csv, phenosrsall.csv \
**Output:** Downloads fMRI data into individual subject directories for 219 subjects.


### 4. s04_assign_training_test_set.R
Randomly assign subjects (n = 219) to training (n = 175) and test (n = 44) set.

**Input:** fileids.csv \
**Output:** fileids_training.txt, fileids_test.txt


### 5. s05_make_train_test_dir.sh
Move training and test data into different parent directories.

**Input:** fileids_training.txt, fileids_test.txt \
**Output:** Moves subject fMRI data into train and test parent directories.


### 6. s06_adjust_SRS.R
Calculate descriptive statistics (Table S11). Adjust the SRS scores using linear regression (Table S12).

**Input:** phenosrsall.csv, fileids_training.txt, fileids_test.txt \
**Output:** srs_train.txt: original and adjusted SRS for train data \
srs_test.txt: original and adjusted SRS for test data


### 7. s07_makeMNI152_brain_mask.sh
Make MNI152 brain mask using $3 \times 3 \times 3$ mm$^3$ voxel template from subject resting state dataset.

**Input:** MNI152_T1_2009c+tlrc, from AFNI files \
s50433_func_preproc.nii.gz, one subject’s fMRI data \
**Output:** MNI152mask_3mm+tlrc: 3mm voxel brain mask, AFNI format \
MNI152maskbin_3mm+tlrc: 3mm voxel binary brain mask, AFNI format \
MNI152maskbin_3mm.nii: 3mm voxel binary brain mask, NIFTI format


### 8. s08_make_binary_mask_vector.m
Make binary mask vector csv file.

**Input:** MNI152maskbin_3mm.nii \
**Output:** MNI152maskbinary3mm.csv


### 9. s09_make_IC_3mm.sh
Resample the independent components (ICs, 4mm voxels) to 3mm voxels to match the resting state data.

**Input:** MNI152mask_3mm+tlrc, melodic_IC.nii (4mm voxels) \
**Output:** melodic_IC_3mm+tlrc, melodic_IC_3mm.nii


### 10. s10_assign_IC_groups.m
Assign each voxel to the independent component of maximal value.

**Input:** melodic_IC_3mm.nii \
**Output:** ICgroups.nii, ICgroups.csv


### 11. s11_define_subcortical_seed.sh
Make a seed consisting of peak voxels of IC17 from Cerliani et al (2015) paper. This is coded as subbrick 16 in the melodic_IC_3mm file. (Figure 5A)

**Input:** melodic_IC_3mm+tlrc.HEAD \
**Output:** subcortical_seed_mask+tlrc, subcortical_seed_mask_binary+tlrc


## 12. s12_make_cortical_mask.sh
Make cortical brain maskof the 2 sensorimotor cortical independent components: \
volume 4: IC5 -- dorsal sensorimotor \
volume 18: IC29 -- ventral sensorimotor \

**Input:** ICgroups.nii \
**Output:** cortical_mask+tlrc


### 13. s13_make_seed_based_connectivity_maps.sh
Make connectivity maps for voxels in cortical mask based on correlation of the 1st principal component (PC) of the subcortical seed voxel time series with all other voxels in the mask. (Figure 5B)


**Input:**  subcortical_seed_mask_binary+tlrc, cortical_mask+tlrc, \
\${subj}_func_preproc.nii.gz: individual subject fMRI data \
**Output:** For each subject (where \${subj} represents subject ID): \
\${subj}_subcortical_pc1_vec.1D: first PC time series \
\${subj}_corr_pc1+tlrc: Pearson correlation map, AFNI format \
\${subj}_corr_pc1_z+tlrc: correlation z-map, AFNI format \
\${subj}_corr_pc1_z.nii: correlation z-map, NIFTI format


### 14. s14_reformat_seed_corr_data.m
Reformat the seed-based correlation data into text files; one for training set (n = 175), one for test set (n = 44). There are 69880 voxels in the whole brain mask, 5476 in the cortical mask. Also standardize the columns for use as predictor matrix.

**Input:** fileids_training.txt, fileids_test.txt, MNI152maskbinary3mm.csv, \${subj}_corr_pc1_z.nii \
**Output:** trainX_5476.txt: training set data \
testX_5476.txt: test set data \
trainXstd_5476.txt: training set standardized data, used as predictor matrix \
testXstd_5476.txt: test set standardized data, used as predictor matrix


### 15. s15_pairwise_voxel_correlation.sh
Calculate pairwise correlation between all voxel time series in the cortical mask for each subject in the training set. Will be used to define voxel groups.

**Input:**  cortical_mask+tlrc, ${subj}_func_preproc.nii.gz (individual subject fMRI data) \
**Output:** \${subj}_pairwise_corr.1D


### 16. s16_calculate_mean_correlation_matrix.m
Take the mean of all pairwise correlation matrices in the training set.

**Input:** fileids_training.txt, \${subj}_pairwise_corr.1D \
**Output:** pairwisecorr_voxel_indices.txt, pairwisecorr_mean.txt


### 17. s17_cluster_voxels.R
Do heirarchical clustering on the average pairwise correlation matrix to form voxel groups, 50 clusters (used in final analysis, Figure 5C) and 100 clusters.

**Input:**  pairwisecorr_mean.txt, pairwisecorr_voxel_indices.txt \
**Output:** diagcorr.txt: indicates which voxels were used (no subjects were missing data) \
voxel_clusters.txt: voxel indices and voxel groups


### 18. s18_reformat_voxel_cluster_data.m
Transform voxel clusters into 3D .nii file for viewing.

**Input:**  diagcorr.txt, voxel_clusters.txt, ICgroups.nii \
**Output:** diagcorr.nii: 3D map of voxels used in the new cortical mask \
cortical_mask_new.txt: text file of new cortical mask (uses only voxels where no training subjects had missing data) \
voxel_clusters_50.nii: 3D map of 50 voxel groups \
voxel_clusters_50.txt: text file of 50 voxel groups


### 19. s19_define_CV_folds.R
Make 10 sets of 5-fold cross validation assignments stratified by adjusted SRS.

**Input:**  srs_train.txt \
**Output:**  folds10.csv


### 20. s20_ridge_and_elastic_net_regression.R
Do ridge and elastic net regression on ABIDE training data. Save ridge coefficients for making adaptive penalty weights. Calculate Sum of Squares Total (SST) for training and test sets. Results reported in Table 3.

**Input:** folds10.csv, srs_train.txt, trainXstd_5476.txt, srs_test.txt, testXstd_5476.txt \
**Output:** betaridge.csv, results for ridge and elastic net regression and SST reported in Table 3


### 21. s21_make_K_matrix.m
Make the K matrix for non-adaptive penalties.

**Input:** cortical_mask_new.txt, voxel_clusters_50.txt, \
makeKmatrix.m: Matlab function to make K matrix \
**Output:**  Kdata.dat: sparse K matrix \
Kmatrix.csv: sparse K matrix \
Kn.csv: additional data on group sizes for fitting FSGL


### 22. s22_make_K_matrix_adaptive.m
Make the K matrix for adaptive penalties. (This is actually the same as the K matrix for non-adaptive penalties. However the function makeKmatrix_adaptive.m also calculates weights.)

**Input:** cortical_mask_new.txt, voxel_clusters_50.txt, \
makeKmatrix_adaptive.m: Matlab function to make K matrix, also calculates weights \
**Output:** NOTE adaptive version of the K matrix is the same as non-adaptive version. \
Kdata_adaptive.dat: sparse K matrix \
Kmatrix_adaptive.csv: sparse K matrix \
Kn.csv: additional data on group sizes for fitting FSGL \
weights.csv: vector of weights for adaptive FSGL


### 23. s23_fsgl_cross_validation_alpha_0_0_gamma_0_8_lambdagrid1.m
Do 5-fold cross-validation on training set for various alpha, gamma combinations, non-adaptive penalties. NOTE: This script is modified to cover all alpha, gamma, and lambda values used in cross-validation. epsilon_abs and epsilon_rel were set to 10 ^ -3.

**Input:**  srs_train.txt, trainXstd_5476.txt, Kdata.dat, Kn.csv, folds10.csv, fsglfit.m \
**Output:** (for the various alpha, gamma, lambdagrid combinations) \
cvresults_alpha_0_0_gamma_0_8_lambdagrid1.mat \
niterations_alpha_0_0_gamma_0_8_lambdagrid1.mat


### 24. s24_fsgl_adaptive_cross_validation_alpha_0_0_gamma_0_8_lambdagrid1.m
Do 5-fold cross-validation on training set for various alpha, gamma combinations, adaptive penalties. NOTE: This script is modified to cover all alpha, gamma, and lambda values used in cross-validation. epsilon_abs and epsilon_rel were set to 10 ^ -3.

**Input:** srs_train.txt \
trainXstd_5476.txt \
Kdata_adaptive.dat \
Kn.csv \
weights.csv \
betaridge.csv \
folds10.csv \
fsglfit_adaptive.m \
**Output:**  (for the various alpha, gamma, lambdagrid combinations) \
cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat \
niterations_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat


### 25. s25_analyze_cv_results.m
Analyze cross-validation results. Produce cross-validation error curve figure (Figure 6A) and optimal lambda and CVMSE results reported in Table 3.

**Input:** cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat \
cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid2.mat \
cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid3.mat \
cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid1.mat \
cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid2.mat \
cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid3.mat \
cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid1.mat \
cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid2.mat \
cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid3.mat \
cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid1.mat \
cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid2.mat \
cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid3.mat \
cvresults_alpha_0_0_gamma_0_8_lambdagrid1.mat \
cvresults_alpha_0_0_gamma_0_8_lambdagrid2.mat \
cvresults_alpha_0_0_gamma_0_8_lambdagrid3.mat \
cvresults_alpha_0_2_gamma_0_8_lambdagrid1.mat \
cvresults_alpha_0_2_gamma_0_8_lambdagrid2.mat \
cvresults_alpha_0_2_gamma_0_8_lambdagrid3.mat \
cvresults_alpha_1_0_gamma_1_0_lambdagrid1.mat \
cvresults_alpha_1_0_gamma_1_0_lambdagrid2.mat \
cvresults_alpha_1_0_gamma_1_0_lambdagrid3.mat \
cvresults_alpha_0_2_gamma_1_0_lambdagrid1.mat \
cvresults_alpha_0_2_gamma_1_0_lambdagrid2.mat \
cvresults_alpha_0_2_gamma_1_0_lambdagrid3.mat \
**Output:**  pdf figure: cross-validation curves \
optimal lambda and CVMSE reported in Table 3


### 26. s26_fit_to_training_and_test_data.m
Fit non-adaptive and adaptive fused sparse group lasso to entire training set at the optimum lambda values and calculate test set predicted values. epsilon_abs and epsilon_rel were set to 10 ^ -10.

**Input:**  srs_train.txt, trainXstd_5476.txt, srs_test.txt, testXstd_5476.txt, Kdata_adaptive.dat, Kn.csv, weights.csv, betaridge.csv \
**Output:** beta_hat.csv: estimated beta coefficients \
columns of beta_hat in order (left to right): \
ridge estimates  \
alpha 1.0 gamma 1.0 \
alpha 0.2 gamma 1.0 \
alpha 0.2 gamma 0.8 \
alpha 0.0 gamma 0.8 \ 
alpha 1.0 gamma 1.0 adaptive \
alpha 0.2 gamma 1.0 adaptive \
alpha 0.2 gamma 0.8 adaptive \
alpha 0.0 gamma 0.8 adaptive \
y_hat.csv: estimated y for test set \
columns of y_hat in order (left to right): \
actual Y (adjusted SRS) \
ridge estimates \ 
alpha 1.0 gamma 1.0 \
alpha 0.2 gamma 1.0 \
alpha 0.2 gamma 0.8 \
alpha 0.0 gamma 0.8 \
alpha 1.0 gamma 1.0 adaptive \
alpha 0.2 gamma 1.0 adaptive \
alpha 0.2 gamma 0.8 adaptive \
alpha 0.0 gamma 0.8 adaptive


### 27. s27_save_betas_as_nifti.m
Save estimated betas for adaptive penalties as .nii files.

**Input:** cortical_mask_new.txt, beta_hat.csv, ICgroups.nii \
**Output:**  betas_alpha_1_0_gamma_1_0_ada_1000.nii \
betas_alpha_0_2_gamma_1_0_ada_1000.nii \
betas_alpha_0_2_gamma_0_8_ada_1000.nii \
betas_alpha_0_0_gamma_0_8_ada_1000.nii


### 28. s28_save_test_dx_data.R
Save diagnosis info for test subjects for plotting.

**Input:** phenosrsall.csv, srs_test.txt \
**Output:** phenotest.txt


### 29. s29_prediction_error_plots.m
Plot Predicted vs. Actual Adjusted SRS for test data for adaptive penalties (Figure 6B). Calculate statistics reported in Table 3.

**Input:** beta_hat.csv, y_hat.csv, phenotest.txt, srs_train.txt, trainXstd_5476.txt, srs_test.txt, testXstd_5476.txt \
**Output:** pdf figure, predicted vs. actual adjusted SRS \
statistics reported in Table 3