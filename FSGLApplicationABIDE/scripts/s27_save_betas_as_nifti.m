%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Save estimated betas for adaptive penalties as .nii files
% (Figure 6C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% cortical_mask_new.txt
% beta_hat.csv
% ICgroups.nii
%%% OUTPUTS:
% betas_alpha_1_0_gamma_1_0_ada_1000.nii
% betas_alpha_0_2_gamma_1_0_ada_1000.nii
% betas_alpha_0_2_gamma_0_8_ada_1000.nii
% betas_alpha_0_0_gamma_0_8_ada_1000.nii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
cd('./ABIDE/data/')
cortical_mask = csvread('cortical_mask_new.txt');
betas = csvread('beta_hat.csv');
% columns of beta_hat in order:
% 1. ridge estimates 
% 2. alpha 1.0 gamma 1.0
% 3. alpha 0.2 gamma 1.0
% 4. alpha 0.2 gamma 0.8
% 5. alpha 0.0 gamma 0.8
% 6. alpha 1.0 gamma 1.0 adaptive
% 7. alpha 0.2 gamma 1.0 adaptive
% 8. alpha 0.2 gamma 0.8 adaptive
% 9. alpha 0.0 gamma 0.8 adaptive


%% alpha 1.0 gamma 1.0 adaptive
% rescale by factor 1000
betas_alpha_1_0_gamma_1_0_ada_1000 = betas(:,6)*1000;
% fill in zeros for unmasked voxels in the volume
betas_alpha_1_0_gamma_1_0_ada_1000_unmasked = cortical_mask;
betas_alpha_1_0_gamma_1_0_ada_1000_unmasked(cortical_mask==1) = ...
    betas_alpha_1_0_gamma_1_0_ada_1000;
% reshape to 3D array
betas_alpha_1_0_gamma_1_0_ada_1000_3D = reshape(betas_alpha_1_0_gamma_1_0_ada_1000_unmasked, 61, 73, 61);
% convert to nii file
% take header information from previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'betas_alpha_1_0_gamma_1_0_ada_1000.nii';
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname; 
% prevent rescaling
HeaderInfo.pinfo = [1;0;0];
HeaderInfo.dt = [64 0];
% use spm_write_vol to write out the new data
% give spm_write_vol the new header information and corresponding data matrix
spm_write_vol(HeaderInfo, betas_alpha_1_0_gamma_1_0_ada_1000_3D);
% clear data
clear betas_alpha_1_0_gamma_1_0_ada_1000
clear betas_alpha_1_0_gamma_1_0_ada_1000_unmasked
clear betas_alpha_1_0_gamma_1_0_ada_1000_3D


%% alpha 0.2 gamma 1.0 adaptive
% rescale by factor 1000
betas_alpha_0_2_gamma_1_0_ada_1000 = betas(:,7)*1000;
% fill in zeros for unmasked voxels in the volume
betas_alpha_0_2_gamma_1_0_ada_1000_unmasked = cortical_mask;
betas_alpha_0_2_gamma_1_0_ada_1000_unmasked(cortical_mask==1) = ...
    betas_alpha_0_2_gamma_1_0_ada_1000;
% reshape to 3D array
betas_alpha_0_2_gamma_1_0_ada_1000_3D = ...
    reshape(betas_alpha_0_2_gamma_1_0_ada_1000_unmasked, 61, 73, 61);
% convert to nii file
% take header information from previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'betas_alpha_0_2_gamma_1_0_ada_1000.nii';
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname; 
% prevent rescaling
HeaderInfo.pinfo = [1;0;0];
HeaderInfo.dt = [64 0];
% use spm_write_vol to write out the new data
% give spm_write_vol the new header information and corresponding data matrix
spm_write_vol(HeaderInfo, betas_alpha_0_2_gamma_1_0_ada_1000_3D); 
% clear data
clear betas_alpha_0_2_gamma_1_0_ada_1000
clear betas_alpha_0_2_gamma_1_0_ada_1000_unmasked
clear betas_alpha_0_2_gamma_1_0_ada_1000_3D


%% alpha 0.2 gamma 0.8 adaptive
% rescale by factor 1000
betas_alpha_0_2_gamma_0_8_ada_1000 = betas(:,8)*1000;
% fill in zeros for unmasked voxels in the volume
betas_alpha_0_2_gamma_0_8_ada_1000_unmasked = cortical_mask;
betas_alpha_0_2_gamma_0_8_ada_1000_unmasked(cortical_mask==1) = ...
    betas_alpha_0_2_gamma_0_8_ada_1000;
% reshape to 3D array
betas_alpha_0_2_gamma_0_8_ada_1000_3D = ...
    reshape(betas_alpha_0_2_gamma_0_8_ada_1000_unmasked, 61, 73, 61);
% convert to nii file
% take header information from previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'betas_alpha_0_2_gamma_0_8_ada_1000.nii';  
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname;  
% prevent rescaling
HeaderInfo.pinfo = [1;0;0];
HeaderInfo.dt = [64 0];
% use spm_write_vol to write out the new data
% give spm_write_vol the new header information and corresponding data matrix
spm_write_vol(HeaderInfo, betas_alpha_0_2_gamma_0_8_ada_1000_3D); 
% clear data
clear betas_alpha_0_2_gamma_0_8_ada_1000
clear betas_alpha_0_2_gamma_0_8_ada_1000_unmasked
clear betas_alpha_0_2_gamma_0_8_ada_1000_3D


%% alpha 0.0 gamma 0.8 adaptive
% rescale by factor 1000
betas_alpha_0_0_gamma_0_8_ada_1000 = betas(:,9)*1000;
% fill in zeros for unmasked voxels in the volume
betas_alpha_0_0_gamma_0_8_ada_1000_unmasked = cortical_mask;
betas_alpha_0_0_gamma_0_8_ada_1000_unmasked(cortical_mask==1) = ...
    betas_alpha_0_0_gamma_0_8_ada_1000;
% reshape to 3D array
betas_alpha_0_0_gamma_0_8_ada_1000_3D = ...
    reshape(betas_alpha_0_0_gamma_0_8_ada_1000_unmasked, 61, 73, 61);
% convert to nii file
% take header information from previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'betas_alpha_0_0_gamma_0_8_ada_1000.nii'; 
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname;  
% prevent rescaling
HeaderInfo.pinfo = [1;0;0];
HeaderInfo.dt = [64 0];
% use spm_write_vol to write out the new data
% give spm_write_vol the new header information and corresponding data matrix
spm_write_vol(HeaderInfo, betas_alpha_0_0_gamma_0_8_ada_1000_3D); 
% clear data
clear betas_alpha_0_0_gamma_0_8_ada_1000
clear betas_alpha_0_0_gamma_0_8_ada_1000_unmasked
clear betas_alpha_0_0_gamma_0_8_ada_1000_3D