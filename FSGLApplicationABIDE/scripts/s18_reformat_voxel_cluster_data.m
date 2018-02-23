%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Transform voxel clusters into 3D nii file for viewing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% diagcorr.txt
% voxel_clusters.txt
% ICgroups.nii
%%% OUTPUTS:
% diagcorr.nii: 3D map of voxels used in the new cortical mask
% cortical_mask_new.txt: text file of new cortical mask
% voxel_clusters_50.nii: 3D map of 50 voxel groups
% voxel_clusters_50.txt: text file of 50 voxel groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('./ABIDE/data/')
diagcorr = readtable('diagcorr.txt');
voxel_clusters = readtable('voxel_clusters.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note about indices created by AFNI program
% 3dAutoTcorrelate:
% 1D index (ijk) is computed from the 3D (i,j,k) indices:
% ijk = i + j*Ni + k*Ni*Nj , with Ni and Nj being the
% number of voxels in the slice orientation and given by:
% 3dinfo -ni -nj YOUR_VOLUME_HERE
% ni = 61, nj = 73, nk = 61, ni*nj*nk = 271633
% turns out we need to add slice in the x dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFORMAT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create zero vectors the size of volume
diagcorr_all = zeros(271633, 1);
voxel_clusters_50 = zeros(271633, 1);

% fill in with the data
% rescale the diagcorr to be between 0 and 1000
diagcorr_all(diagcorr{:,1},1) = diagcorr{:,2}*1000;
voxel_clusters_50(voxel_clusters{:,1},1) = voxel_clusters{:,2};

% reshape into a 3D array
diagcorr_3d = reshape(diagcorr_all,61,73,61);
% add two slices 73*61 y*z onto the side;
% and remove two slices from the other side
diagcorr_3d_new = zeros(63,73,61);
diagcorr_3d_new(3:63,:,:) = diagcorr_3d; 
diagcorr_3d_new = diagcorr_3d_new(2:62,:,:);

% reshape into a 3D array
voxel_clusters_50_3d = reshape(voxel_clusters_50,61,73,61);
% add two slices 73*61 y*z onto the side;
% and remove two slices from the other side
voxel_clusters_50_new = zeros(63,73,61);
voxel_clusters_50_new(3:63,:,:) = voxel_clusters_50_3d; 
voxel_clusters_50_new = voxel_clusters_50_new(2:62,:,:);

% convert diagcorr_3d_new to nii file
% take header information from a previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'diagcorr.nii';  
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname;  
% write out the new data
% give spm_write_vol the new header information
% and corresponding data matrix
spm_write_vol(HeaderInfo, diagcorr_3d_new);

% convert voxel_clusters_50_new to nii file
% take header information from a previous file with similar dimensions 
% and voxel sizes and change the filename in the header
HeaderInfo = spm_vol('ICgroups.nii');
% fill in the new filename
HeaderInfo.fname = 'voxel_clusters_50.nii';  
% replace the old filename in another location within the header
HeaderInfo.private.dat.fname = HeaderInfo.fname; 
% write out the new data
% give spm_write_vol the new header information
% and corresponding data matrix
spm_write_vol(HeaderInfo, voxel_clusters_50_new);

% save data as text files
% new cortical mask
cortical_mask = reshape(voxel_clusters_50_new, 1, []) > 0;
dlmwrite('cortical_mask_new.txt', cortical_mask, 'delimiter', ',');
% voxel clusters -- mask out zeros first
clust50 = reshape(voxel_clusters_50_new, 1, []);
clust50 = clust50(cortical_mask==1);
dlmwrite('voxel_clusters_50.txt', clust100, 'delimiter', ',');