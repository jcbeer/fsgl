%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Make binary mask vector csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% MNI152maskbin_3mm.nii
%%% OUTPUTS:
% MNI152maskbinary3mm.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set directory
cd('./data')

% read the 3mm binary MNI152 mask map 
V=spm_vol('MNI152maskbin_3mm.nii');
mask=spm_read_vols(V(1));
maskvect = reshape(mask, [1,271633]);
tabulate(maskvect)

% save to a csv file
csvwrite('MNI152maskbinary3mm.csv', maskvect);