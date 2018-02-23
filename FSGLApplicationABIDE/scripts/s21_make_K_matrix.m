%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Make the K matrix for non-adaptive penalties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% cortical_mask_new.txt
% voxel_clusters_50.txt
% makeKmatrix.m: Matlab function to make K matrix
%%% OUTPUTS:
% Kdata.dat: sparse K matrix
% Kmatrix.csv: sparse K matrix
% Kn.csv: additional data on group sizes for fitting FSGL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('./ABIDE/data/')
% read in cortical mask
cortical_mask = csvread('cortical_mask_new.txt');
% read in voxel group data
voxelgps50 = csvread('voxel_clusters_50.txt');

% make K matrix
[K, nj, nd, ngroups, groupsizes] = ...
    makeKmatrix(61, 73, 61, cortical_mask, voxelgps50);

% save Kmatrix as csv
[i, j, val] = find(K);
Kdata = [i, j, val];
% add row at the bottom specifying size 
Kdata = [Kdata ; [size(K, 1), size(K, 2), 0]];
% reconvert using K = spconvert(Kdata);
writetable(array2table(Kdata), 'Kdata.dat', 'WriteVariableNames', false)
csvwrite('Kmatrix.csv', Kdata)

% save Kn as csv
Kn = [nj ; nd ; ngroups ; groupsizes];
dlmwrite('Kn.csv', Kn, 'delimiter', ',', 'precision', 9);