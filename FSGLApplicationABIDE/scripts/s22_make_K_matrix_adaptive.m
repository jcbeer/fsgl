%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Make the K matrix for adaptive penalties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% cortical_mask_new.txt
% voxel_clusters_50.txt
% makeKmatrix_adaptive.m: Matlab function to make K matrix, 
% also calculates weights
%%% OUTPUTS:
% NOTE adaptive version of the K matrix is the same as non-adaptive version
% Kdata_adaptive.dat: sparse K matrix
% Kmatrix_adaptive.csv: sparse K matrix
% Kn.csv: additional data on group sizes for fitting FSGL
% weights.csv: vector of weights for adaptive FSGL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('./ABIDE/data/')
% read in cortical mask
cortical_mask = csvread('cortical_mask_new.txt');
% read in voxel group data
voxelgps50 = csvread('voxel_clusters_50.txt');
% read in the weights
betaridge = csvread('betaridge.csv');

% make K matrix
[K, nj, nd, ngroups, groupsizes, weights] = makeKmatrix_adaptive(...
    61, 73, 61, cortical_mask, voxelgps50, betaridge);

% save Kmatrix as csv
[i, j, val] = find(K);
Kdata = [i, j, val];
% add row at the bottom specifying size 
Kdata = [Kdata ; [size(K, 1), size(K, 2), 0]];
% reconvert using K = spconvert(Kdata);
writetable(array2table(Kdata), 'Kdata_adaptive.dat', 'WriteVariableNames', false)
csvwrite('Kmatrix_adaptive.csv', Kdata)

% save Kn as csv
Kn = [nj ; nd ; ngroups ; groupsizes];
dlmwrite('Kn.csv', Kn, 'delimiter', ',', 'precision', 9);
dlmwrite('weights.csv', weights, 'delimiter', ',', 'precision', 9);