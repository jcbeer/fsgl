%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Reformat the seed-based correlation data into text files
% One for training set (n = 175), one for test set (n = 44)
% 5476 voxels in the cortical brain mask
% Also standardize the columns for use as predictor matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% fileids_training.txt
% fileids_test.txt
% cortical_mask_new.txt
% ${subj}_corr_pc1_z.nii
%%% OUTPUTS:
% trainX_5476.txt
% testX_5476.txt
% trainXstd_5476.txt
% testXstd_5476.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('./ABIDE/data/')
% read in training ids
fileID = fopen('fileids_training.txt', 'r');
trainids = textscan(fileID, '%s');
fclose(fileID);
% read in test ids
fileID = fopen('fileids_test.txt', 'r');
testids = textscan(fileID, '%s');
fclose(fileID);
% read in the new cortical mask (5476 voxels)
cortical_mask = csvread('cortical_mask_new.txt');

%% TRAINING SET
% navigate to subject data directory
cd('./subjectdata/train/')
% create empty array
trainX = zeros(175, 5476);
for i = 1:length(trainids{1})
    % define filename of subject's seed connectivity z-stat matrix
    filename = fullfile(trainids{1}(i), strcat(trainids{1}(i), '_corr_pc1_z.nii'));
    % read in subject's data
    V = spm_vol(filename);
    map = spm_read_vols(V{1});
    mapvector = reshape(map, [1, 271633]);
    maskedmapvector = mapvector(1,cortical_mask==1);
    trainX(i,:) = maskedmapvector;
end

%% TEST SET
% navigate to subject data directory
cd('./subjectdata/test/')
% create empty array
testX = zeros(44, 5476);
for i = 1:length(testids{1})
    % define filename of subject's seed connectivity z-stat matrix
    filename = fullfile(testids{1}(i), strcat(testids{1}(i), '_corr_pc1_z.nii'));
    % read in subject's data
    V = spm_vol(filename);
    map = spm_read_vols(V{1});
    mapvector = reshape(map, [1, 271633]);
    maskedmapvector = mapvector(1,cortical_mask==1);
    testX(i,:) = maskedmapvector;
end

%% standardize the columns
% based on the training set data means and sds
traincolmeans = mean(trainX, 1);
traincolsds = std(trainX, 0, 1);
%% check if close to test column means and sds
% testcolmeans = mean(testX, 1);
% testcolsds = std(testX, 0, 1);
% meandiff = traincolmeans - testcolmeans;
% sddiff = traincolsds - testcolsds;
% histogram(meandiff);
% histogram(sddiff);
%% standardize
trainXstd = zscore(trainX);
traincolmeanmatrix = repmat(traincolmeans, 44, 1);
traincolsdmatrix = repmat(traincolsds, 44, 1);
testXstd = (testX - traincolmeanmatrix) ./ traincolsdmatrix;

%% save to csv file
% navigate to data directory
cd('/ABIDE/data/')
dlmwrite('trainX_5476.txt', trainX, 'delimiter', ',', 'precision', '%.6f');
dlmwrite('testX_5476.txt', testX, 'delimiter', ',', 'precision', '%.6f');
dlmwrite('trainXstd_5476.txt', trainXstd, 'delimiter', ',', 'precision', '%.6f');
dlmwrite('testXstd_5476.txt', testXstd, 'delimiter', ',', 'precision', '%.6f');