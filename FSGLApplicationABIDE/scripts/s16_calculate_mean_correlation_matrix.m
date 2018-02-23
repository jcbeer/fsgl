%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Take the mean of all pairwise correlation matrices 
% in the training set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% fileids_training.txt
% ${subj}_pairwise_corr.1D
%%% OUTPUTS:
% pairwisecorr_voxel_indices.txt
% pairwisecorr_mean.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('./ABIDE/data/')
% read in training ids
fileID = fopen('fileids_training.txt', 'r');
trainids = textscan(fileID, '%s');
fclose(fileID);
% navigate to subject data directory
cd('./subjectdata/train/')
% create empty array
pairwisecorr = zeros(6630, 6630, 175);
for i = 1:length(trainids{1})
    % define filename of subject's pairwise voxel correlation matrix
    filename = fullfile(trainids{1}(i), strcat(trainids{1}(i), '_pairwise_corr.1D'));
    % read in subject's data
    fileID = fopen(char(filename), 'r');
    data = textscan(fileID, '%f', 'HeaderLines', 5);
    fclose(fileID);
    % convert to matrix
    datamatrix = vec2mat(data{1}, 6631);
    % omit first column (row names) and
    % store in the big array
    pairwisecorr(:,:,i) = datamatrix(:,2:6631);
end

clear data;
% save the indices to file
dlmwrite('pairwisecorr_voxel_indices.txt', datamatrix(:,1), 'delimiter', ',', 'precision', 9);
clear datamatrix

% calculate mean of pairwisecorr
% WARNING: computationally / memory intensive 
% probably more efficient method exists, e.g. accumarray function
pairwisecorr_mean = mean(pairwisecorr, 3);

% save to a csv file
% navigate to subject data directory
cd('/ABIDE/data/')
dlmwrite('pairwisecorr_mean.txt', pairwisecorr_mean, 'delimiter', ',', 'precision', '%.6f');