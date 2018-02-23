%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Save mat phenotype data as csv
% for ABIDE subjects used in Cerliani et al paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% pheno_359.mat 
% downloaded from https://github.com/sblnin/rsfnc
%%% OUTPUTS:
% phenomat.csv: matrix of phenotype data for 359 subjects
% phenolabels.csv: header labels for phenomat.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ./data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% export .mat phenotype data as .csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load pheno_359.mat
csvwrite('phenomat.csv', phenomat)
% export header labels
fid = fopen('phenolab.csv', 'w');
fprintf(fid, '%s,', phenolabels{1,1:end-1});
fprintf(fid, '%s\n', phenolabels{1,end});
fclose(fid);