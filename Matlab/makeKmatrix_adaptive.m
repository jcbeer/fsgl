%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKEKMATRIX_ADAPTIVE
% Function to create K matrix and generate adaptive 
% weights for fused sparse group lasso adaptive
% penalties. Provides input for fsglfit_adaptive.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methodology described in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
% dim1: dimension 1 (x) of 3D image volume
% dim2: dimension 2 (y) of 3D image volume
% dim3: dimension 3 (z) of 3D image volume
% mask: vector of length dim1*dim2*dim3 that has a 1 for each predictor 
%   voxel in the 3D volume and a 0 for each voxel that is not 
%   of interest in the 3D volume (e.g. voxels outside of brain). 
%   The sum of the mask vector equals to p, the number of predictors.
%   NOTE: The 3D volume is assumed to be concatenated into a vector 
%   in the same way the reshape function works in Matlab, i.e., 
%   counting fastest along dim1, then dim2, then dim3.
% groups: vector of length p that assigns each predictor to a unique 
%   group, using integers 1, 2, etc.
% betaestimates: predetermined estimates for beta, e.g. from OLS or 
%   ridge regression results, used for defining the adaptive weights
%%% OUTPUTS
% K: K matrix 
%   (input for fsglfit_adaptive function)
% nj: number of predictors, p 
%   (input for fsglfit_adaptive function)
% nd: number of rows of D matrix 
%   (input for fsglfit_adaptive function)
% ngroups: number of coefficient groups 
%   (input for fsglfit_adaptive function) 
% groupsizes: vector of coefficient group sizes 
%   (input for fsglfit_adaptive function) 
% weights: adaptive penalty weights to rescale lambda  
%   (input for fsglfit_adaptive function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, nj, nd, ngroups, groupsizes, weights] = ...
    makeKmatrix(dim1, dim2, dim3, mask, groups, betaestimates)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make J matrix for lasso penalty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    masksize = sum(mask);
    J = sparse(1:masksize, 1:masksize, repelem(1, masksize), masksize, masksize);
    % save the number of rows of J
    nj = masksize;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make D matrix for fused penalty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dnrow = 3*dim1*dim2*dim3 - dim3*dim2 - dim3*dim1 - dim1*dim2;
    Dncol = dim1*dim2*dim3;
    % create column indices
    % ONES
    % to bind columns
    bindcolones = 1:dim1*dim2*dim3;
    bindcolonesremove = (1:(dim2*dim3))*dim1;
    bindcolones(bindcolonesremove) = [];
    % to bind rows
    bindrowones = 1:dim1*dim2*dim3;
    a = (1:dim1) + (dim1*(dim2-1));
    A = repmat(a', 1, dim3);
    b = [0, (1:(dim3-1)) * (dim1*dim2)];
    c = bsxfun(@plus, A, b);
    bindrowonesremove = reshape(c, [1, numel(c)]);
    bindrowones(bindrowonesremove) = [];
    % to bind slices
    bindsliceones = 1:dim1*dim2*dim3;
    bindsliceonesremove = ((dim1*dim2*dim3 - dim1*dim2) + 1):dim1*dim2*dim3;
    bindsliceones(bindsliceonesremove) = [];
    % NEGATIVE ONES
    bindcolnegones = bindcolones + 1;
    bindrownegones = bindrowones + dim1;
    bindslicenegones = bindsliceones + (dim1*dim2);
    % put all column indices together
    Dcolindices = [bindcolones bindrowones bindsliceones bindcolnegones bindrownegones bindslicenegones];
    % make unmasked D matrix
    Dunmasked = sparse([1:Dnrow 1:Dnrow], Dcolindices, [repelem(1, Dnrow) repelem(-1, Dnrow)], Dnrow, Dncol);
    % clear some temporary variables
    clear bindcolones bindcolnegones bindcolonesremove bindrowones bindrowneg ones bindrowonesremove bindsliceones bindslicenegones bindsliceonesremove a A b c Dncol Dnrow
    % apply the mask to the D matrix
    % remove any row that has an entry at a zero position in the mask
    % then remove all columns corresponding to zero positions in the mask
    Dmaskzerocols = abs(Dunmasked(:,mask==0));
    Dmaskzerocolsrowsum = sum(Dmaskzerocols, 2);
    % tabulate(Dmaskedcolsrowsum)
    D = Dunmasked(Dmaskzerocolsrowsum==0, mask==1);
    % save rows of D
    nd = size(D, 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make G matrix for groups
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for a complete partition of the data as in this case
    % the G matrix is just a reordering of the rows of the J matrix
    grouptab = tabulate(groups);
    % save the number of groups
    ngroups = size(grouptab, 1);
    % save the number of voxels in each group;
    groupsizes = grouptab(:,2);
    % create an ordering of the indices
    [~, groupindices] = sort(groups);
    % make the G matrix
    G = sparse(1:masksize, groupindices, repelem(1, masksize), masksize, masksize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % combine J D and G to form K
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = [J; D; G];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the weight vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate lasso weights
    % take absolute value
    absbetaestimates = abs(betaestimates);
    % take inverse
    absbetaestimatesinv = 1./absbetaestimates;
    % rescale so coefficients add to p
    lassowts = absbetaestimatesinv*(nj/sum(absbetaestimatesinv));     
    % calculate fusion weights
    fusionwts = (nj/sum(absbetaestimatesinv))./(abs(D * betaestimates));
    % calculate group weights
    groupwts = zeros(ngroups, 1);
    for i = 1:ngroups
        groupwts(i) = (nj/sum(absbetaestimatesinv))./norm(betaestimates(groups==i));
    end
    weights = [lassowts;  fusionwts; groupwts];
end