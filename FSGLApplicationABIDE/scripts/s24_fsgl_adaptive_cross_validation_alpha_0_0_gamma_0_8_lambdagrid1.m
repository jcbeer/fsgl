%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Do 5-fold cross-validation on training set
% for various alpha, gamma combinations
% ADAPTIVE PENALTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% srs_train.txt
% trainXstd_5476.txt
% Kdata_adaptive.dat
% Kn.csv
% weights.csv
% betaridge.csv
% folds10.csv
% fsglfit_adaptive.m
%%% OUTPUTS:
% (for various alpha, gamma, lambdagrid combinations)
% cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat
% niterations_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
% load outcome Y -- adjusted SRS
srstrain = readtable('srs_train.txt');
Y = table2array(srstrain(:,3));
clear srstrain
% load predictor matrix X -- seed connectivity
X = table2array(readtable('trainXstd_5476.txt'));
% load K matrix
load Kdata_adaptive.dat
Kmatrix = spconvert(Kdata_adaptive);
clear Kdata_adaptive
% load other info
Kn = importdata('Kn.csv', ',');
nj = Kn(1); 
nd = Kn(2);
ngroups = Kn(3); 
groupsizes = Kn(4:size(Kn, 1));
weights = importdata('weights.csv', ',');
initialbeta = importdata('betaridge.csv', ','); % used to initialize beta
folds = importdata('folds10.csv', ',');

%% define alpha, gamma, and lambda
alphagamma = zeros(1,2);
% uncomment one of the alpha, gamma combinations 
alphagamma(1,:) = [0.0 0.8]; 
% alphagamma(1,:) = [0.2 0.8]; 
% alphagamma(1,:) = [0.2 1.0]; 
% alphagamma(1,:) = [1.0 1.0]; 
% lambda
exponent = [0:(7/47):7 7:(3/23):10];
% complete grid of lambda values
lambdagrid = (repelem(exp(1), 72) .^ exponent)';
% can divide into sets of 24 and run in parallel
% first set
lambdagrid1 = lambdagrid(1:24);
% second set
% lambdagrid2 = lambdagrid(25:48);
% third set
% lambdagrid3 = lambdagrid(49:72);

%% divide training set into 5 CV folds
folds1 = folds(:,1);

%% create empty array to store results
% dim 1 is alpha gamma values
% dim 2 is lambda values
% dim 3 are the 5 CV folds
cvresults = zeros(size(alphagamma,1), length(lambdagrid), 5);
niterations = zeros(size(alphagamma,1), length(lambdagrid), 5);

%% begin loop for CV
% runs in parallel
% set up pool
parpool('local', 24);
% can loop over alphagamma if desired 
% (currently alphagamma has only one row)
for i = 1:size(alphagamma,1)
    alpha = alphagamma(i,1);
    gamma = alphagamma(i,2);
    % loop over lambda (done in parallel)
    parfor j = 1:length(lambdagrid1) % change lambdagrid as needed
    % for j = 1:length(lambdagrid1)
        lambda = lambdagrid1(j); % change lambdagrid as needed
        % loop over folds
        for k = 1:5
            % divide data based on CV fold
            Xfoldk = X(folds1==k,:);
            Yfoldk = Y(folds1==k);
            Xfoldnotk = X(folds1~=k,:);
            Yfoldnotk = Y(folds1~=k);
            % find column means for training folds
            Xfoldnotk_means = mean(Xfoldnotk, 1);
            Yfoldnotk_mean = mean(Yfoldnotk);
            % center Xs and Ys according to mean of training folds
            Xfoldnotk_centered = detrend(Xfoldnotk, 'constant');
            Yfoldnotk_centered = detrend(Yfoldnotk, 'constant');
            Xfoldk_centered = Xfoldk - repmat(Xfoldnotk_means, length(Yfoldk), 1);
            Yfoldk_centered = Yfoldk - Yfoldnotk_mean; 
            % clear variables to save memory
            % (can't use clear inside parfor)
            Xfoldk = []; 
            Yfoldk = []; 
            Xfoldnotk = [];
            Yfoldnotk = []; 
            Xfoldnotk_means = []; 
            Yfoldnotk_mean = [];
            % fit model to all folds but k
            % and get estimated beta
            [~, ~, beta, ~, ~, niter] = fsglfit_adaptive(...
                Xfoldnotk_centered,...
                Yfoldnotk_centered,...
                Kmatrix,...
                nj,... 
                nd,... 
                ngroups,... 
                groupsizes',...
                weights,...
                initialbeta,...
                alpha,... 
                gamma,...
                lambda,...
                1000,...    % rho (step size)
                6000);      % max ADMM iterations);
            % calculate the mse on the fold left out
            cv_mse_y = (1/35)*sum((Yfoldk_centered - Xfoldk_centered * beta).^2);
            % save cross validation error
            cvresults(i,j,k) = cv_mse_y;
            % save number of ADMM iterations
            niterations(i,j,k) = niter;
            % print progress to console
            % ijk = [i j k];
            % disp(ijk)
        end
    end
end
delete(gcp('nocreate'))


%% save output 
% change filename to match alpha, gamma, lambdagrid used
save('cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid1.mat', 'cvresults')
save('niterations_adaptive_alpha_0_2_gamma_0_8_lambdagrid1.mat', 'niterations')

