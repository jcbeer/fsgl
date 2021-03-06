%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Fit non-adaptive and adaptive FSGL 
% to entire training set 
% at the optimum lambda values
% and calculate test set predicted values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% srs_train.txt
% trainXstd_5476.txt
% srs_test.txt
% testXstd_5476.txt
% Kdata_adaptive.dat
% Kn.csv
% weights.csv
% betaridge.csv
%%% OUTPUTS:
% beta_hat.csv: estimated beta coefficients
% y_hat.csv: estimated y for test set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
cd('./ABIDE/data/')
% TRAINING DATA: load outcome Y -- adjusted SRS
srstrain = readtable('srs_train.txt');
Y = table2array(srstrain(:,3));
clear srstrain
% TRAINING DATA: load predictor matrix X -- seed connectivity
X = table2array(readtable('trainXstd_5476.txt'));
% TEST DATA: load outcome Y -- adjusted SRS
srstest = readtable('srs_test.txt');
Ytest = table2array(srstest(:,3));
clear srstest
% TEST DATA: load predictor matrix X -- seed connectivity
Xtest = table2array(readtable('testXstd_5476.txt'));
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
initialbeta = importdata('betaridge.csv', ',');


%% run FSGL fit NON-ADAPTIVE penalties
% define alpha, gamma, and lambda
%% lasso
alpha = 1.0;
gamma = 1.0;
lambda = 1847.8;
% run FSGL fit
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit(...
    X,...
    Y,...
    Kmatrix,...
    nj,...
    nd,...
    ngroups,...
    groupsizes',...
    alpha,...
    gamma,...
    lambda,...
    1000,... % rho (initial ADMM step size)
    6000);   % max number of ADMM iterations
% toc
% save estimated beta
beta_1_0_1_0 = beta;
% save predicted y for test data
pred_y_test_1_0_1_0 = Xtest * beta;
% save output 
save('alpha_1_0_gamma_1_0_lambda_1847.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% sparse group lasso
% define alpha, gamma, and lambda
alpha = 0.2;
gamma = 1.0;
lambda = 520.8;
% run FSGL fit
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit(...
    X,...
    Y,...
    Kmatrix,...
    nj,...
    nd,...
    ngroups,...
    groupsizes',...
    alpha,...
    gamma,...
    lambda,...
    1000,... % rho (initial ADMM step size)
    6000);   % max number of ADMM iterations
% toc
% save estimated beta
beta_0_2_1_0 = beta;
% save predicted y for test data
pred_y_test_0_2_1_0 = Xtest * beta;
% save output 
save('alpha_0_2_gamma_1_0_lambda_521.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% fused sparse group lasso
% define alpha, gamma, and lambda
alpha = 0.2;
gamma = 0.8;
lambda = 604.4;
% run FSGL fit
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit(...
    X,...
    Y,...
    Kmatrix,...
    nj,...
    nd,...
    ngroups,...
    groupsizes',...
    alpha,...
    gamma,...
    lambda,...
    1000,... % rho (initial ADMM step size)
    6000);   % max number of ADMM iterations
% toc
% save estimated beta
beta_0_2_0_8 = beta;
% save predicted y for test data
pred_y_test_0_2_0_8 = Xtest * beta;
% save output 
save('alpha_0_2_gamma_0_8_lambda_604.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% fused group lasso
% define alpha, gamma, and lambda
alpha = 0.0;
gamma = 0.8;
lambda = 604.4;
% run FSGL fit
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit(...
    X,...
    Y,...
    Kmatrix,...
    nj,...
    nd,...
    ngroups,...
    groupsizes',...
    alpha,...
    gamma,...
    lambda,...
    1000,... % rho (initial ADMM step size)
    6000);   % max number of ADMM iterations
% toc
% save estimated beta
beta_0_0_0_8 = beta;
% save predicted y for test data
pred_y_test_0_0_0_8 = Xtest * beta;
% save output 
save('alpha_0_0_gamma_0_8_lambda_604.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% run FSGL fit ADAPTIVE penalties
% define alpha, gamma, and lambda
%% adaptive lasso
alpha = 1.0;
gamma = 1.0;
lambda = 814.1;
% fit estimator
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit_adaptive(...
    X,...
    Y,...
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
    100,...  % rho (initial ADMM step size)
    10000);  % max number of ADMM iterations
% toc
% save estimated beta
beta_1_0_1_0_ada = beta;
% save predicted y for test data
pred_y_test_1_0_1_0_ada = Xtest * beta;
% save output 
save('alpha_1_0_gamma_1_0_lambda_814_adaptive.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% adaptive sparse group lasso
% define alpha, gamma, and lambda
alpha = 0.2;
gamma = 1.0;
lambda = 4041.4;
% fit estimator
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit_adaptive(...
    X,...
    Y,...
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
    100,...  % rho (initial ADMM step size)
    10000);  % max number of ADMM iterations
% toc
% save estimated beta
beta_0_2_1_0_ada = beta;
% save predicted y for test data
pred_y_test_0_2_1_0_ada = Xtest * beta;
% save output 
save('alpha_0_2_gamma_1_0_lambda_4041_adaptive.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% adaptive fused sparse group lasso
% define alpha, gamma, and lambda
alpha = 0.2;
gamma = 0.8;
lambda = 1096.6;
% fit estimator
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit_adaptive(...
    X,...
    Y,...
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
    100,...  % rho (initial ADMM step size)
    10000);  % max number of ADMM iterations
% toc
% save estimated beta
beta_0_2_0_8_ada = beta;
% save predicted y for test data
pred_y_test_0_2_0_8_ada = Xtest * beta;
% save output 
save('alpha_0_2_gamma_0_8_lambda_1097_adaptive.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% adaptive fused group lasso
% define alpha, gamma, and lambda
alpha = 0.0;
gamma = 0.8;
lambda = 1423.5;
% fit estimator
% tic
[pred_y, mse_y, beta, theta, mu, niter] = fsglfit_adaptive(...
    X,...
    Y,...
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
    100,...  % rho (initial ADMM step size)
    10000);  % max number of ADMM iterations
% toc
% save estimated beta
beta_0_0_0_8_ada = beta;
% save predicted y for test data
pred_y_test_0_0_0_8_ada = Xtest * beta;
% save output 
save('alpha_0_0_gamma_0_8_lambda_1424_adaptive.mat', 'alpha', 'gamma',...
    'lambda', 'pred_y', 'mse_y', 'beta', 'theta', 'mu', 'niter')


%% save estimated betas and initialbeta (ridge) as csv
% columns of beta_hat in order:
% ridge estimates 
% alpha 1.0 gamma 1.0
% alpha 0.2 gamma 1.0
% alpha 0.2 gamma 0.8
% alpha 0.0 gamma 0.8
% alpha 1.0 gamma 1.0 adaptive
% alpha 0.2 gamma 1.0 adaptive
% alpha 0.2 gamma 0.8 adaptive
% alpha 0.0 gamma 0.8 adaptive
beta_hat = [...
    initialbeta...
    beta_1_0_1_0...
    beta_0_2_1_0...
    beta_0_2_0_8...
    beta_0_0_0_8...
    beta_1_0_1_0_ada...
    beta_0_2_1_0_ada...
    beta_0_2_0_8_ada...
    beta_0_0_0_8_ada];
dlmwrite('beta_hat.csv', beta_hat, 'delimiter', ',', 'precision', 8);

%% save training y and predicted y as csv
% columns of y_hat in order:
% actual Y
% ridge estimates 
% alpha 1.0 gamma 1.0
% alpha 0.2 gamma 1.0
% alpha 0.2 gamma 0.8
% alpha 0.0 gamma 0.8
% alpha 1.0 gamma 1.0 adaptive
% alpha 0.2 gamma 1.0 adaptive
% alpha 0.2 gamma 0.8 adaptive
% alpha 0.0 gamma 0.8 adaptive
pred_y_test_ridge = Xtest * initialbeta;
yhat = [...
    Ytest...
    pred_y_test_ridge...
    pred_y_test_1_0_1_0...
    pred_y_test_0_2_1_0...
    pred_y_test_0_2_0_8...
    pred_y_test_0_0_0_8...
    pred_y_test_1_0_1_0_ada...
    pred_y_test_0_2_1_0_ada...
    pred_y_test_0_2_0_8_ada...
    pred_y_test_0_0_0_8_ada];
dlmwrite('y_hat.csv', yhat, 'delimiter', ',', 'precision', 8);
