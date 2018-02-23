%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Analyze cross-validation results
% Plot CV error curves (Figure 6A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat
% cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid2.mat
% cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid3.mat
% cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid1.mat
% cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid2.mat
% cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid3.mat
% cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid1.mat
% cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid2.mat
% cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid3.mat
% cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid1.mat
% cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid2.mat
% cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid3.mat
% cvresults_alpha_0_0_gamma_0_8_lambdagrid1.mat
% cvresults_alpha_0_0_gamma_0_8_lambdagrid2.mat
% cvresults_alpha_0_0_gamma_0_8_lambdagrid3.mat
% cvresults_alpha_0_2_gamma_0_8_lambdagrid1.mat
% cvresults_alpha_0_2_gamma_0_8_lambdagrid2.mat
% cvresults_alpha_0_2_gamma_0_8_lambdagrid3.mat
% cvresults_alpha_1_0_gamma_1_0_lambdagrid1.mat
% cvresults_alpha_1_0_gamma_1_0_lambdagrid2.mat
% cvresults_alpha_1_0_gamma_1_0_lambdagrid3.mat
% cvresults_alpha_0_2_gamma_1_0_lambdagrid1.mat
% cvresults_alpha_0_2_gamma_1_0_lambdagrid2.mat
% cvresults_alpha_0_2_gamma_1_0_lambdagrid3.mat
%%% OUTPUTS:
% pdf figure: cross-validation curves
% optimal lambda and CVMSE reported in Table 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data for adaptive penalty
cd('./ABIDE/data/')
% load each result file and store the mean and SD
% cross validation error
%% alpha 0.0, gamma 0.8, adaptive
load('cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid1.mat')
row1 = mean(cvresults, 3);
row1sd = std(cvresults, [], 3);
load('cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid2.mat')
row1 = [row1 mean(cvresults, 3)];
row1sd = [row1sd std(cvresults, [], 3)];
load('cvresults_adaptive_alpha_0_0_gamma_0_8_lambdagrid3.mat')
row1 = [row1 mean(cvresults, 3)];
row1sd = [row1sd std(cvresults, [], 3)];
%% alpha 0.2, gamma 0.8, adaptive
load('cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid1.mat')
row2 = mean(cvresults, 3);
row2sd = std(cvresults, [], 3);
load('cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid2.mat')
row2 = [row2 mean(cvresults, 3)];
row2sd = [row2sd std(cvresults, [], 3)];
load('cvresults_adaptive_alpha_0_2_gamma_0_8_lambdagrid3.mat')
row2 = [row2 mean(cvresults, 3)];
row2sd = [row2sd std(cvresults, [], 3)];
%% alpha 1.0, gamma 1.0, adaptive
load('cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid1.mat')
row3 = mean(cvresults, 3);
row3sd = std(cvresults, [], 3);
load('cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid2.mat')
row3 = [row3 mean(cvresults, 3)];
row3sd = [row3sd std(cvresults, [], 3)];
load('cvresults_adaptive_alpha_1_0_gamma_1_0_lambdagrid3.mat')
row3 = [row3 mean(cvresults, 3)];
row3sd = [row3sd std(cvresults, [], 3)];
%% alpha 0.2, gamma 1.0, adaptive
load('cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid1.mat')
row4 = mean(cvresults, 3);
row4sd = std(cvresults, [], 3);
load('cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid2.mat')
row4 = [row4 mean(cvresults, 3)];
row4sd = [row4sd std(cvresults, [], 3)];
load('cvresults_adaptive_alpha_0_2_gamma_1_0_lambdagrid3.mat')
row4 = [row4 mean(cvresults, 3)];
row4sd = [row4sd std(cvresults, [], 3)];
%% save in one matrix
% please note that this reorders the rows
cve = [row3; row4; row2; row1];
cve_sd = [row3sd; row4sd; row2sd; row1sd];
clear row1 row2 row3 row4 row1sd row2sd row3sd row4sd

%% load data for non-adaptive penalty
cd('./ABIDE/data/')
% load each result file and store the mean and SD
% cross validation error
%% alpha 0.0, gamma 0.8, non-adaptive
load('cvresults_alpha_0_0_gamma_0_8_lambdagrid1.mat')
row1 = mean(cvresults, 3);
row1sd = std(cvresults, [], 3);
load('cvresults_alpha_0_0_gamma_0_8_lambdagrid2.mat')
row1 = [row1 mean(cvresults, 3)];
row1sd = [row1sd std(cvresults, [], 3)];
load('cvresults_alpha_0_0_gamma_0_8_lambdagrid3.mat')
row1 = [row1 mean(cvresults, 3)];
row1sd = [row1sd std(cvresults, [], 3)];
%% alpha 0.2, gamma 0.8, non-adaptive
load('cvresults_alpha_0_2_gamma_0_8_lambdagrid1.mat')
row2 = mean(cvresults, 3);
row2sd = std(cvresults, [], 3);
load('cvresults_alpha_0_2_gamma_0_8_lambdagrid2.mat')
row2 = [row2 mean(cvresults, 3)];
row2sd = [row2sd std(cvresults, [], 3)];
load('cvresults_alpha_0_2_gamma_0_8_lambdagrid3.mat')
row2 = [row2 mean(cvresults, 3)];
row2sd = [row2sd std(cvresults, [], 3)];
%% alpha 1.0, gamma 1.0, non-adaptive
load('cvresults_alpha_1_0_gamma_1_0_lambdagrid1.mat')
row3 = mean(cvresults, 3);
row3sd = std(cvresults, [], 3);
load('cvresults_alpha_1_0_gamma_1_0_lambdagrid2.mat')
row3 = [row3 mean(cvresults, 3)];
row3sd = [row3sd std(cvresults, [], 3)];
load('cvresults_alpha_1_0_gamma_1_0_lambdagrid3.mat')
row3 = [row3 mean(cvresults, 3)];
row3sd = [row3sd std(cvresults, [], 3)];
%% alpha 0.2, gamma 1.0, non-adaptive
load('cvresults_alpha_0_2_gamma_1_0_lambdagrid1.mat')
row4 = mean(cvresults, 3);
row4sd = std(cvresults, [], 3);
load('cvresults_alpha_0_2_gamma_1_0_lambdagrid2.mat')
row4 = [row4 mean(cvresults, 3)];
row4sd = [row4sd std(cvresults, [], 3)];
load('cvresults_alpha_0_2_gamma_1_0_lambdagrid3.mat')
row4 = [row4 mean(cvresults, 3)];
row4sd = [row4sd std(cvresults, [], 3)];
%% save in one matrix
% please note that this reorders the rows
cve = [row3; row4; row2; row1; cve];
cve_sd = [row3sd; row4sd; row2sd; row1sd; cve_sd];
clear row1 row2 row3 row4 row1sd row2sd row3sd row4sd

%% define lambda grid 1, 2, 3 
% lambda grid 1 and 2
exponent = 0:(7/47):7;
lambdagrid12 = (repelem(exp(1), 48) .^ exponent)';
% lambda grid 3
exponent = 7:(3/23):10;
lambdagrid3 = (repelem(exp(1), 24) .^ exponent)';
% put all together
lambdagrid = [lambdagrid12; lambdagrid3];
clear exponent lambdagrid12 lambdagrid3


%% plot CVE curves
% mean CV MSE -- lambda on log scale
h = figure, hold on
plot(log(lambdagrid), cve(1,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(2,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(3,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(4,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(5,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(6,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(7,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(8,:), 'LineWidth', 2, 'Color',[1, 0.0784, 0.5765])
lgd = legend(...
    {'\alpha = 1.0, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 0.8',...
    '\alpha = 0.0, \gamma = 0.8',...
    '\alpha = 1.0, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 0.8, adaptive',...
    '\alpha = 0.0, \gamma = 0.8, adaptive'},...
    'Location', 'Southwest');
lgd.FontSize = 14;
axis([0 10 1300 1800]);
xlabel('log(\lambda)', 'FontSize', 14, 'FontWeight','bold')
ylabel('Cross-validation MSE', 'FontSize', 14, 'FontWeight','bold')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
% save the figure as pdf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'cvecurves2','-dpdf','-r0')


%% Non adaptive curves only
% mean CV MSE -- lambda on log scale
figure, hold on
plot(log(lambdagrid), cve(1,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(2,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(3,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(4,:), 'LineWidth', 2)
lgd = legend(...
    {'\alpha = 1.0, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 0.8',...
    '\alpha = 0.0, \gamma = 0.8',...
    },...
    'Location', 'Southwest');
lgd.FontSize = 14;
axis([0 10 1660 1780]);
xlabel('log(\lambda)', 'FontSize', 14)
ylabel('Cross-validation MSE', 'FontSize', 14)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)

%% Adaptive curves only
% mean CV MSE -- lambda on log scale
figure, hold on
plot(log(lambdagrid), cve(5,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(6,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(7,:), 'LineWidth', 2)
plot(log(lambdagrid), cve(8,:), 'LineWidth', 2)
lgd = legend(...
    {'\alpha = 1.0, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 0.8, adaptive',...
    '\alpha = 0.0, \gamma = 0.8, adaptive'},...
    'Location', 'Southwest');
lgd.FontSize = 14;
xlabel('log(\lambda)', 'FontSize', 14)
ylabel('Cross-validation MSE', 'FontSize', 14)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)


%% plot standard error curves
% mean CV MSE -- lambda on log scale
h = figure, hold on
plot(log(lambdagrid), cve_sd(1,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(2,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(3,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(4,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(5,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(6,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(7,:), 'LineWidth', 2)
plot(log(lambdagrid), cve_sd(8,:), 'LineWidth', 2, 'Color',[1, 0.0784, 0.5765])
lgd = legend(...
    {'\alpha = 1.0, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 1.0',...
    '\alpha = 0.2, \gamma = 0.8',...
    '\alpha = 0.0, \gamma = 0.8',...
    '\alpha = 1.0, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 1.0, adaptive',...
    '\alpha = 0.2, \gamma = 0.8, adaptive',...
    '\alpha = 0.0, \gamma = 0.8, adaptive'},...
    'Location', 'Southwest');
lgd.FontSize = 14;
% axis([0 10 1300 1800]);
xlabel('log(\lambda)', 'FontSize', 14)
ylabel('Cross-validation MSE standard deviation', 'FontSize', 14)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)


%% check location of minima
% optimal lambda and CVMSE reported in Table 3
% alpha = 1.0, gamma = 1.0
[val,ind]=min(cve(1,:))
lambdagrid(ind)
% alpha = 0.2, gamma = 1.0
[val,ind]=min(cve(2,:))
lambdagrid(ind)
% alpha = 0.2, gamma = 0.8
[val,ind]=min(cve(3,:))
lambdagrid(ind)
% alpha = 0.0, gamma = 0.8
[val,ind]=min(cve(4,:))
lambdagrid(ind)
% alpha = 1.0, gamma = 1.0, adaptive
[val,ind]=min(cve(5,:))
lambdagrid(ind)
% alpha = 0.2, gamma = 1.0, adaptive
[val,ind]=min(cve(6,:))
lambdagrid(ind)
% alpha = 0.2, gamma = 0.8, adaptive
[val,ind]=min(cve(7,:))
lambdagrid(ind)
% alpha = 0.0, gamma = 0.8, adaptive
[val,ind]=min(cve(8,:))
lambdagrid(ind)








