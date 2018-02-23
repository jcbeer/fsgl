%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused Sparse Group Lasso ABIDE Application
% Plot Predicted vs. Actual Adjusted SRS for test data
% for adaptive penalties (Figure 6B)
% Calculate statistics reported in Table 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script used for analyses reported in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%% INPUTS: 
% beta_hat.csv
% y_hat.csv
% phenotest.txt
% srs_train.txt
% trainXstd_5476.txt
% srs_test.txt
% testXstd_5476.txt
%%% OUTPUTS:
% pdf figure, predicted vs. actual adjusted SRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data
beta = csvread('beta_hat.csv');
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
yhat = csvread('y_hat.csv');
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
phenotest = readtable('phenotest.txt');
dxgp = table2array(phenotest(:,3));
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



%% plot predicted Y test versus actual Y test
% COLOR BY DX GROUP
f = figure
% alpha = 1, gamma = 1
subplot(2, 2, 1)
plot(yhat(:,1), yhat(:,7), 'o', 'MarkerSize', 1, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
axis([-100 100 -50 50]);
pbaspect([1 1 1])
l = lsline
set(l, 'LineWidth', 2, 'Color', 'k');
title({'Adaptive'; '\alpha = 1.0, \gamma = 1.0, \lambda = 814'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
text(-20,-40,'r = 0.41, p = 0.006', 'FontSize', 12)
hold on 
p1 = plot(yhat(dxgp==1,1), yhat(dxgp==1,7), 'o', 'MarkerSize', 5, 'LineWidth', 2)
p2 = plot(yhat(dxgp==2,1), yhat(dxgp==2,7), 'x', 'MarkerSize', 7, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'LineWidth', 2)
legend([p1 p2], 'ASD', 'TD', 'Location','northwest')
% alpha = 0.2, gamma = 1
subplot(2, 2, 2)
plot(yhat(:,1), yhat(:,8), 'o', 'MarkerSize', 1, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
axis([-100 100 -50 50]);
pbaspect([1 1 1])
l = lsline
set(l, 'LineWidth', 2, 'Color', 'k');
title({'Adaptive'; '\alpha = 0.2, \gamma = 1.0, \lambda = 4041'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
text(-20,-40,'r = 0.40, p = 0.008', 'FontSize', 12)
hold on 
plot(yhat(dxgp==1,1), yhat(dxgp==1,8), 'o', 'MarkerSize', 5, 'LineWidth', 2)
plot(yhat(dxgp==2,1), yhat(dxgp==2,8), 'x', 'MarkerSize', 7, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'LineWidth', 2)
% alpha = 0.2, gamma = 0.8
subplot(2, 2, 3)
plot(yhat(:,1), yhat(:,9), 'o', 'MarkerSize', 1, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
axis([-100 100 -50 50]);
pbaspect([1 1 1])
l = lsline
set(l, 'LineWidth', 2, 'Color', 'k');
title({'Adaptive'; '\alpha = 0.2, \gamma = 0.8, \lambda = 1097'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
text(-20,-40,'r = 0.44, p = 0.003', 'FontSize', 12)
hold on 
plot(yhat(dxgp==1,1), yhat(dxgp==1,9), 'o', 'MarkerSize', 5, 'LineWidth', 2)
plot(yhat(dxgp==2,1), yhat(dxgp==2,9), 'x', 'MarkerSize', 7, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'LineWidth', 2)
% alpha = 0, gamma = 0.8
subplot(2, 2, 4)
plot(yhat(:,1), yhat(:,10), 'o', 'MarkerSize', 1, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w')
axis([-100 100 -50 50]);
pbaspect([1 1 1])
l = lsline
set(l, 'LineWidth', 2, 'Color', 'k');
title({'Adaptive'; '\alpha = 0.0, \gamma = 0.8, \lambda = 1424'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
text(-20,-40,'r = 0.39, p = 0.008', 'FontSize', 12)
hold on 
plot(yhat(dxgp==1,1), yhat(dxgp==1,10), 'o', 'MarkerSize', 5, 'LineWidth', 2)
plot(yhat(dxgp==2,1), yhat(dxgp==2,10), 'x', 'MarkerSize', 7, 'MarkerFaceColor', [0, 0.4470, 0.7410], 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'LineWidth', 2)
% add axis title
bkclr = get(f,'Color')  % The figure's background color
bkax = axes; % this creates a figure-filling new axes
set(bkax,'Xcolor',[1 1 1],'Ycolor', [1 1 1],'box','off','Color','none') % make it 'disappear'
xlab = xlabel('Actual Adjusted Social Responsiveness Scale Score','Color','k', 'FontWeight','bold')
ylab = ylabel('Predicted Adjusted Social Responsiveness Scale Score','Color','k', 'FontWeight','bold')
set(xlab, 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
set(ylab, 'Units', 'Normalized', 'Position', [-0.08, 0.5, 0]);
uistack(bkax,'bottom'); % this moves it to the background.
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)

% save the figure as pdf
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'pred_vs_actual.pdf','-dpdf','-r0')


%% calculate MSE and correlation of actual vs. predicted Y
%% alpha = 1.0, gamma = 1.0 (lasso)
% calculate y hat
yhat_train_1_0_1_0 = X*beta(:,2);
yhat_test_1_0_1_0 = Xtest*beta(:,2);
% MSE and corr training
sum((yhat_train_1_0_1_0 - Y).^2)/175
[r, p] = corr(Y, yhat_train_1_0_1_0)
% MSE and corr test
sum((yhat_test_1_0_1_0 - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_1_0_1_0)

%% alpha = 0.2, gamma = 1.0 (sparse group lasso)
% calculate y hat
yhat_train_0_2_1_0 = X*beta(:,3);
yhat_test_0_2_1_0 = Xtest*beta(:,3);
% MSE and corr training
sum((yhat_train_0_2_1_0 - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_2_1_0)
% MSE and corr test
sum((yhat_test_0_2_1_0 - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_2_1_0)

%% alpha = 0.2, gamma = 0.8 (fused sparse group lasso)
% calculate y hat
yhat_train_0_2_0_8 = X*beta(:,4);
yhat_test_0_2_0_8 = Xtest*beta(:,4);
% MSE and corr training
sum((yhat_train_0_2_0_8 - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_2_0_8)
% MSE and corr test
sum((yhat_test_0_2_0_8 - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_2_0_8)

%% alpha = 0.0, gamma = 0.8 (fused group lasso)
% calculate y hat
yhat_train_0_0_0_8 = X*beta(:,5);
yhat_test_0_0_0_8 = Xtest*beta(:,5);
% MSE and corr training
sum((yhat_train_0_0_0_8 - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_0_0_8)
% MSE and corr test
sum((yhat_test_0_0_0_8 - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_0_0_8)

%% alpha = 1.0, gamma = 1.0 (adaptive lasso)
% calculate y hat
yhat_train_1_0_1_0_ada = X*beta(:,6);
yhat_test_1_0_1_0_ada = Xtest*beta(:,6);
% MSE and corr training
sum((yhat_train_1_0_1_0_ada - Y).^2)/175
[r, p] = corr(Y, yhat_train_1_0_1_0_ada)
% MSE and corr test
sum((yhat_test_1_0_1_0_ada - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_1_0_1_0_ada)

%% alpha = 0.2, gamma = 1.0 (adaptive sparse group lasso)
% calculate y hat
yhat_train_0_2_1_0_ada = X*beta(:,7);
yhat_test_0_2_1_0_ada = Xtest*beta(:,7);
% MSE and corr training
sum((yhat_train_0_2_1_0_ada - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_2_1_0_ada)
% MSE and corr test
sum((yhat_test_0_2_1_0_ada - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_2_1_0_ada)

%% alpha = 0.2, gamma = 0.8 (adaptive fused sparse group lasso)
% calculate y hat
yhat_train_0_2_0_8_ada = X*beta(:,8);
yhat_test_0_2_0_8_ada = Xtest*beta(:,8);
% MSE and corr training
sum((yhat_train_0_2_0_8_ada - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_2_0_8_ada)
% MSE and corr test
sum((yhat_test_0_2_0_8_ada - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_2_0_8_ada)

%% alpha = 0.0, gamma = 0.8 (adaptive fused group lasso)
% calculate y hat
yhat_train_0_0_0_8_ada = X*beta(:,9);
yhat_test_0_0_0_8_ada = Xtest*beta(:,9);
% MSE and corr training
sum((yhat_train_0_0_0_8_ada - Y).^2)/175
[r, p] = corr(Y, yhat_train_0_0_0_8_ada)
% MSE and corr test
sum((yhat_test_0_0_0_8_ada - Ytest).^2)/44
[r, p] = corr(Ytest, yhat_test_0_0_0_8_ada)

