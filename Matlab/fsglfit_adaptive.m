%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FSGLFIT_ADAPTIVE
% Function to fit fused sparse group lasso
% for adaptive penalties.
% Uses output from function makeKmatrix_adaptive.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Methodology described in the manuscript
% "Incorporating Prior Information with Fused Sparse Group Lasso:
% Application to Prediction of Clinical Measures from Neuroimages"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
% X: n*p predictor matrix
% Y: n*1 scalar outcome 
% K: K matrix 
%   (output from makeKmatrix_adaptive function)
% nj: number of predictors, p 
%   (output from makeKmatrix_adaptive function)
% nd: number of rows of D matrix 
%   (output from makeKmatrix_adaptive function)
% ngroups: number of coefficient groups 
%   (output from makeKmatrix_adaptive function) 
% groupsizes: vector of coefficient group sizes 
%   (output from makeKmatrix_adaptive function) 
% weights: adaptive penalty weights to rescale lambda  
%   (output from makeKmatrix_adaptive function)
% initialbeta: initial coefficient values 
% alpha: tuning parameter that controls balance between group (alpha = 0)
%   and L1 lasso (alpha = 1) penalty terms 
% gamma: tuning parameter that controls balance between fusion (gamma = 0)
%   and sparsity (gamma = 1) penalty terms 
% lambda: tuning parameter that controls overall level of regularization
% rho: initial ADMM step-size
% Niter: maximum number of ADMM iterations
%%% OUTPUTS
% pred_y: predicted Y values
% mse_y: mean squared error of the predicted Y values
% beta: final value for beta (estimated beta)
% theta: final value for theta
% mu: final value for mu
% niter: number of iterations before stopping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pred_y, mse_y, beta, theta, mu, niter]...
    = fsglfit(...
    X,...
    Y,...
    K,...
    nj,... 
    nd,... 
    ngroups,... 
    groupsizes,... 
    weights,...
    initialbeta,...
    alpha,... 
    gamma,...
    lambda,...
    rho,...
    Niter)
    % define soft threshold function
    function b = softthresh(a, kappa)
        if sum(a(:)==0) == size(a(:))
            b = a;
        else
            b = max(0, (1 - kappa/norm(a))) * a;
        end
    end
    % define epsilons
    % NOTE: These can be changed to match desired tolerance
    epsilon_abs = 10 ^ -3; 
    epsilon_rel = 10 ^ -3;
    % center columns of X and Y
    Xcentered = detrend(X, 'constant');
    Ycentered = detrend(Y, 'constant');
    % make a cell array with indices of each group for theta update
    ni = [repelem(1, nj) repelem(1, nd) groupsizes]; % vector of group sizes
    cumsumni = [0 cumsum(ni)];
    N = nj + nd + ngroups;
    groupindices = cell(1, N);
    for i = 2:(N + 1)
        groupindices{i-1} = (cumsumni(i-1) + 1):cumsumni(i);
    end
    % calculate beta update factor 
    % which does not change over loop iterations
    beta_update_factor = Xcentered' * Xcentered + rho * (K' * K);
    % calculate beta update term 
    % which does not change over loop iterations
    beta_update_term1 = Xcentered' * Ycentered;
    % set initial values
    beta0 = initialbeta;
    theta0 = repelem(0, size(K, 1));
    mu0 = repelem(0, size(K, 1));
    % initialize parameters
    currentbeta = beta0';
    currenttheta = theta0';
    currentmu = mu0';
    clear beta0 theta0 mu0
    % define lambda vector
    Lambda = [...
        repelem(alpha * gamma * lambda, nj)...
        repelem((1 - gamma) * lambda, nd)...
        repelem((1 - alpha) * gamma * lambda, ngroups)];
    % LOOP over ADMM iterations
    for t = 1:Niter
        % UPDATE PARAMETERS
        % update beta
        beta_update_term2 = K' * (currentmu - rho * currenttheta);
        currentbeta = beta_update_factor \ (beta_update_term1 - beta_update_term2);
        % update theta
        eta = K * currentbeta + currentmu / rho;
        previoustheta = currenttheta;
        for i = 1:N
            currenttheta(groupindices{i}) = softthresh(eta(groupindices{i}),...
                (Lambda(i) * weights(i)) / rho);
        end
        % update mu
        currentmu = currentmu + rho * (K * currentbeta - currenttheta);
        % CHECK STOPPING CONDITIONS
        % calculate dual residual
        epsilon_dual = sqrt(size(currenttheta, 1)) * epsilon_abs...
            + epsilon_rel * norm(K' * currentmu);
        current_s = norm(rho * K' * (previoustheta - currenttheta));
        % calculate primal residual
        epsilon_primal = sqrt(size(currentbeta, 1)) * epsilon_abs... 
            + epsilon_rel * max(norm(K * currentbeta), norm(currenttheta));
        current_r = norm(K * currentbeta - currenttheta);
        % break loop if stopping conditions are met
        if (current_s < epsilon_dual) && (current_r < epsilon_primal)
            break
        end
        % update rho 
        if (current_r > 10 * current_s)
            rho = 2 * rho;
        elseif (current_s > 10 * current_r)
            rho = rho / 2;
        else rho = rho;
        end
    end % end loop over ADMM iterations
    % save outputs
    pred_y = Xcentered * currentbeta;
    mse_y = mean((Ycentered - pred_y).^2);
    beta = currentbeta;
    theta = currenttheta;
    mu = currentmu;
    niter = t;
end
