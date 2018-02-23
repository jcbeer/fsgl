#' Fused Sparse Group Lasso
#'
#' @description Fits a regression model for one set of tuning parameter values by minimizing a penalized least squares loss function.
#' @param X a n*p matrix of predictor variables with observations in rows.
#' @param Y a n*1 vector of the response variable.
#' @param K a (nj + nd + ng)*p matrix encoding the lasso penalty (first nj rows), the graph structure for the fused penalty (the next nd rows), and the group structure for the group penalty (last ng rows). Can be made with function \code{makeKmatrix}.
#' @param nj number of rows of K that encode the lasso penalty. If lasso penalty is applied to all coefficients then this will equal p.
#' @param nd number of rows of K that encode the graph structure for the fused penalty.
#' @param ngroups number of groups for the group penalty.
#' @param groupsizes a vector of length ngroups that gives the size of each group in the order they appear in the K matrix. Sum should equal ng.
#' @param alpha tuning parameter that controls the degree of group (alpha = 0) vs L1 (alpha=1) sparsity. 0 <= alpha <= 1
#' @param gamma tuning parameter that controls the degree of sparsity (gamma=1) vs fusion (gamma=0) penalty. 0 <= gamma <= 1
#' @param lambda tuning parameter that controls the overall degree of regularization. 
#' @param beta0 starting values for beta. Defaults to zero vector.
#' @param theta0 starting values for theta. Defaults to zero vector.
#' @param mu0 starting values for mu. Defaults to zero vector.
#' @param beta_update_factor by default function will set this equal to \code{solve(t(Xcentered) \%*\% Xcentered + rho * t(K) \%*\% K)}. It may take a long time to compute and does not depend on tuning parameters, so it may be provided as an argument if fitting the model at many sets of tuning parameters. 
#' @param beta_update_term by default function will set this equal toset equal to \code{t(Xcentered) \%*\% Ycentered}. Does not depend on tuning parameters, so it may be provided as an argument if fitting the model at many sets of tuning parameters. 
#' @param rho step size for ADMM algorithm. Defaults to 1.
#' @param Niter number of ADMM iterations. Defaults to 2000.
#' @param epsilon_abs absolute tolerance for ADMM convergence. Defaults to 10^-3.
#' @param epsilon_rel relative tolerance for ADMM convergence. Defaults to 10^-3.
#' @param verbose if TRUE will print number of ADMM iterations used. Defaults to FALSE.
#' @keywords lasso
#' @export
#' @examples 
#' Kmatrix <- makeKmatrix2d(dim1=3, dim2=3, groups=c(1, 2, 3, 2, 2, 1, 3, 3, 1))
#' x <- matrix(rnorm(54), nrow=6)
#' y <- x %*% c(0, 0, 10, 0, 0, 0, 10, 10, 0)
#' fsgl.fit(X=x, Y=y, K=Kmatrix$K, nj=Kmatrix$nj, nd=Kmatrix$nd, ngroups=Kmatrix$ngroups, groupsizes=Kmatrix$groupsizes, alpha=0.1, gamma=0.8, lambda=5)


fsgl.fit <-
  function(X,
           Y,
           K,
           nj,
           nd,
           ngroups,
           groupsizes,
           alpha,
           gamma,
           lambda,
           beta0 = NULL,
           theta0 = NULL,
           mu0 = NULL,
           beta_update_factor = NULL,
           beta_update_term1 = NULL,
           rho = 1,
           Niter = 2000,
           epsilon_abs = 10 ^ -3,
           epsilon_rel = 10 ^ -3, 
           verbose = FALSE) {
    # center the columns of X and y
    # so we won't have to estimate an intercept
    Xcentered <- scale(X, center = TRUE, scale = FALSE)
    Ycentered <- scale(Y, center = TRUE, scale = FALSE)
    # make a list of 'group' indices for theta update step
    # where 'group' includes all lasso, fused lasso, and group lasso terms
    group_indices <- list()
    ni <- c(rep(1, nj), rep(1, nd), groupsizes) # vector of group sizes
    cumsumni <- cumsum(ni)
    N <- nj + nd + ngroups
    for (i in 1:N){
      group_indices[[i]] <- (max(0, cumsumni[i-1]) + 1):cumsumni[i]
    }
    # if not supplied as arguments then
    # calculate beta update terms
    # which don't change over loop iterations
    if (is.null(beta_update_factor))
      beta_update_factor <-
      solve(t(Xcentered) %*% Xcentered + rho * t(K) %*% K)
    if (is.null(beta_update_term1))
      beta_update_term1 <- t(Xcentered) %*% Ycentered
    # set initial values to zero if not specified
    if (is.null(beta0))
      beta0 <- rep(0, ncol(K))
    if (is.null(theta0))
      theta0 <- rep(0, nrow(K))
    if (is.null(mu0))
      mu0 <- rep(0, nrow(K))
    # initialize parameters
    current_beta <- beta0
    current_theta <- theta0
    current_mu <- mu0
    # define lambda vector
    Lambda <- c(rep(alpha * gamma * lambda, nj),
                rep((1 - gamma) * lambda, nd),
                rep((1 - alpha) * gamma * lambda, ngroups))
    # LOOP over ADMM iterations
    for (t in 1:Niter) {
      # UPDATE PARAMETERS
      # update beta
      beta_update_term2 <-
        t(K) %*% matrix((current_mu - rho * current_theta), ncol = 1)
      current_beta <-
        beta_update_factor %*% (beta_update_term1 - beta_update_term2)
      # update theta
      eta <- K %*% current_beta + current_mu / rho
      previous_theta <- current_theta
      for (i in 1:N) {
        current_theta[group_indices[[i]]] <-
          softthresh(eta[group_indices[[i]]], kappa = (Lambda[i] * sqrt(length(eta[group_indices[[i]]]))) /
                       rho)
      }
      # update mu
      current_mu <- current_mu + rho * (K %*% current_beta - current_theta)
      # CHECK STOPPING CONDITIONS
      # calculate dual residual
      epsilon_dual <-
        sqrt(length(current_theta)) * epsilon_abs + epsilon_rel * sqrt(sum((t(K) %*% current_mu) ^ 2))
      current_s <-
        sqrt(sum((rho * t(K) %*% (previous_theta - current_theta)) ^ 2))
      # calculate primal residual
      epsilon_primal <-
        sqrt(length(current_beta)) * epsilon_abs + epsilon_rel * max(sqrt(sum((K %*% current_beta) ^ 2)), sqrt(sum(current_theta ^ 2)))
      current_r <- sqrt(sum((K %*% current_beta - current_theta) ^ 2))
      # break loop if conditions are met
      if (current_s < epsilon_dual &
          current_r < epsilon_primal)
        break
    } # end loop over ADMM iterations
    if (verbose == TRUE) print(t)
    # save outputs
    pred.y <- Xcentered %*% current_beta
    mse.y <- mean((Ycentered - pred.y) ^ 2)
    output <-
      list(
        pred.y = pred.y,
        mse.y = mse.y,
        beta = current_beta,
        theta = current_theta,
        mu = current_mu,
        niter = t
      )
    return(output)
  }
