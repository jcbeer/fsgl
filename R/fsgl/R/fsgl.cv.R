#' Cross-validation for Fused Sparse Group Lasso
#'
#' @description Does k-fold validation for \code{fsgl.fit}.
#' @param X a n*p matrix of predictor variables with observations in rows.
#' @param Y a n*1 vector of the response variable.
#' @param K a (nj + nd + ng)*p matrix encoding the lasso penalty (first nj rows), the graph structure for the fused penalty (the next nd rows), and the group structure for the group penalty (last ng rows). Can be made with function \code{makeKmatrix}.
#' @param nj number of rows of K that encode the lasso penalty. If lasso penalty is applied to all coefficients then this will equal p.
#' @param nd number of rows of K that encode the graph structure for the fused penalty.
#' @param ngroups number of groups for the group penalty.
#' @param groupsizes a vector of length ngroups that gives the size of each group in the order they appear in the K matrix. Sum should equal ng.
#' @param alphagamma a two-column matrix with rows containing alpha, gamma pairs at which cross-validation will be done for a range of lambda values. alpha is a tuning parameter that controls the degree of group (alpha = 0) vs L1 (alpha=1) sparsity. gamma is a tuning parameter that controls the degree of sparsity (gamma=1) vs fusion (gamma=0) penalty. 0 <= alpha, gamma <= 1
#' @param lambda vector of decreasing lambda values at which k-fold cross-validation will be done. tuning parameter that controls the overall degree of regularization. 
#' @param k number of folds for cross-validation.
#' @param folds a vector encoding the fold assignments. By default folds will be randomly assigned.
#' @param beta0 starting values for beta. Defaults to zero vector.
#' @param theta0 starting values for theta. Defaults to zero vector.
#' @param mu0 starting values for mu. Defaults to zero vector.
#' @param rho step size for ADMM algorithm. Defaults to 1.
#' @param Niter number of ADMM iterations used for each fsgl fit. Defaults to 2000.
#' @param epsilon_abs absolute tolerance for ADMM convergence. Defaults to 10^-3.
#' @param epsilon_rel relative tolerance for ADMM convergence. Defaults to 10^-3.
#' @param verbose if TRUE will print progress through alpha, gamma combinations for each fold. Defaults to FALSE.
#' @keywords lasso
#' @export
#' @examples 
#' Kmatrix <- makeKmatrix2d(dim1=3, dim2=3, groups=c(1, 2, 3, 2, 2, 1, 3, 3, 1))
#' x <- matrix(rnorm(54), nrow=6)
#' y <- x %*% c(0, 0, 10, 0, 0, 0, 10, 10, 0)
#' # define lambda
#' Lambda <- 10 ^ seq(3, -3, length = 50)
#' # define alpha and gamma 
#' alphas <- c(rep(0,5), rep(0.2, 4), rep(0.5, 4), rep(0.8, 4), rep(1, 4))
#' gammas <- c(0, rep(c(0.2, 0.5, 0.8, 1), 5))
#' AlphaGamma <- cbind(alphas, gammas)
#' cv <- fsgl.cv(X=x, Y=y, K=Kmatrix$K, nj=Kmatrix$nj, nd=Kmatrix$nd, ngroups=Kmatrix$ngroups, groupsizes=Kmatrix$groupsizes, alphagamma=AlphaGamma, lambda=Lambda, k=3, verbose=TRUE)


fsgl.cv <-
  function(X,
           Y,
           K,
           nj,
           nd,
           ngroups,
           groupsizes,
           alphagamma,
           lambda,
           k,
           folds = NULL,
           beta0 = NULL,
           theta0 = NULL,
           mu0 = NULL,
           rho = 1,
           Niter = 2000,
           epsilon_abs = 10 ^ -3,
           epsilon_rel = 10 ^ -3, 
           verbose = FALSE) {
    if (is.null(folds)) folds <- sample(rep_len(1:k, nrow(X)))
    # create array to save results
    results <- array(dim=c(nrow(alphagamma), length(lambda), k))
    # LOOP over folds
    for (fold in 1:k){
      if (verbose == TRUE) print(paste0('Starting Fold ', fold))
      # select in and out of CV fold samples and center the data
      curX <- scale(X[folds != fold,], center=TRUE, scale=FALSE)
      curY <- scale(Y[folds != fold], center=TRUE, scale=FALSE)
      curXout <- scale(X[folds == fold,], center=TRUE, scale=FALSE)
      curYout <- scale(Y[folds == fold], center=TRUE, scale=FALSE)
      # calculate beta update terms for CV sample
      beta_update_factor_cv <- solve(t(curX) %*% curX + rho * t(K) %*% K)
      beta_update_term1_cv <- t(curX) %*% curY
      # LOOP over alpha, gamma pairs
      for (i in 1:nrow(alphagamma)){
        # set alpha and gamma 
        alpha <- alphagamma[i, 1]
        gamma <- alphagamma[i, 2]
        if (verbose == TRUE) print(paste0('alpha = ', alpha, ', gamma = ', gamma))
        # initialize parameters
        # start all parameters at zero when lambda is largest
        if (is.null(beta0))
          beta0 <- rep(0, ncol(K))
        if (is.null(theta0))
          theta0 <- rep(0, nrow(K))
        if (is.null(mu0))
          mu0 <- rep(0, nrow(K))
        current_beta <- beta0
        current_theta <- theta0
        current_mu <- mu0
        # LOOP over lambda values
        for (j in 1:length(lambda)){
          # fit estimator to the current CV sample
          fit <-
            fsgl.fit(
              X=curX,
              Y=curY,
              K=K,
              nj=nj,
              nd=nd,
              ngroups=ngroups,
              groupsizes=groupsizes,
              alpha=alpha,
              gamma=gamma,
              lambda=lambda[j],
              beta0 = current_beta,
              theta0 = current_theta,
              mu0 = current_mu,
              beta_update_factor = beta_update_factor_cv,
              beta_update_term1 = beta_update_term1_cv,
              rho = 1,
              Niter = 2000,
              epsilon_abs = 10 ^ -3,
              epsilon_rel = 10 ^ -3,
              verbose = TRUE
            )
          # calculate and save the out of sample MSE
          pred.y.cv <- curXout %*% fit$beta
          results[i, j, fold] <- mean((curYout - pred.y.cv) ^ 2)
          # update the parameters for warm start at next lambda value
          current_beta <- fit$beta
          current_theta <- fit$theta
          current_mu <- fit$mu
        } # end loop over lambda values
      } # end loop over alpha, gamma pairs
    } # end loop over folds
    # determine optimal lambda for each alpha, gamma pair
    # calculate mean and sd cve across folds
    mean.cve <- apply(results, c(1, 2), mean)
    sd.cve <- apply(results, c(1, 2), sd)
    # find which lambda is minimum for each alpha, gamma pair
    min.cve.index <- apply(mean.cve, 1, which.min)
    opt.lambda <- lambda[min.cve.index]
    opt.mean.cve <- apply(mean.cve, 1, function(x) x[which.min(x)])
    opt.sd.cve <- sd.cve[cbind(1:nrow(sd.cve), min.cve.index)]
    optimal.lambdas <- cbind(alphagamma, opt.lambda, opt.mean.cve, opt.sd.cve)
    return(list(optimal.lambdas=optimal.lambdas, mean.cve=mean.cve, sd.cve=sd.cve))
  }
