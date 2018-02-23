####################################################################
# Fused Sparse Group Lasso Simulation Study
# Run Simulation
####################################################################
# Variations of this script were used to produce the
# simulation results reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# none
### OUTPUTS:
# csv file of results for one simulation scenario
# results include training set MSE Y, MSE beta, estimated beta
####################################################################

######################################################
# strategy for choosing lambdas
# by 5-fold cross validation:
# we fix alpha and gamma and choose optimal lambda
# lambda1 = alpha*gamma*lambda
# lambda2 = (1-gamma)*lambda
# lambda3 = (1-alpha)*gamma*lambda
# alpha = 0, 02, 0.5, 0.8, 1
# gamma = 0, 0.2, 0.5, 0.8, 1
# lambda = 10^x where x is grid of 50 values from 3 to -3
######################################################

######################################################
# LOAD PACKAGES
######################################################
# fields package for image.plot function
# (optional; used to plot coefficient group structures)
library(fields)
# Matrix package for creating block diagonal matrix
library(Matrix)

######################################################
# SET OUTPUT DIRECTORY
# change if needed
######################################################
outputdir <- paste0(getwd(), '/')

######################################################
# DEFINE FUNCTIONS
######################################################

########################################
# Soft-thresholding function
# a is a scalar or vector
########################################
softthresh <- function(a, kappa){
  if(sum(a == 0) == length(a)){
    return(a)
  } else {
    return(max(0, (1 - kappa/sqrt(sum(a^2))))*a)
  }
}

########################################
# Function to fit FSGLasso by ADMM algorithm
# for one lambda and its corresponding set 
# of lambda1, lambda2, lambda3 values
### INPUTS
# Xmatrix: n*p predictor matrix
# Yvector: n-vector of outputs
# Kmatrix: matrix with p columns encoding spatial and group structure of coefficients
# beta0: initial beta values
# theta0: initial theta values
# mu0: initial mu values
# beta_update_factor: equal to solve(t(Xmatrix) %*% Xmatrix + rho * t(Kmatrix) %*% Kmatrix)
# beta_update_term1: equal to t(Xmatrix) %*% Yvector
# trueBetas: p-vector of true coefficient values
# groupindex: list of indices of elements for each effective group
# lambda: N-vector of lambda1, lambda2, lambda3, values, one for each effective group
# rho: ADMM step size
# p: number of predictors 
# N: number of effective groups
# Niter: maximum number of ADMM iterations
# epsilon_abs: to set absolute tolerance
# epsilon_rel: to set relative tolerance
### OUTPUTS
# stats: MSE for coefficients and predicted Y
# beta: final value of theta beta at that lambda
# beta_matrix: value of beta at each ADMM iteration
# theta: final value of theta at that lambda
# mu: final value of mu at that lambda
########################################
fsgl.fit <- function(Xmatrix ,Yvector, Kmatrix, beta0, theta0, mu0, beta_update_factor, beta_update_term1, trueBetas, groupindex, lambda, rho=1, p, N, Niter=2000, epsilon_abs=10^-3, epsilon_rel=10^-3){
  # INITIALIZE PARAMETERS
  current_beta <- beta0
  current_theta <- theta0
  current_mu <- mu0
  # create matrix to save betas for each iteration
  beta_matrix <- matrix(nrow=Niter, ncol=p)
  # BEGIN LOOP OVER ITERATIONS
  for (t in 1:Niter){
    # store current beta
    beta_matrix[t,] <- current_beta
    # UPDATE PARAMETERS
    # update beta
    beta_update_term2 <- t(Kmatrix) %*% matrix((current_mu - rho*current_theta), ncol=1)
    current_beta <- beta_update_factor %*% (beta_update_term1 - beta_update_term2)
    # update theta
    eta <- Kmatrix %*% current_beta + current_mu/rho
    previous_theta <- current_theta
    for (j in 1:N){
      current_theta[group_indices[[j]]] <- softthresh(eta[group_indices[[j]]], kappa=(lambda[j]*sqrt(length(eta[group_indices[[j]]])))/rho)
    }
    # update mu
    current_mu <- current_mu + rho*(Kmatrix %*% current_beta - current_theta)
    # CHECK STOPPING CONDITIONS
    # calculate dual residual
    epsilon_dual <- sqrt(length(current_theta))*epsilon_abs + epsilon_rel*sqrt(sum((t(Kmatrix) %*% current_mu)^2))
    current_s <- sqrt(sum((rho*t(Kmatrix) %*% (previous_theta - current_theta))^2))
    # calculate primal residual
    epsilon_primal <- sqrt(length(current_beta))*epsilon_abs + epsilon_rel*max(sqrt(sum((Kmatrix %*% current_beta)^2)), sqrt(sum(current_theta^2)))
    current_r <- sqrt(sum((Kmatrix %*% current_beta - current_theta)^2))
    # break loop if conditions are met
    if(current_s < epsilon_dual & current_r < epsilon_primal) break
  }
  # CALCULATE OUTPUTS
  mse.betas <- mean((current_beta - trueBetas)^2)
  pred.y <- Xmatrix %*% current_beta
  mse.y <- mean((Yvector - pred.y)^2)
  output <- list(stats=c(mse.betas=mse.betas, mse.y=mse.y), beta=current_beta, beta_matrix=beta_matrix[1:t,], theta=current_theta, mu=current_mu)
  return(output)
}
########################################
# END fsgl.fit function
########################################

######################################################
# SET SOME SIMULATION PARAMETERS
######################################################
# number of subjects
n <- 50
# 20*20 grid = 400 pixels
# dim 1 is number of rows of image
dim1 <- 20
# dim 2 is number of columns of image
dim2 <- 20
p <- dim1*dim2

######################################################
# DEFINE GROUP STRUCTURE
# 16 groups of 25 pixels
# NOTE: Uncomment one of three possible group structures, A, B, or C
######################################################
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
Groups <- c(rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5),
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 4,
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 8,
            rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 12)
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
# numbers <- cbind(1:16, c(4:16, 1:3), c(7:16, 1:6), c(10:16, 1:9))
# numbers.list <- apply(numbers, 1, list)
# # function to create subblock (5*5 square)
# subblock <- function(x){
#   x <- unlist(x)
#   data <- c(rep(c(rep(x[1], 3), rep(x[2], 2)), 2), rep(x[1], 3), rep(x[4], 2), rep(c(rep(x[2], 2), x[4], rep(x[3], 2)), 2))
#   return(matrix(data, nrow=5, byrow=TRUE))
# }
# # create subblocks
# subblocks <- lapply(numbers.list, subblock)
# # combine subblocks in rows
# row1 <- cbind(subblocks[[1]], subblocks[[2]], subblocks[[3]], subblocks[[4]])
# row2 <- cbind(subblocks[[5]], subblocks[[6]], subblocks[[7]], subblocks[[8]])
# row3 <- cbind(subblocks[[9]], subblocks[[10]], subblocks[[11]], subblocks[[12]])
# row4 <- cbind(subblocks[[13]], subblocks[[14]], subblocks[[15]], subblocks[[16]])
# Groups <- as.vector(rbind(row1, row2, row3, row4))
# # remove stuff we don't need
# rm('row1', 'row2', 'row3', 'row4', 'numbers', 'numbers.list', 'subblocks')
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
# Groups <- c(rep(1:16, 25))
########################################
# OPTIONAL STEP: Visualize the groups
########################################
# image.plot(matrix(Groups, nrow=20, byrow=TRUE), main='Group Structure')

######################################################
# DEFINE TRUE COEFFICIENTS
# NOTE: Uncomment one of nine possible sets of true coefficients
# Should match group structure chosen above
######################################################
# set coefficient magnitude
beta.value <- 3  
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
### TRUE COEFFICIENTS 1A: Complete Group
trueBetas <- as.numeric(Groups == 7)*beta.value
### TRUE COEFFICIENTS 4A: Sparse Group 
# trueBetas <- as.numeric(Groups == 7)*beta.value
# trueBetas[c(113, 114, 115, 133, 134, 135, 154, 155, 175, 195)] <- 0
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
### TRUE COEFFICIENTS 2B: Complete Group
# trueBetas <- as.numeric(Groups == 10)*beta.value
### TRUE COEFFICIENTS 5B: Sparse Group
# trueBetas <- as.numeric(Groups == 10)*beta.value
# trueBetas[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 7B: Extra Sparse Group 
# trueBetas <- as.numeric(Groups == 10)*beta.value
# trueBetas[c(44, 45, 111, 112, 113, 131, 132, 133, 151, 152, 153, 209, 210, 229, 230, 267, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 8B: Misspecified Group 
# trueBetas <- as.numeric(Groups == 10)*beta.value
# rotate <- function(x) t(apply(x, 2, rev))
# trueBetas <- as.vector(rotate(matrix(trueBetas, nrow=20)))
### TRUE COEFFICIENTS 9B: Misspecified Sparse Group 
# trueBetas <- as.numeric(Groups == 10)*beta.value
# trueBetas[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
# rotate <- function(x) t(apply(x, 2, rev))
# trueBetas <- as.vector(rotate(matrix(trueBetas, nrow=20)))
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
### TRUE COEFFICIENTS 3C: Complete Group
# trueBetas <- as.numeric(Groups == 7)*beta.value
### TRUE COEFFICIENTS 6C: Sparse Group
# trueBetas <- as.numeric(Groups == 7)*beta.value
# set.seed(1)
# trueBetas[which(trueBetas==3)[sample(1:25, 10)]] <- 0
########################################
# OPTIONAL STEP: Visualize the true coefficients
########################################
# image.plot(matrix(trueBetas, nrow=20, byrow=TRUE), main='True Coefficients')

######################################################
# MAKE THE K MATRIX
# dimension N*p  where N = p + d + G
# p is number of pixels (400)
# d is number of rows of D matrix
# G is number of groups
######################################################
# create the J matrix 
# (p*p identity matrix)
########################################
J <- diag(p)
########################################
# create the D matrix 
# for 2D fused lasso
########################################
# matrix block to fuse one row
D.fuse.row.vector <- rep(0, dim2*(dim2-1))
D.fuse.row.vector[(0:(dim2-2)*(dim2) + 1:(dim2-1))] <- -1
D.fuse.row.vector[(0:(dim2-2)*(dim2) + 2:(dim2))] <- 1
D.fuse.row <- matrix(D.fuse.row.vector, nrow=(dim2-1), ncol=dim2, byrow=TRUE)
# make this into a block diagonal matrix
D.fuse.row.list <- paste0('list(', paste(rep('D.fuse.row', dim1), collapse=', '), ')')
big.D.row <- bdiag(eval(parse(text=D.fuse.row.list)))
# matrix block to fuse columns
big.D.column <- matrix(0, nrow=(dim1-1)*dim2, ncol=p)
big.D.column.neg.ones <- 0:(dim1-2)*dim2 + rep(1:dim2, each=(dim1-1))
big.D.column.ones <- 1:(dim1-1)*dim2 + rep(1:dim2, each=(dim1-1))
for(i in 1:dim(big.D.column)[1]){
  big.D.column[i, big.D.column.neg.ones[i]] <- -1
  big.D.column[i, big.D.column.ones[i]] <- 1
}
# now bind these together
D <- rbind(big.D.row, big.D.column)
D <- as.matrix(D)
# remove stuff we don't need
rm('D.fuse.row.vector', 'D.fuse.row', 'D.fuse.row.list', 'big.D.row', 'big.D.column', 'big.D.column.neg.ones', 'big.D.column.ones')
########################################
# create the G matrix
########################################
for(i in 1:16){
  x <- diag(Groups==i)
  assign(paste0('group', i), x[rowSums(x)==1,]*1)
}
# make G matrix
G <- rbind(group1, group2, group3, group4, group5, group6, group7, group8, group9, group10, group11, group12, group13, group14, group15, group16)
########################################
# put together the K matrix
########################################
K <- rbind(J, D, G)

######################################################
# DEFINE SOME MORE QUANTITIES USED IN ALL FITS
# NOTE: For fused sparse group lasso, 
# each individual coefficient is treated as an effective group 
# (there are p of these);
# each difference defined in the D matrix (each row) 
# is treated as an effective group
# (there are d of these);
# and each actual group of coefficients is treated as an effective group
# (there are G of these).
# Total number of effective groups is p + d + G = N
# Size of each group is nj for j in {1,...,N}
######################################################
# make a vector of effective group sample sizes
nj <- c(rep(1, nrow(J)), rep(1, nrow(D)), nrow(group1), nrow(group2), nrow(group3), nrow(group4), nrow(group5), nrow(group6), nrow(group7), nrow(group8), nrow(group9), nrow(group10), nrow(group11), nrow(group12), nrow(group13), nrow(group14), nrow(group15), nrow(group16))
# save N, the total number of effective groups
N <- length(nj)
# save number of rows of J and D
nrowsJ <- nrow(J)
nrowsD <- nrow(D)
# save ngroups, number of actual groups
ngroups <- 16
# store group indices in list
cumsumnj <- cumsum(nj)
group_indices <- list()
for (j in 1:N){
  group_indices[[j]] <- (max(0, cumsumnj[j-1]) + 1):cumsumnj[j]
}
# remove stuff we don't need
rm('group1', 'group2', 'group3', 'group4', 'group5', 'group6', 'group7', 'group8', 'group9', 'group10', 'group11', 'group12', 'group13', 'group14', 'group15', 'group16', 'x', 'i', 'J', 'D', 'G')
# set initial parameters
beta_0 <- rep(0, p)
theta_0 <- rep(0, sum(nj))
mu_0 <- rep(0, sum(nj)) 
# set ADMM step size rho
rho <- 1
# assign subjects to folds
folds <- rep(1:5, 10)
# define lambda
exponent <- seq(3, -3, length = 50)
lambda.grid <- 10 ^ exponent
# define alpha and gamma 
alphavals <- c(rep(0,5), rep(0.2, 4), rep(0.5, 4), rep(0.8, 4), rep(1, 4))
gammavals <- c(0, rep(c(0.2, 0.5, 0.8, 1), 5))
alphagamma <- cbind(alphavals, gammavals)

######################################################
# BEGIN SIMULATION
# For each of the 9 true coefficient scenarios,
# loop over random seeds from 1 to 100.
# For each random seed,
# loop over 21 alpha/gamma combinations.
######################################################
# make an empty results matrix
results <- matrix(NA, ncol=408)
colnames(results) <- c('seed', 'alpha', 'gamma', 'opt.lambda', 'mean.cve', 'sd.cve', 'mse.betas', 'mse.y', paste0('beta', 1:400))

########################################
# GENERATE THE DATA
# NOTE: It takes some time to run simulation 
# for all 100 random seeds.
# Try doing 1 seed first to see how long it takes.
# Can parallelize the seeds.
########################################
start.time <- Sys.time()
# seeds <- 1:100
seeds <- 1
########################################
# BEGIN LOOP OVER RANDOM SEEDS
########################################
for (seed in 1:length(seeds)) {
  set.seed(seeds[seed])
  # generate X matrix
  # iid random standard normal
  X.sample <- matrix(rnorm(n * p), ncol = p)
  # generate the errors
  epsilon <- rnorm(n, mean = 0, sd = 2)
  # generate the Ys
  y.sample <- X.sample %*% trueBetas + epsilon
  # center the columns of the X.sample matrix and y.sample
  # so we won't have to estimate an intercept
  X.all <- scale(X.sample, center = TRUE, scale = FALSE)
  y.all <- scale(y.sample, center = TRUE, scale = FALSE)
  # calculate beta update terms for the full sample
  XtX_KtK <- t(X.all) %*% X.all + rho * t(K) %*% K
  Beta_update_factor <- solve(XtX_KtK)
  Beta_update_term1 <- t(X.all) %*% y.all
  
  ########################################
  # BEGIN LOOP OVER ALPHA/GAMMA VALUES
  ########################################
  for (k in 1:21) {
    # set alpha
    alpha <- alphagamma[k, 1]
    # set gamma
    gamma <- alphagamma[k, 2]
    
    ########################################
    # DO 5-FOLD CROSS VALIDATION
    ########################################
    cve.grid <- matrix(nrow = length(lambda.grid), ncol = 5)
    for (fold in 1:5) {
      # current CV sample
      curX <- X.sample[folds != fold,]
      curY <- y.sample[folds != fold]
      # out of sample fold
      curXout <- X.sample[folds == fold,]
      curYout <- y.sample[folds == fold]
      # center the columns of the curX matrix and curY
      # so we won't have to estimate an intercept
      X <- scale(curX, center = TRUE, scale = FALSE)
      y <- scale(curY, center = TRUE, scale = FALSE)
      curXoutcentered <- scale(curXout, center = TRUE, scale = FALSE)
      curYoutcentered <- scale(curYout, center = TRUE, scale = FALSE)
      # calculate beta update terms for the current CV sample
      XtX_KtK <- t(X) %*% X + rho * t(K) %*% K
      Beta_update_factor_cv <- solve(XtX_KtK)
      Beta_update_term1_cv <- t(X) %*% y
      # start all parameters at zero when lambda is largest
      # thereafter we will use warm starts as we decrease lambda
      current_beta <- beta_0
      current_theta <- theta_0
      current_mu <- mu_0
      # loop over lambda
      for (j in 1:length(lambda.grid)) {
        # define lambda vector
        Lambda <- c(rep(alpha * gamma * lambda.grid[j], nrowsJ),
            rep((1 - gamma) * lambda.grid[j], nrowsD),
            rep((1 - alpha) * gamma * lambda.grid[j], ngroups))
        # fit estimator to the current CV sample
        fit <- fsgl.fit(
          Xmatrix = X,
          Yvector = y,
          Kmatrix = K,
          beta0 = current_beta,
          theta0 = current_theta,
          mu0 = current_mu,
          beta_update_factor = Beta_update_factor_cv,
          beta_update_term1 = Beta_update_term1_cv,
          trueBetas = trueBetas,
          groupindex = group_indices,
          lambda = Lambda,
          rho = rho,
          p = p,
          N = N,
          Niter = 2000,
          epsilon_abs = 10^-3,
          epsilon_rel = 10^-3
          )
        # calculate and save the out of sample MSE
        pred.y.cv <- curXoutcentered %*% fit$beta
        cve.grid[j, fold] <- mean((curYoutcentered - pred.y.cv)^2)
        # update the parameters for warm start at next lambda value
        current_beta <- fit$beta
        current_theta <- fit$theta
        current_mu <- fit$mu
      } # END LOOP OVER LAMBDA
    } # END LOOP OVER CROSS VALIDATION FOLDS
    
    # calculate mean and standard deviation CVE for each lambda
    mean.cve <- rowMeans(cve.grid)
    sd.cve <- apply(cve.grid, 1, sd)
    # calculate mean / sd CVE at optimal lambda
    opt.lambda <- lambda.grid[which.min(mean.cve)]
    opt.mean.cve <- mean.cve[which.min(mean.cve)]
    opt.sd.cve <- sd.cve[which.min(mean.cve)]
    
    ########################################
    # FOR GIVEN ALPHA/GAMMA
    # FIT MODEL TO THE ENTIRE SAMPLE 
    # AT OPTIMAL LAMBDA (minimum mean CVE)
    ########################################
    Lambda <-c(rep(alpha * gamma * lambda.grid[which.min(mean.cve)], nrowsJ),
      rep((1 - gamma) * lambda.grid[which.min(mean.cve)], nrowsD),
      rep((1 - alpha) * gamma * lambda.grid[which.min(mean.cve)], ngroups)
    )
    fit <- fsgl.fit(
      Xmatrix = X.all,
      Yvector = y.all,
      Kmatrix = K,
      beta0 = beta_0,
      theta0 = theta_0,
      mu0 = mu_0,
      beta_update_factor = Beta_update_factor,
      beta_update_term1 = Beta_update_term1,
      trueBetas = trueBetas,
      groupindex = group_indices,
      lambda = Lambda,
      rho = rho,
      p = p,
      N = N,
      Niter = 2000,
      epsilon_abs = 10^-4,
      epsilon_rel = 10^-4
    )
    
    # save outputs
    curResult <- c(seeds[seed], alpha, gamma, opt.lambda, opt.mean.cve, opt.sd.cve, fit$stats[[1]], fit$stats[[2]], fit$beta)
    results <- rbind(results, curResult)
    # NOTE: change filename 'simresults.csv' to reflect the group structure, 
    # true coefficient structure, and random seeds used
    write.csv(results[-1,], paste0(outputdir, 'simresults.csv'), row.names = FALSE)
    # OPTIONAL: print info to track progress
    print(c(seeds[seed], alphagamma[k,]))
  } # END LOOP OVER ALPHA/GAMMA VALUES
} # END LOOP OVER SEEDS

stop.time <- Sys.time()
stop.time - start.time
