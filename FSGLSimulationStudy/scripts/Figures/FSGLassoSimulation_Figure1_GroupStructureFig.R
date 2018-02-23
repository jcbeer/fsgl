####################################################################
# Fused Sparse Group Lasso Simulation Study
# Make Figure 1: Group Structures and True Coefficients
####################################################################
# This script is used to report simulation methods
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# none
### OUTPUTS:
# pdf figure 'GroupStructTrueCoef.pdf'
####################################################################

library(fields) # for tim.colors() to make rainbow
# set the working directory (where figure will be saved)
# setwd()

########################################
# 9 scenarios:
# 1A. completely aggregated
# 2B. partially aggregated
# 3C. completely distributed
# 4A. sparse group completely aggregated
# 5B. sparse group parially aggregated
# 6C. sparse group completely distributed
# 7B. partially aggregated extra sparse
# 8B. partially aggregated misspecified
# 9B. partially aggregated misspecified sparse
########################################

########################################
# SET SOME SIMULATION PARAMETERS
########################################
# number of subjects
n <- 50
# 20*20 grid = 400 pixels
# dim 1 is number of rows of image
dim1 <- 20
# dim 2 is number of columns of image
dim2 <- 20
p <- dim1*dim2

######################################################
# DEFINE GROUP STRUCTURES
# 16 groups of 25 pixels
######################################################
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
GroupsA <- c(rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5),
             rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 4,
             rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 8,
             rep(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)), 5) + 12)
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
numbers <- cbind(1:16, c(4:16, 1:3), c(7:16, 1:6), c(10:16, 1:9))
numbers.list <- apply(numbers, 1, list)
# function to create subblock (5*5 square)
subblock <- function(x){
  x <- unlist(x)
  data <- c(rep(c(rep(x[1], 3), rep(x[2], 2)), 2), rep(x[1], 3), rep(x[4], 2), rep(c(rep(x[2], 2), x[4], rep(x[3], 2)), 2))
  return(matrix(data, nrow=5, byrow=TRUE))
}
# create subblocks
subblocks <- lapply(numbers.list, subblock)
# combine subblocks in rows
row1 <- cbind(subblocks[[1]], subblocks[[2]], subblocks[[3]], subblocks[[4]])
row2 <- cbind(subblocks[[5]], subblocks[[6]], subblocks[[7]], subblocks[[8]])
row3 <- cbind(subblocks[[9]], subblocks[[10]], subblocks[[11]], subblocks[[12]])
row4 <- cbind(subblocks[[13]], subblocks[[14]], subblocks[[15]], subblocks[[16]])
GroupsB <- as.vector(rbind(row1, row2, row3, row4))
# remove stuff we don't need
rm('row1', 'row2', 'row3', 'row4', 'numbers', 'numbers.list', 'subblocks')
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
GroupsC <- c(rep(1:16, 25))
########################################
# OPTIONAL STEP: Visualize the groups
########################################
# image.plot(matrix(GroupsA, nrow=20, byrow=TRUE), main='Group Structure')

######################################################
# DEFINE TRUE COEFFICIENTS
######################################################
# set coefficient magnitude
beta.value <- 3  
########################################
### GROUP STRUCTURE A: Completely Aggregated
# 16 groups of 5*5 blocks
########################################
### TRUE COEFFICIENTS 1A: Complete Group
trueBetas1A <- as.numeric(GroupsA == 7)*beta.value
### TRUE COEFFICIENTS 4A: Sparse Group 
trueBetas4A <- as.numeric(GroupsA == 7)*beta.value
trueBetas4A[c(113, 114, 115, 133, 134, 135, 154, 155, 175, 195)] <- 0
########################################
### GROUP STRUCTURE B: Partially Aggregated
# 16 groups (1) 3*3 + (3) 2*2 + (2) 1*2 blocks
########################################
### TRUE COEFFICIENTS 2B: Complete Group
trueBetas2B <- as.numeric(GroupsB == 10)*beta.value
### TRUE COEFFICIENTS 5B: Sparse Group
trueBetas5B <- as.numeric(GroupsB == 10)*beta.value
trueBetas5B[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 7B: Extra Sparse Group 
trueBetas7B <- as.numeric(GroupsB == 10)*beta.value
trueBetas7B[c(44, 45, 111, 112, 113, 131, 132, 133, 151, 152, 153, 209, 210, 229, 230, 267, 364, 365, 384, 385)] <- 0
### TRUE COEFFICIENTS 8B: Misspecified Group 
trueBetas8B <- as.numeric(GroupsB == 10)*beta.value
rotate <- function(x) t(apply(x, 2, rev))
trueBetas8B <- as.vector(rotate(matrix(trueBetas8B, nrow=20)))
### TRUE COEFFICIENTS 9B: Misspecified Sparse Group 
trueBetas9B <- as.numeric(GroupsB == 10)*beta.value
trueBetas9B[c(44, 45, 209, 210, 229, 230, 364, 365, 384, 385)] <- 0
rotate <- function(x) t(apply(x, 2, rev))
trueBetas9B <- as.vector(rotate(matrix(trueBetas9B, nrow=20)))
########################################
#### GROUP STRUCTURE C: Completely Distributed
# 16 groups distributed in 1*1 blocks
########################################
### TRUE COEFFICIENTS 3C: Complete Group
trueBetas3C <- as.numeric(GroupsC == 7)*beta.value
### TRUE COEFFICIENTS 6C: Sparse Group
trueBetas6C <- as.numeric(GroupsC == 7)*beta.value
set.seed(1)
trueBetas6C[which(trueBetas6C == 3)[sample(1:25, 10)]] <- 0
########################################
# OPTIONAL STEP: Visualize the true coefficients
########################################
# image.plot(matrix(trueBetas1A, nrow=20, byrow=TRUE), main='True Coefficients')


####################################################################
# MAKE FIGURE: GROUP STRUCTURES  
# AND TRUE COEFFICIENTS
####################################################################
pdf('GroupStructTrueCoef.pdf', width=6, height=10, title='FSGLasso Simulation Group Structure and True Coefficients')
# plot group structures
par(mfrow=c(4, 3), mar=c(1, 1, 4, 1), pty='s')
# ROW 1 Group structures
image(matrix(GroupsA, nrow=20, byrow=TRUE), main='Group Structure A:\nCompletely Aggregated', yaxt='n', xaxt='n', col=tim.colors(16))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(GroupsB, nrow=20, byrow=TRUE), main='Group Structure B:\nPartially Aggregated', yaxt='n', xaxt='n', col=tim.colors(16))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(GroupsC, nrow=20, byrow=TRUE), main='Group Structure C:\nCompletely Distributed', yaxt='n', xaxt='n', col=tim.colors(16))
grid(20, 20, lty=1, col='white', lwd=0.5)
# ROW 2 True coefficients complete groups
image(matrix(trueBetas1A, nrow=20, byrow=TRUE), main='True Coefficients 1A:\nComplete Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas2B, nrow=20, byrow=TRUE), main='True Coefficients 2B:\nComplete Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas3C, nrow=20, byrow=TRUE), main='True Coefficients 3C:\nComplete Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
# ROW 3 True coefficients sparse groups
image(matrix(trueBetas4A, nrow=20, byrow=TRUE), main='True Coefficients 4A:\nSparse Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas5B, nrow=20, byrow=TRUE), main='True Coefficients 5B:\nSparse Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas6C, nrow=20, byrow=TRUE), main='True Coefficients 6C:\nSparse Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
# ROW 4 True coefficients additional
image(matrix(trueBetas7B, nrow=20, byrow=TRUE), main='True Coefficients 7B:\nExtra Sparse Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas8B, nrow=20, byrow=TRUE), main='True Coefficients 8B:\nMisspecified Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
image(matrix(trueBetas9B, nrow=20, byrow=TRUE), main='True Coefficients 9B:\nMisspecified Sparse Group', yaxt='n', xaxt='n', col=c('white', 'black'))
grid(20, 20, lty=1, col='white', lwd=0.5)
dev.off()

