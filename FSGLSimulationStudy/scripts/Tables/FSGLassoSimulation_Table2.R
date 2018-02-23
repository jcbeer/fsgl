####################################################################
# Fused Sparse Group Lasso Simulation Study
# Make Table 2: Combinations of (alpha, gamma) yielding the most frequent 
# lowest error out of 100 simulation replications
####################################################################
# This script is used to analyze the simulation results
# reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
## original simulation statistics plus test set MSE
# compaggplot.csv, partaggplot.csv, compdistplot.csv,
# spcompaggplot.csv, sppartaggplot.csv, spcompdistplot.csv, 
# exspplot.csv, misspplot.csv, misspspplot.csv
### OUTPUTS:
# text printed to console that can be cut and pasted into LaTeX document
####################################################################

library(xtable) # creates LaTeX tables

########################################
# LOAD DATASETS
########################################
# set the working and data directories
# setwd()
datadir1 <- paste0(getwd(), '/DataForPlots/')
compagg <- read.csv(paste0(datadir1, 'compaggplot.csv'))
partagg <- read.csv(paste0(datadir1, 'partaggplot.csv'))
compdist <- read.csv(paste0(datadir1, 'compdistplot.csv'))
spcompagg <- read.csv(paste0(datadir1, 'spcompaggplot.csv'))
sppartagg <- read.csv(paste0(datadir1, 'sppartaggplot.csv'))
spcompdist <- read.csv(paste0(datadir1, 'spcompdistplot.csv'))
exsp <- read.csv(paste0(datadir1, 'exspplot.csv'))
missp <- read.csv(paste0(datadir1, 'misspplot.csv'))
misspsp <- read.csv(paste0(datadir1, 'misspspplot.csv'))

########################################
# function to determine optimal alpha, gamma combo 
# based on most frequent lowest error
########################################
optparam <- function(data, errorvar){
  lowest <- rep(NA, 100)
  for (i in 1:100){
    seeddata <- data[data$seed==i,]
    seeddata$alphagamma <- as.character(interaction(seeddata$alpha, seeddata$gamma, sep=', '))
    lowest[i] <- seeddata$alphagamma[which.min(seeddata[,errorvar])]
  }
  mostfreqlowest <- names(which.max(table(lowest)))
  return(mostfreqlowest)
}

########################################
# determine optimal alpha, gamma combos 
# based on most frequent lowest error
# mean CVE
########################################
meancve <- c(
  paste0('(', optparam(compagg, 'mean.cve'), ')'),
  paste0('(', optparam(partagg, 'mean.cve'), ')'),
  paste0('(', optparam(compdist, 'mean.cve'), ')'),
  paste0('(', optparam(spcompagg, 'mean.cve'), ')'),
  paste0('(', optparam(sppartagg, 'mean.cve'), ')'),
  paste0('(', optparam(spcompdist, 'mean.cve'), ')'),
  paste0('(', optparam(exsp, 'mean.cve'), ')'),
  paste0('(', optparam(missp, 'mean.cve'), ')'),
  paste0('(', optparam(misspsp, 'mean.cve'), ')')
)

########################################
# determine optimal alpha, gamma combos 
# based on most frequent lowest error
# mse.beta
########################################
msebeta <- c(
  paste0('(', optparam(compagg, 'mse.betas'), ')'),
  paste0('(', optparam(partagg, 'mse.betas'), ')'),
  paste0('(', optparam(compdist, 'mse.betas'), ')'),
  paste0('(', optparam(spcompagg, 'mse.betas'), ')'),
  paste0('(', optparam(sppartagg, 'mse.betas'), ')'),
  paste0('(', optparam(spcompdist, 'mse.betas'), ')'),
  paste0('(', optparam(exsp, 'mse.betas'), ')'),
  paste0('(', optparam(missp, 'mse.betas'), ')'),
  paste0('(', optparam(misspsp, 'mse.betas'), ')')
)

########################################
# determine optimal alpha, gamma combos 
# based on most frequent lowest error
# mse y.test
########################################
mseytest <- c(
  paste0('(', optparam(compagg, 'mse.y.test'), ')'),
  paste0('(', optparam(partagg, 'mse.y.test'), ')'),
  paste0('(', optparam(compdist, 'mse.y.test'), ')'),
  paste0('(', optparam(spcompagg, 'mse.y.test'), ')'),
  paste0('(', optparam(sppartagg, 'mse.y.test'), ')'),
  paste0('(', optparam(spcompdist, 'mse.y.test'), ')'),
  paste0('(', optparam(exsp, 'mse.y.test'), ')'),
  paste0('(', optparam(missp, 'mse.y.test'), ')'),
  paste0('(', optparam(misspsp, 'mse.y.test'), ')')
)

########################################
# scenarios
########################################
scenarios <- c(
  '1A. Completely aggregated',
  '2B. Partially aggregated',
  '3C. Completely distributed',
  '4A. Sparse completely aggregated',
  '5B. Sparse partially aggregated',
  '6C. Sparse completely distributed',
  '7B. Extra sparse partially aggregated',
  '8B. Misspecified partially aggregated',
  '9B. Misspecified sparse partially aggregated'
)

########################################
# make table
# NOTE: There was a tie for 7B test MSE
# (See Supplementary Table S8)
# (0.8, 0.8) needs to be added manually
########################################
header <- c('True coefficient scenario', 'Mean CVE', 'MSE Beta', 'MSE Y Test')
table2 <- cbind(scenarios, meancve, msebeta, mseytest)
colnames(table2) <- header
print(xtable(table2, caption='Combinations of ($\\alpha$, $\\gamma$) yielding the most frequent lowest error out of 100 simulation replications', label='tab:simresults', align=c('l', 'l', 'r', 'r', 'r')), include.rownames=FALSE)
