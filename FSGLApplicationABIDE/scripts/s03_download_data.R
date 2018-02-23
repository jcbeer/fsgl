####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Download fMRI data from ABIDE website
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# fileids.csv
# phenosrsall.csv
### OUTPUTS:
# downloads fMRI data into individual subject directories 
# for 219 subjects
####################################################################

######################################################
# DEFINE INPUT AND OUTPUT DIRECTORIES
######################################################
# set working directory where input files are
# setwd()
# set top level output directory
outputdir <- ''

######################################################
# LOAD FILE IDs
######################################################
# load files
fileids <- read.csv('fileids.csv', stringsAsFactors=FALSE)
phenosrsall <- read.csv('phenosrsall.csv')
subjids <- phenosrsall[,3]

######################################################
# ABIDE fMRI data URL format
######################################################
# 'https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/[pipeline]/[strategy]/[derivative]/[file identifier]_[derivative].[ext]'
# [pipeline] = ccs | cpac | dparsf | niak 
# [strategy] = filt_global | filt_noglobal | nofilt_global | nofilt_noglobal
# [file identifier] = the FILE_ID value from the summary spreadsheet
# [derivative] = alff | degree_binarize | degree_weighted | dual_regression | ... 
#                eigenvector_binarize | eigenvector_weighted | falff | func_mask | ... 
#                func_mean | func_preproc | lfcd | reho | rois_aal | rois_cc200 | ... 
#                rois_cc400 | rois_dosenbach160 | rois_ez | rois_ho | rois_tt | vmhc
# [ext] = 1D | nii.gz

######################################################
# NOTE: 'Yale' and 'Leuven' sites did not use all
# caps in their filenames, so replace these parts
######################################################
fileids2 <- as.data.frame(sapply(fileids, gsub, pattern="YALE", replacement="Yale"))
fileids3 <- as.data.frame(sapply(fileids2, gsub, pattern="LEUVEN", replacement="Leuven"))
ids <- cbind(fileids3, subjids)

######################################################
# Create list of the URLs
######################################################
fileurls <- apply(ids, 1, function(x) paste0('https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/ccs/filt_noglobal/func_preproc/', x[1], '_func_preproc.nii.gz'))

######################################################
# LOOP over subjects
######################################################
for(i in 1:219){
  # define name of directory for the subject
  subjdir <- paste0(outputdir, '/s', ids[i,2])
  # create directory for the subject
  dir.create(subjdir)
  # define name of file for the subject
  subjfile <- paste0(subjdir, '/s', ids[i,2], '_func_preproc.nii.gz')
  # download the subject's fMRI data
  download.file(fileurls[i], destfile=subjfile, method='curl')
}


