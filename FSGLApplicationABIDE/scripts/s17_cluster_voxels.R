####################################################################
# Fused Sparse Group Lasso ABIDE Application
# Do heirarchical clustering on the average pairwise correlation matrix 
# 50 clusters (used in final analysis) (Figure 5C)
# 100 clusters
####################################################################
# Script used for analyses reported in the manuscript
# "Incorporating Prior Information with Fused Sparse Group Lasso:
# Application to Prediction of Clinical Measures from Neuroimages"
### INPUTS: 
# pairwisecorr_mean.txt
# pairwisecorr_voxel_indices.txt
### OUTPUTS:
# diagcorr.txt: indicates which voxels were used (no subjects were missing data)
# voxel_clusters.txt: voxel indices and voxel groups
####################################################################

######################################################
# LOAD DATA
######################################################
setwd('./data')
library(cluster)
# read in data
data <- as.matrix(read.csv('pairwisecorr_mean.txt', header=FALSE))
# read in indices
indices <- read.table('pairwisecorr_voxel_indices.txt', header=FALSE)

######################################################
# CLUSTER VOXELS
######################################################
# visualize
hist(data)
d <- diag(data)
plot(d, type='l')
sum(d==1) 
# 5476 of 6630 voxels have 1 on the diagonal

# remove voxels from the analysis 
# where diagonal is not equal to one
# since some subjects don't have data at this location
data_new <- data[d==1, d==1]
indices_new <- indices[d==1,]

# convert Pearson correlation to 
# Euclidean distance matrix
distance <- sqrt(2*(1 - data_new))
hist(distance)
# make dissimilarity matrix 
dissim <- as.dist(distance)

# do hierarchical clustering
clusters <- hclust(dissim, method='ward.D2', members=NULL)

# plot dendogram
pdf('cluster_plot.pdf')
plot(clusters, labels=FALSE)
dev.off()

# cut dendogram
clust50 <- cutree(clusters, k=50)
clust100 <- cutree(clusters, k=100)

# save d for later 3D viewing of voxels which were excluded
diagcorr <- data.frame(index=indices, diagcorr=d, voxels_kept=(d==1)*1)
colnames(diagcorr) <- c('index', 'diagcorr', 'voxels_kept')
write.table(diagcorr, 'diagcorr.txt', col.names=TRUE, row.names = FALSE)

# save clustering groups 
clustgps <- data.frame(index=indices_new, clust50=clust50, clust100=clust100)
write.table(clustgps, 'voxel_clusters.txt', col.names=TRUE, row.names = FALSE)
