#' Make 3D Fused K Matrix
#'
#' @description Builds the K matrix for fused sparse group lasso \code{fsgl.fit} function where graph structure is 3D fused. Requires Matrix package.
#' @param dim1 dimension along x axis.
#' @param dim2 dimension along y axis.
#' @param dim3 dimension along z axis.
#' @param mask binary vector indicating masked areas (1s) zeros (0s). 
#' @param groups vector encoding group membership, with groups indexed as 1, 2, etc. If there is a mask, put zeros or other placeholder at location of masked voxels. Note: columns of each slice of image are concatenated from left to right, group membership coding should correspond.
#' @keywords lasso
#' @export
#' @examples 
#' Kmatrix <- makeKmatrix3d(dim1=2, dim2=3, dim3=4, mask=c(rep(0, 4), rep(1, 16), rep(0, 4)), groups=c(rep(0, 4), rep(1, 4), rep(2, 4),  rep(3, 4), rep(4, 4), rep(0, 4)))


makeKmatrix3d <- function(dim1, dim2, dim3, mask=NULL, groups){
  ### make J matrix for lasso penalty
  if(!is.null(mask)){
    J <- diag(sum(mask))
  } else {
    J <- diag(dim1*dim2*dim3)
  }
  ### make D matrix for 3D fused lasso
  # matrix block to fuse one colunm
  D.fuse.column.vector <- rep(0, dim1*(dim1-1))
  D.fuse.column.vector[(0:(dim1-2)*(dim1) + 1:(dim1-1))] <- -1
  D.fuse.column.vector[(0:(dim1-2)*(dim1) + 2:(dim1))] <- 1
  D.fuse.column <- matrix(D.fuse.column.vector, nrow=(dim1-1), ncol=dim1, byrow=TRUE)
  # make this into a block diagonal matrix
  D.fuse.column.list <- paste0('list(', paste(rep('D.fuse.column', dim2), collapse=', '), ')')
  big.D.column <- Matrix::bdiag(eval(parse(text=D.fuse.column.list)))
  # remove some stuff
  rm(list=c('D.fuse.column.vector', 'D.fuse.column', 'D.fuse.column.list'))
  # matrix block to fuse rows
  big.D.row <- matrix(0, nrow=(dim2-1)*dim1, ncol=dim1*dim2)
  big.D.row.neg.ones <- 0:(dim2-2)*dim1 + rep(1:dim1, each=(dim2-1))
  big.D.row.ones <- 1:(dim2-1)*dim1 + rep(1:dim1, each=(dim2-1))
  for(i in 1:dim(big.D.row)[1]){
    big.D.row[i, big.D.row.neg.ones[i]] <- -1
    big.D.row[i, big.D.row.ones[i]] <- 1
  }
  # remove some stuff
  rm(list=c('big.D.row.neg.ones', 'big.D.row.ones'))
  # now bind these together and make a block diagonal
  D.row.col <- rbind(big.D.column, big.D.row)
  rm(list=c('big.D.column', 'big.D.row'))
  D.row.col.list <- paste0('list(', paste(rep('D.row.col', dim3), collapse=', '), ')')
  big.D.row.col <- Matrix::bdiag(eval(parse(text=D.row.col.list)))
  # remove some stuff
  rm(list=c('D.row.col', 'D.row.col.list'))
  # matrix block to fuse slices
  D.slice <- diag(dim1*dim2*dim3)[1:(dim1*dim2*(dim3-1)),]
  neg.ones <- (1 + dim1*dim2):(dim1*dim2*dim3)
  D.slice[cbind(1:(dim1*dim2*(dim3-1)), neg.ones)] <- -1
  # add this to existing D matrix
  D <- rbind(big.D.row.col, D.slice)
  # remove some stuff
  rm(list=c('big.D.row.col', 'D.slice', 'neg.ones'))
  D <- as.matrix(D)
  # remove non-masked 
  if(!is.null(mask)){
    keep.rows <- which(rowSums(abs(D[,mask==0])) == 0)
    D <- D[keep.rows, mask==1]
  }
  ### make G matrix for groups
  if(!is.null(mask)){
    groups <- groups[mask==1]
  }
  G <- diag(groups==1)*1
  for(i in 2:max(groups)){
    G <- rbind(G, diag(groups==i)*1)
  }
  # remove zero rows
  G <- G[rowSums(G)==1,]
  # put together into K matrix
  K <- rbind(J, D, G)
  return(list(K=K, nj=nrow(J), nd=nrow(D), ngroups=max(groups), groupsizes=as.vector(table(groups))))
} # end function definition