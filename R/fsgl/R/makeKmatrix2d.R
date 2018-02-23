#' Make 2D Fused K Matrix
#'
#' @description Builds the K matrix for fused sparse group lasso \code{fsgl.fit} function where graph structure is 2D fused. Requires Matrix package.
#' @param dim1 number of rows of image.
#' @param dim2 number of columns of image.
#' @param groups vector encoding group membership, with groups indexed as 1, 2, etc. Note: columns of image are concatenated from left to right, group membership coding should correspond.
#' @keywords lasso
#' @export
#' @examples 
#' makeKmatrix2d(dim1=3, dim2=3, groups=c(1, 2, 3, 2, 2, 1, 3, 3, 1))


makeKmatrix2d <- function(dim1, dim2, groups){
  ### make J matrix for lasso
  p <- dim1*dim2
  J <- diag(p)
  ### make D matrix for 2D fused lasso
  # matrix block to fuse one colunm
  D.fuse.column.vector <- rep(0, dim1*(dim1-1))
  D.fuse.column.vector[(0:(dim1-2)*(dim1) + 1:(dim1-1))] <- -1
  D.fuse.column.vector[(0:(dim1-2)*(dim1) + 2:(dim1))] <- 1
  D.fuse.column <- matrix(D.fuse.column.vector, nrow=(dim1-1), ncol=dim1, byrow=TRUE)
  # make this into a block diagonal matrix
  D.fuse.column.list <- paste0('list(', paste(rep('D.fuse.column', dim2), collapse=', '), ')')
  big.D.column <- Matrix::bdiag(eval(parse(text=D.fuse.column.list)))
  # matrix block to fuse rows
  big.D.row <- matrix(0, nrow=(dim2-1)*dim1, ncol=p)
  big.D.row.neg.ones <- 0:(dim2-2)*dim1 + rep(1:dim1, each=(dim2-1))
  big.D.row.ones <- 1:(dim2-1)*dim1 + rep(1:dim1, each=(dim2-1))
  for(i in 1:dim(big.D.row)[1]){
    big.D.row[i, big.D.row.neg.ones[i]] <- -1
    big.D.row[i, big.D.row.ones[i]] <- 1
  }
  # now bind these together
  D <- rbind(big.D.column, big.D.row)
  D <- as.matrix(D)
  ### make G matrix for groups
  G <- diag(groups==1)
  for(i in 2:max(groups)){
    G <- rbind(G, diag(groups==i))
  }
  # remove zero rows
  G <- G[rowSums(G)==1,]
  # put together into K matrix
  K <- rbind(J, D, G)
  return(list(K=K, nj=nrow(J), nd=nrow(D), ngroups=max(groups), groupsizes=as.vector(table(groups))))
}