#' Unmask function to return to original space.
#'
#' Transforms data from R^m to R^t, where m is the dimension of the masked space and t is the dimension of the total space. Input is the data (m-dimensional vector or a m by q dimensional matrix) to be transformed and a t-dimensional binary mask vector.
#' @param data a m-dimensional vector or a m by q dimensional matrix.
#' @param mask t-dimensional binary mask vector. Zeros indicate the entires or columns to be omitted, ones indicate entries or rows that will be kept.
#' @keywords mask
#' @export
#' @examples 
#' x <- matrix(1:32, nrow=2)
#' Mask <- c(rep(0, 2), rep(1, 8), rep(0, 2), rep(1, 8), rep(0,4))
#' unmask(x, Mask)

unmask <- function(data, mask){
  maskmatrix <- diag(mask)[mask==1,]
  unmaskmatrix <- solve(t(maskmatrix) %*% maskmatrix) %*% t(maskmatrix)
}