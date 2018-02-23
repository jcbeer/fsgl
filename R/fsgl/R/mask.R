#' Mask function to reduce dimension.
#'
#' Transforms data from R^t to R^m, where t is the dimension of the total space and m is the dimension of the masked space. Input is the data (t-dimensional vector or a t by q dimensional matrix) to be transformed and a t-dimensional binary mask vector.
#' @param data a t-dimensional vector or a t by q dimensional matrix.
#' @param mask t-dimensional binary mask vector. Zeros indicate the entires or columns to be omitted, ones indicate entries or rows that will be kept.
#' @keywords mask
#' @export
#' @examples 
#' x <- c(rep(0,4), 1:16, rep(0,4))
#' Mask <- c(rep(0, 4), rep(1, 16), rep(0, 4))
#' mask(x, Mask)

mask <- function(data, mask){
  maskmatrix <- diag(mask)[mask==1,]
  maskmatrix %*% data
}