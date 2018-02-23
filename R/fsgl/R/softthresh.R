#' Soft-thresholding function
#'
#' Soft-threshold a scalar or vector.
#' @param a a scalar or vector.
#' @param kappa soft-thresholding parameter.
#' @keywords soft thresholding
#' @export
#' @examples 
#' softthresh(c(3, 4), 1)

softthresh <- function(a, kappa){
  if(sum(a == 0) == length(a)){
    return(a)
  } else {
    return(max(0, (1 - kappa/sqrt(sum(a^2))))*a)
  }
}