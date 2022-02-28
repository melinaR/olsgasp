#' @title Mater5_2 kernel
#' @description
#' Compute the matern5_2 kernel
#' @param d the distance.
#' @param beta the inverse of range parameter.
#' @examples
#' 
#'matern_5_2_kernel(0.5, beta = 2)
#'  
#' @export
#'

matern_5_2_kernel<-function(d,beta){
  res=(1+sqrt(5)*beta*d + 5*beta^2*d^2/3 )*exp(-sqrt(5)*beta*d)
  res
}