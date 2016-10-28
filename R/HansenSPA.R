#' Execute the Hansen, P. R. (2005) Superior Predictive Ability Test
#'
#' @param Dmat A matrix with the relative performance variables (NxK)
#' @param bVec A vector with the relative performance benchmark (Nx1)
#' @param typeFunc Type of function to be used:
#' 0 - max(0,x)
#' 1 - x*1_{x >= -sqrt{(w^2/n)2log log n}}
#' 2 - x
#' @param B total of boostraps samples.
#' @return Statistic of the test and p-value
#' @examples
#' add(1, 1)
#' add(10, 1)
hansen.spa <- function(Dmat,bVec,typeFunc,B) {

}
