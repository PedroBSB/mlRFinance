#' @name CSVML1
#' @title C-SVM L1 - Support Vector Regression with C cost and L1 regularization.
#' @description Optimize the Lagrange multiplier for the C-SVM L1:
#'
#' Min (1/2)u^{t}Qu-1^{t}u
#' s.t.
#' 0<=u<=C1
#'
#' where d is the vector of dependent variable,
#' and Q=K.*(d*t(d))=DKD. C is the Cost parameter.
#'
#' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
#' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
#' @param C Cost parameter. Should be C>0.
#' @param kernel Name of the kernel that will be used.
#' @param parms Parameters associated with chosen kenel.
#' @return List Support Vectors, Kernel used and parameters.
#' If the results for the Support Vectors are NaN it means that
#' there is no Support Vector and the Quadratic Programming Problem
#' is unfeasible.
#' @examples
#'
#' A<-matrix(c(1,2,5,6,
#' 2,4,1,2),nrow=4,ncol=2)
#' d<-c(-1,-1,+1,-1)
#' svm1<- CSVML1(d, A, 1, "Gaussian", c(0.5))
#'
#' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
svm<-function(y, X, C, kernel, parms, biasTerm=TRUE){
  #Stop cases
  stopifnot(inherits(y, "numeric"))
  stopifnot(inherits(X, "matrix"))
  stopifnot(inherits(C, "numeric"))
  stopifnot(inherits(kernel, "character"))
  stopifnot(inherits(parms, "numeric"))
  stopifnot(inherits(biasTerm, "logical"))
  #y must have the same length as X
  if(length(y)!=nrow(X)){
    stop('The argument "y" must have the same length as the argument "X"')
  }
#  if(n != round(n)){
#    warning('The "factorial" function was used')
#    return(factorial(n))
#  }

  #Calculate the CSVML1
  res<- mlRFinance:::CSVML1(y, X, C, kernel, biasTerm)
  #Define the class
  class(res)<-c("list","mlRFinance")
  return(res)
}




