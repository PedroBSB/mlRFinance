#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
#include "progress.hpp"
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]



/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/

//Define the KernelMatrix function
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,
                                        std::string stringValue,
                                        Eigen::RowVectorXd parms);

//Nearest positive semidefinite matrix (Matrix::nearPD)
Eigen::MatrixXd nearPDefinite(Eigen::MatrixXd mat, int maxit, double eigtol, double conv_tol, double posd_tol, bool keepDiagonal);

/***********************************************************************************************/
/*********************************     KPCA FUNCTIONS    ***************************************/
/***********************************************************************************************/

//' @name KPCA
//' @title Kernel PCA
//' @description Compute Kernel PCA
//'
//'
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List with Principal Componentes, Variability, Eigenvectors and EigenValues
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML1(d, A, 1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Eigen::MatrixXd KPCAMatrix(Eigen::MatrixXd X, std::string kernel, Eigen::RowVectorXd parms){
  //Compute the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Identity
  Eigen::MatrixXd Eye = Eigen::MatrixXd::Identity(X.rows(),X.rows());
  //Ones
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(X.rows(),X.rows());
  Ones.fill(1.0/X.rows());
  //M matrix
  Eigen::MatrixXd M(X.rows(),X.rows());
  M=Eye-Ones;
  //Center the Kernel Matrix
  Eigen::MatrixXd Kmod(X.rows(),X.rows());
  Kmod = M*K*M;

  //Return the results
  return Kmod;
}

