#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>

/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/
//Define the Solver for Quadratic Programming
Eigen::VectorXd rcppeigen_quadratic_solve(Eigen::MatrixXd & G,
                                          Eigen::VectorXd & g0,
                                          const Eigen::MatrixXd & CE,
                                          const Eigen::VectorXd & ce0,
                                          const Eigen::MatrixXd & CI,
                                          const Eigen::VectorXd & ci0);



/***********************************************************************************************/
/*********************************     QPFS AlgoRithm   ****************************************/
/***********************************************************************************************/

//' @name QPFS
//' @title QPFS - Quadratic Programming Feature Selection.
//' @description Feature Selection:
//'
//' @param Q Similarity Matrix. Dimension equal PxP.
//' @param f Relevance vector. Dimension equal Px1
//' @param alpha Weight between Similarity and Relevance (Keep NA to automatic choose alpha)
//' @return Weights for each variable presented in matrix X
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML1(d, A, 1, 0.1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Eigen::VectorXd QPFS(Eigen::MatrixXd Q, Eigen::VectorXd f, double alpha){

  //nearPositiveDefinite(Q,1e-10);
  Q = nearPDefinite(Q, 1e+6, 1e-06, 1e-07, 1e-08, true);

  //Compute the automatic alpha
  if(std::isnan(alpha)){
    double fbar = f.mean();
    double qbar = Q.mean();
    alpha = qbar/(qbar+fbar);
  }
  Q = (1.0-alpha)*Q;
  f = alpha*f;

  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(Q.cols(),Q.cols());

  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(Q.cols());

  //RHS equality 2
  Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(1);
  ce0.fill(-1.0);
  //LHS equality 2
  Eigen::MatrixXd CE0(1,Q.cols());
  CE0.fill(1.0);

  //Get the weights
  Eigen::VectorXd vecWeights = rcppeigen_quadratic_solve(Q, f, CE0.transpose(),ce0, CI1.transpose(), ci1);

  //Return the results
  return (vecWeights);
}
