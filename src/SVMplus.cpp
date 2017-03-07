#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]


/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/

//Define the KernelMatrix function
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,
                                        std::string stringValue,
                                        Eigen::RowVectorXd parms);
//Define the Solver for Quadratic Programming
Eigen::VectorXd rcppeigen_quadratic_solve(Eigen::MatrixXd & G,
                                          Eigen::VectorXd & g0,
                                          const Eigen::MatrixXd & CE,
                                          const Eigen::VectorXd & ce0,
                                          const Eigen::MatrixXd & CI,
                                          const Eigen::VectorXd & ci0);
//Test if the matrix is Positive Definite
bool IsPositiveDefinite(Eigen::MatrixXd mat);
//nearest positive semidefinite matrix in terms of Frobenius norm
void nearPositiveDefinite(Eigen::MatrixXd &mat,double noise);
//Nearest positive semidefinite matrix (Matrix::nearPD)
Eigen::MatrixXd nearPDefinite(Eigen::MatrixXd mat, int maxit, double eigtol, double conv_tol, double posd_tol, bool keepDiagonal);
//Add some noise to the matrix diagonal
void addNoise(Eigen::MatrixXd &mat,double noise);
//Print Object at Console
void PrintObject(Eigen::MatrixXd mat);
void PrintObject(Eigen::VectorXd vec);


/************************************ C-SVM+ L1 *************************************************/
// [[Rcpp::export]]
Rcpp::List CSVMplusL1(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, double C, double gamma, double kappa, std::string kernel, Eigen::RowVectorXd parms,std::string kernelStar, Eigen::RowVectorXd parmsStar, bool biasTerm){
  //Support Vectors
  Eigen::VectorXd SV(2*y.size());
  //Create g0
  Eigen::VectorXd g1 = Eigen::VectorXd::Ones(y.size());
  g1=(-1.0)*g1;
  Eigen::VectorXd g2 = Eigen::VectorXd::Zero(y.size());
  Eigen::VectorXd g0(2.0*y.size());
  g0 << g1, g2;

  //Diagonal matrix
  Eigen::MatrixXd yDiag = y.asDiagonal();

  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create the Kernel Matrix
  Eigen::MatrixXd Kstar = KernelMatrixComputation(Z,kernelStar,parmsStar);
  Eigen::MatrixXd KstarY =yDiag*Kstar*yDiag;
  //Crate H matrix
  Eigen::MatrixXd H11 = yDiag*K*yDiag + gamma*KstarY;
  Eigen::MatrixXd H12 = -gamma*KstarY;
  Eigen::MatrixXd H = Eigen::MatrixXd(H11.rows()+H12.cols(),H11.cols()+H12.cols());
  H << H11, H12,
       H12.transpose(),-H12;

  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(Q,1e-10);
  H = nearPDefinite(H, 1e+6, 1e-06, 1e-07, 1e-08, true);

  //Equality Constraint
  Eigen::VectorXd col1 =  Eigen::VectorXd(2*y.size());
  Eigen::VectorXd col2 =  Eigen::VectorXd(2*y.size());
  col1 << y, g2;
  col2 << g2, y;
  Eigen::MatrixXd CE = Eigen::MatrixXd(2*y.size(),2);
  CE << col1, col2;
  //RHS equality
  Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(2);

  //Inequality Constraint
  Eigen::VectorXd ci0 = Eigen::VectorXd::Zero(2*y.size());
  Eigen::VectorXd ci11 = Eigen::VectorXd(y.size());
  ci11.fill(-kappa*C);
  Eigen::VectorXd ci12 = Eigen::VectorXd(y.size());
  ci12.fill(-C);
  ci0<<ci11,ci12;
  //LHS: Inequality 1
  Eigen::MatrixXd CI = Eigen::MatrixXd::Identity(4*y.size(),2*y.size());
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(2*y.size(),2*y.size());
  CI <<  CI1,
        -CI1;

  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(H, g0, CE.transpose(), ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}
