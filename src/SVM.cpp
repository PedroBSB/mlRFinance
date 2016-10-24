#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
using namespace Rcpp;

/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/

//Define the KernelMatrix function
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,
                                        std::string stringValue,
                                        arma::vec parms);
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
//Add some noise to the matrix diagonal
void addNoise(Eigen::MatrixXd &mat,double noise);
//Print Object at Console
void PrintObject(Eigen::MatrixXd mat);
void PrintObject(Eigen::VectorXd vec);


/***********************************************************************************************/
/*********************************     SVM FUNCTIONS    ****************************************/
/***********************************************************************************************/

//' @name CSVML1
//' @title C-SVM L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVM L1:
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' 0<=u<=C1
//'
//' where d is the vector of dependent variable,
//' and Q=K.*(d*t(d))=DKD. C is the Cost parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
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
Rcpp::List CSVML1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, arma::vec parms){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0;
  //LHS equality
  Eigen::MatrixXd CE;
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(y.size(),y.size());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(y.size());
  ci2.fill(C);
  //Append RHS
  Eigen::VectorXd ci0(2.0*y.size());
  ci0 << ci1, ci2;
  //Append LHS
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
  //Diagonal matrix
  Eigen::VectorXd me(y.size());
  me.fill(-1.0);
  Eigen::MatrixXd mI = me.asDiagonal();
  //Vertical concatenation
  CI << CI1,
        mI;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  nearPositiveDefinite(Q,1e-10);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE.transpose(),ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//' @name CSVML2
//' @title C-SVM L2 - Support Vector Regression with C cost and L2 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVM L2:
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' u>=0
//'
//' where d is the vector of dependent variable,
//' and Q=(K+I/C).*(t(d)*d). C is the Cost parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML2(d, A, 1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List CSVML2(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, arma::vec parms){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0;
  //LHS equality
  Eigen::MatrixXd CE;
  //RHS: Inequality 1
  Eigen::VectorXd ci0 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI = Eigen::MatrixXd::Identity(y.size(),y.size());
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  K=K+(CI/C);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  nearPositiveDefinite(Q,1e-10);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE,ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//' @name nuSVM
//' @title nu-SVM - Support Vector Regression with nu parameter.
//' @description The m-support vector classification (Scholkopf, Smola,
//' Williamson, & Bartlett, 2000) uses a new parameter nu which controls
//' the number of support vectors and training errors. The
//' parameter nu in (0, 1] is an upper bound on the fraction of training
//' errors and a lower bound of the fraction of support vectors.
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' d^{t}*u=0
//' nu <=1^t*u
//' 0<=u<=1/l
//'
//' where d is the vector of dependent variable,
//' and Q=K.*(t(d)*d). nu is the parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- nuSVM(d, A, 0.2, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine @chang2001training
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List nuSVM(Eigen::VectorXd y, Eigen::MatrixXd X, double nu, std::string kernel, arma::vec parms){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(1);
  //LHS equality
  Eigen::MatrixXd CE(1,y.size());
  CE = y;
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(y.size(),y.size());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(y.size());
  ci2.fill(1.0/y.size());
  //Append RHS
  Eigen::VectorXd ci0(2.0*y.size()+1);
  ci0 << ci1, ci2, -nu;
  //Append LHS
  Eigen::MatrixXd CI3(1,y.size());
  CI3.fill(1);
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows()+1, CI1.cols());
  CI <<        CI1,
        (-1.0)*CI1,
               CI3;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  nearPositiveDefinite(Q,1e-10);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE,ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


// [[Rcpp::export]]
Eigen::VectorXd solveTest(Eigen::MatrixXd Dmat, Eigen::VectorXd dvec, Eigen::MatrixXd Amat,Eigen::VectorXd bvec, Eigen::MatrixXd CE, Eigen::VectorXd ce){
  //Convert to the quadprog
  Eigen::VectorXd d0=-dvec;
  Eigen::VectorXd b0=-bvec;
  Eigen::VectorXd ce0=-ce;
  //Get the solution
  Eigen::VectorXd x = rcppeigen_quadratic_solve(Dmat,d0, CE,ce0, Amat, b0);
  return(x);
}


// [[Rcpp::export]]
Eigen::VectorXd solveTest2(Eigen::MatrixXd Dmat, Eigen::VectorXd dvec, Eigen::MatrixXd Amat,Eigen::VectorXd bvec){
  Eigen::MatrixXd CE;
  Eigen::VectorXd ce;
  //Convert to the quadprog
  Eigen::VectorXd d0=-dvec;
  Eigen::VectorXd b0=-bvec;
  Eigen::VectorXd ce0=-ce;
  //Get the solution
  Eigen::VectorXd x = rcppeigen_quadratic_solve(Dmat,d0, CE,ce0, Amat, b0);
  return(x);
}
