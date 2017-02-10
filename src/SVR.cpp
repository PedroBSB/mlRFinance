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
                                        Eigen::RowVectorXd parms);
Eigen::MatrixXd KernelMatrixComputationPred(Eigen::MatrixXd datMat,
                                        Eigen::RowVectorXd predMat,
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


/***********************************************************************************************/
/*********************************     SVR FUNCTIONS    ****************************************/
/***********************************************************************************************/

//' @name CSVRL1
//' @title C-SVR L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVR L1:
//'
//' Min (1/2)u^{t}Qu+g^{t}u
//' s.t.
//' 0<=u<=C1
//'
//' where u=(lambda*,lambda), g=(e-y,e+y)
//' and Q=|K -K|.
//'       |-K K|
//' C is the Cost parameter and e (epsilon) is the insentitive band
//'
//' @param y Vector with dependent variables. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param epsilon Insentitive band. Should be epsilon>0.
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
//' svm1<- CSVML1(d, A, 1, 0.1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List CSVRL1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, double epsilon, std::string kernel, Eigen::RowVectorXd parms){
  //Support Vectors
  Eigen::VectorXd SV(2*y.size());
  //Create the one vector 2Nx1
  Eigen::VectorXd yfull = Eigen::VectorXd(2*y.size());
  yfull<< (-1.0)*y, (+1.0)*y;
  Eigen::VectorXd evec = Eigen::VectorXd(2*y.size());
  evec.fill((+1.0)*epsilon);
  Eigen::VectorXd g = Eigen::VectorXd(2*y.size());
  g = yfull+evec;
  //RHS equality
  Eigen::VectorXd ce0(1);
  ce0.fill(-0.0);
  //LHS equality
  Eigen::MatrixXd Ones = Eigen::MatrixXd::Ones(1, y.size());
  Eigen::MatrixXd CE(1,2*y.size());
  CE <<   Ones, -Ones;
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(2*y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(2*y.size(),2*y.size());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(2*y.size());
  ci2.fill(C);
  //Append RHS
  Eigen::VectorXd ci0(4.0*y.size());
  ci0 << -ci1, ci2;
  //Append LHS
  Eigen::MatrixXd CI(2*CI1.rows(), CI1.cols());
  CI <<  CI1,
        -CI1;

  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);

  //Create matrix Q
  Eigen::MatrixXd Q = Eigen::MatrixXd(2*y.size(),2*y.size());
  Q<< K,-K,
     -K, K;

  //Nearest positive semidefinite
  //nearPositiveDefinite(Q,1e-5);
  Q = nearPDefinite(Q, 1e+2, 1e-06, 1e-07, 1e-08, true);

  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,g, CE.transpose(),ce0, CI.transpose(), ci0);

  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms,
                            Rcpp::Named("Epsilon") = epsilon);
}

//' @name Predicted CSVRL1
//' @title C-SVR L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVR L1:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVRL1 List of Results of the CSVRL1
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param Xprev Numeric matrix with the explanatory variables (predicted). Dimension equal MxP
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return Eigen::VectorXd with the predicted values for Xpred
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
Eigen::VectorXd PredictedCSVRL1(Rcpp::List CSVRL1, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev){
  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVRL1["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVRL1["Kernel"]);

  //Get the parameters
  Eigen::RowVectorXd parms = as<Eigen::RowVectorXd> (CSVRL1["Parameters"]);

  //Get the epsilon band
  double epsilon = as<double> (CSVRL1["Epsilon"]);

  //Total number of observations
  int size = Xprev.rows();
  Eigen::VectorXd predVec(size);
  //Separating the SV
  Eigen::VectorXd diffLambda = SV.head(X.rows()) - SV.tail(X.rows());

  for(int i=0;i<size;i++){
    //Create the Kernel Matrix
    Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
    Eigen::VectorXd F = diffLambda.array() *K.array();
    predVec(i) = F.sum();
  }

  //Get the new bias term
  Eigen::VectorXd gamma = y.array()-predVec.array()-epsilon;
  double bGamma = gamma.mean();

  return(predVec.array()-bGamma);
}



//' @name Pseudo R2 - Predicted CSVRL1
//' @title C-SVR L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVR L1:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVRL1 List of Results of the CSVRL1
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return Eigen::VectorXd with the Pseudo R2 for each variable.
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
Eigen::VectorXd R2PredictedCSVRL1(Rcpp::List CSVRL1, Eigen::VectorXd y, Eigen::MatrixXd X){
  //Results
  Eigen::VectorXd R2vec(X.cols());
  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVRL1["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVRL1["Kernel"]);

  //Get the parameters
  Eigen::RowVectorXd parms = as<Eigen::RowVectorXd> (CSVRL1["Parameters"]);

  //Get the epsilon band
  double epsilon = as<double> (CSVRL1["Parameters"]);

  //Prediction for the full model
  int size = X.rows();
  Eigen::VectorXd predVec(size);
  predVec = PredictedCSVRL1(CSVRL1, y,  X, X);
  //Sum of squared errors
  double SSE = predVec.squaredNorm();

  //For each variable:
  for(int v=0;v<X.cols();v++){
    //Zero columns
    Eigen::MatrixXd Xprev = X;
    //Zero the variable
    Xprev.col(v).fill(0.0);
    //Total number of observations
    int size = Xprev.rows();
    Eigen::VectorXd predVec2(size);
    predVec2 = PredictedCSVRL1(CSVRL1, y,  X, Xprev);
    double SSEvar = predVec2.squaredNorm();
    R2vec(v) = SSEvar/SSE;
  }
  return(R2vec);
}


//' @param y Vector with dependent variables. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param epsilon Insentitive band. Should be epsilon>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.

// [[Rcpp::export]]
Eigen::MatrixXd minimumCSVRL1(Eigen::VectorXd y, Eigen::MatrixXd X, double epsilon, std::string kernel, Eigen::RowVectorXd parms){
  Eigen::VectorXd SV(2*y.size());
  //Create the one vector 2Nx1
  Eigen::VectorXd yfull = Eigen::VectorXd(2*y.size());
  yfull<< (-1.0)*y, (+1.0)*y;
  Eigen::VectorXd evec = Eigen::VectorXd(2*y.size());
  evec.fill(epsilon);
  Eigen::VectorXd g = Eigen::VectorXd(2*y.size());
  g = evec+yfull;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix Q
  Eigen::MatrixXd Q = Eigen::MatrixXd(2*y.size(),2*y.size());
  Q<< K,-K,
     -K, K;
  //Nearest positive semidefinite matrix in terms of Frobenius norm

//  nearPositiveDefinite(Q,1e-10);
//  std::cout<<Q<<std::endl;
  //Calculate the Minimum Objective function
  SV = -Q.inverse()*g;
  return(Q);
}
