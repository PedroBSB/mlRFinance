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


// [[Rcpp::export]]
NumericVector cpp_order(const NumericVector & x, bool desc = false) {
  auto n = x.size();
  NumericVector idx = no_init(n);
  auto begin = std::begin(idx), end = std::end(idx);
  std::iota(begin, end, static_cast<size_t>(0));
  auto comparator = [&](const size_t & a, const size_t & b){ return x[a] < x[b]; };
  std::sort(begin, end, comparator);
  // In R indicies start at 1
  auto plus_1 = [](const size_t & v){ return v + 1; };
  std::transform(begin, end, begin, plus_1);
  if (desc) std::reverse(begin, end);
  return idx;
}


/***********************************************************************************************/
/*********************************     SVR FUNCTIONS    ****************************************/
/***********************************************************************************************/

//' @name CSVWQR
//' @title C-SVWQR - Support Vector Weighted Quantile Regression with C cost and L1 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVWQR:
//'
//' Xu, Q., Zhang, J., Jiang, C., Huang, X., & He, Y. (2015).
//' Weighted quantile regression via support vector machine.
//' Expert Systems with Applications, 42(13), 5441-5451.
//'
//'
//' Min (1/2)u^{t}Qu+g^{t}u
//' s.t.
//' 0<=u<= tau*C*qi
//' 0<=u^{*}<= (1-tau)*C*qi
//'
//' where u=(lambda*,lambda), g=(e-y,e+y)
//' and Q=|K -K|.
//'       |-K K|
//' C is the Cost parameter and e (epsilon) is the insentitive band
//'
//' @param y Vector with dependent variables. Dimension equal Nx1.
//' @param rank Vector with the rank of the dependent variables. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param tau Quantile of interest. Should be 0<tau<1.
//' @gamma gamma Weight smooth based on Cao and Gu (2002),
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
Rcpp::List CSVWQR(Eigen::VectorXd y, Eigen::MatrixXd X, double C, double tau, double gamma, std::string kernel, Eigen::RowVectorXd parms){
  //Support Vectors
  Eigen::VectorXd SV(2*y.size());
  //Create the one vector 2Nx1
  Eigen::VectorXd g = Eigen::VectorXd(2*y.size());
  g<< (-1.0)*y, (+1.0)*y;
  //RHS equality
  Eigen::VectorXd ce0;
  //LHS equality
  Eigen::MatrixXd CE;
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(2*y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(2*y.size(),2*y.size());

  //Find the rank
  Rcpp::NumericVector temp = as<Rcpp::NumericVector>(wrap(y));
  NumericVector rankTemp = cpp_order(temp,true);
  Eigen::VectorXd rank = as<Eigen::VectorXd>(wrap(rankTemp));

  //Calculate the q according  Cao and Gu (2002)
  rank = gamma-2.0*rank.array()/rank.size();
  rank = rank.array().exp();
  Eigen::VectorXd q = 2.0/(1.0+rank.array());

  //RHS: Inequality 2 tau*C*qi
  Eigen::VectorXd ci21(y.size());
  ci21 = tau*C*q.array();

  //RHS: Inequality 2 (1-tau)*C*qi
  Eigen::VectorXd ci22(y.size());
  ci22 = (1.0-tau)*C*q.array();

  //RHS: Inequality
  Eigen::VectorXd ci2(2*y.size());
  ci2 << ci21, ci22;

  //Append RHS
  Eigen::VectorXd ci0(4.0*y.size());
  ci0 << ci1, ci2;
  //Append LHS
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
  //Diagonal matrix
  Eigen::VectorXd me(2*y.size());
  me.fill(-1.0);
  Eigen::MatrixXd mI = me.asDiagonal();
  //Vertical concatenation
  CI << CI1,
        mI;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix Q
  Eigen::MatrixXd Q = Eigen::MatrixXd(2*y.size(),2*y.size());
  Q<< K,-K,
     -K, K;
  //Nearest positive semidefinite
  //nearPositiveDefinite(Q,1e-10);
  Q = nearPDefinite(Q, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,g, CE.transpose(),ce0, CI.transpose(), ci0);

  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//' @name Predicted CSVWQR
//' @title C-SVWQR - Support Vector Weighted Quantile Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVWQR:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVWQR List of Results of the CSVWQR
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
Eigen::VectorXd PredictedCSVWQR(Rcpp::List CSVWQR, Eigen::MatrixXd X, Eigen::MatrixXd Xprev){
  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVWQR["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVWQR["Kernel"]);

  //Get the parameters
  Eigen::RowVectorXd parms = as<Eigen::RowVectorXd> (CSVWQR["Parameters"]);

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
return(predVec);
}



//' @name Pseudo R2 - Predicted CSVWQR
//' @title C-SVWQR - Support Vector Weighted Quantile Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVWQR:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVWQR List of Results of the CSVWQR
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
Eigen::VectorXd R2PredictedCSVWQR(Rcpp::List CSVWQR, Eigen::MatrixXd X){
  //Results
  Eigen::VectorXd R2vec(X.cols());
  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVWQR["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVWQR["Kernel"]);

  //Get the parameters
  Eigen::RowVectorXd parms = as<Eigen::RowVectorXd> (CSVWQR["Parameters"]);

  //Prediction for the full model
  //Total number of observations
  int size = X.rows();
  Eigen::VectorXd predVec(size);
  //Separating the SV
  Eigen::VectorXd diffLambda = SV.head(X.rows()) - SV.tail(X.rows());

  for(int i=0;i<size;i++){
    //Create the Kernel Matrix
    Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
    Eigen::VectorXd F = diffLambda.array() *K.array();
    predVec(i) = F.sum();
  }
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
    //Separating the SV
    Eigen::VectorXd diffLambda = SV.head(X.rows()) - SV.tail(X.rows());

    for(int i=0;i<size;i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
      Eigen::VectorXd F = diffLambda.array() *K.array();
      predVec2(i) = F.sum();
    }
    double SSEvar = predVec2.squaredNorm();
    R2vec(v) = SSEvar/SSE;
  }
  return(R2vec);
}
