#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "Utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
using namespace Rcpp;

//Mean Square Error
double MSEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  double res = (y-yPred).squaredNorm();
  return(res);
}

// [[Rcpp::export]]
Rcpp::List ErrorMeasures(Eigen::VectorXd y, Eigen::VectorXd yPred){
  //Calculate Mean Square Error
  double mse = MSEfunction(y, yPred);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("MSE") = mse);
}
