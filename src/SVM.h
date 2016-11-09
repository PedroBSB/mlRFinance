// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef SVM_H
#define SVM_H

Rcpp::List CSVML1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, arma::vec parms);
Eigen::VectorXd PredictedCSVML1(Rcpp::List CSVML1,Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev, std::string kernel, arma::vec parms, int typePredict);
// This is the end of the header guard
#endif


