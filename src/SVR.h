// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef SVR_H
#define SVR_H

Rcpp::List CSVRL1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, double epsilon, std::string kernel, Eigen::VectorXd parms);
Eigen::VectorXd PredictedCSVRL1(Rcpp::List CSVRL1, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev);
// This is the end of the header guard
#endif


