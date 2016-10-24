// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef UTILS_H
#define UTILS_H

bool IsPositiveDefinite(Eigen::MatrixXd mat);

void PrintTime();

void PrintObject(arma::mat mat);

void PrintObjectLine(arma::uvec mat);

void PrintObjectLine(arma::vec mat);

void PrintObject(Eigen::MatrixXd mat);

void PrintObject(Eigen::VectorXd mat);

void nearPositiveDefinite(Eigen::MatrixXd &mat);

void addNoise(Eigen::MatrixXd &mat,double noise);

// This is the end of the header guard
#endif


