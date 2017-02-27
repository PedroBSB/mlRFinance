// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef UTILS_H
#define UTILS_H

bool IsPositiveDefinite(Eigen::MatrixXd mat);

void PrintTime();

void PrintObject(Eigen::MatrixXd mat);

void PrintObject(Eigen::VectorXd mat);

void nearPositiveDefinite(Eigen::MatrixXd &mat);

Eigen::MatrixXd nearPDefinite(Eigen::MatrixXd mat, int maxit, double eigtol, double conv_tol, double posd_tol, bool keepDiagonal);

void addNoise(Eigen::MatrixXd &mat,double noise);

// This is the end of the header guard
#endif


