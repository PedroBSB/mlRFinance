// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H

Eigen::MatrixXd KernelMatrixComputation(arma::mat datMat,std::string stringValue, arma::vec parms);

Eigen::MatrixXd KernelMatrixComputationPred(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat,std::string stringValue, arma::vec parms);

// This is the end of the header guard
#endif


