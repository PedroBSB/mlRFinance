// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef KERNELMATRIX_H
#define KERNELMATRIX_H

Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,std::string stringValue, Eigen::RowVectorXd parms);

Eigen::MatrixXd KernelMatrixComputationPred(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat,std::string stringValue, Eigen::RowVectorXd parms);

double KernelMatrixComputationValue(Eigen::RowVectorXd datMat, Eigen::RowVectorXd predMat, std::string stringValue, Eigen::RowVectorXd parms);

// This is the end of the header guard
#endif


