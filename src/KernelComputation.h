// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef KERNELCOMPUTATION_H
#define KERNELCOMPUTATION_H

double CauchyKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double ChiSquareKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double ExponentialKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double GaussianKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double GeneralizedTStudentKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double HyperbolicTangentKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double InverseMultiquadraticKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double LaplacianoKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double LinearKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double LogLinearKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double MultiquadraticKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double PowerKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double RationalQuadraticKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double WaveletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double BesselKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

double ArccosKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

// This is the end of the header guard
#endif


