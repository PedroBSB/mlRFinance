// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef KERNELCOMPUTATION_H
#define KERNELCOMPUTATION_H

double CauchyKernel(arma::vec x,arma::vec y,arma::vec parms);

double ChiSquareKernel(arma::vec x,arma::vec y,arma::vec parms);

double ExponentialKernel(arma::vec x,arma::vec y,arma::vec parms);

double GaussianKernel(arma::vec x,arma::vec y,arma::vec parms);

double GeneralizedTStudentKernel(arma::vec x,arma::vec y,arma::vec parms);

double HyperbolicTangentKernel(arma::vec x,arma::vec y,arma::vec parms);

double InverseMultiquadraticKernel(arma::vec x,arma::vec y,arma::vec parms);

double LaplacianoKernel(arma::vec x,arma::vec y,arma::vec parms);

double LinearKernel(arma::vec x,arma::vec y,arma::vec parms);

double LogLinearKernel(arma::vec x,arma::vec y,arma::vec parms);

double MultiquadraticKernel(arma::vec x,arma::vec y,arma::vec parms);

double PowerKernel(arma::vec x,arma::vec y,arma::vec parms);

double RationalQuadraticKernel(arma::vec x,arma::vec y,arma::vec parms);

double WaveletKernel(arma::vec x,arma::vec y,arma::vec parms);

double BesselKernel(arma::vec x,arma::vec y,arma::vec parms);

double ArccosKernel(arma::vec x,arma::vec y,arma::vec parms);

// This is the end of the header guard
#endif


