#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "KernelComputation.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
#include <functional>

using namespace Rcpp;

//TODO: Reescrever os kernels usando Eigen ao invés de arma e remover esse casting da linha 32
//Typecasting between Eigen::MatrixXd datMat and arma::mat
arma::mat convertEigenToArma(Eigen::MatrixXd datMat){
  arma::mat res(datMat.rows(), datMat.cols());
  for(int i=0;i<datMat.rows();i++){
    for(int j=0;j<datMat.cols();j++){
      res(i,j)=datMat(i,j);
    }
  }
  return(res);
}

//Typecasting between Eigen::MatrixXd datMat and arma::mat
arma::vec convertEigenToArma(Eigen::RowVectorXd datMat){
  arma::vec res(datMat.size());
  for(int i=0;i<datMat.size();i++){
      res(i)=datMat(i);
  }
  return(res);
}

// Create the Kernel matrix
// @param datMat  Matrix with the data
// @param function Kernel Function
// @param parms vector of parameters fot the kernel
// @return Kernel Matrix
Eigen::VectorXd KernelMatrix(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat,const std::function<double (arma::vec, arma::vec, arma::vec)> kernel, arma::vec parms){
  //Get the number of rows
  int rows=datMat.rows();
  //Typecasting
  arma::mat datMat2 = convertEigenToArma(datMat);
  arma::vec datPred2 = convertEigenToArma(predMat);
  //Initialize the matriz
  Eigen::VectorXd matKernel = Eigen::VectorXd::Zero(rows);
  for(unsigned int c1=0;c1<rows;c1++){
      //First column with variables
      arma::vec vec1 = datMat2.row(c1).t();
      //Calculate the kernel value
      double val= kernel(vec1,datPred2,parms);
      //Store the kernel value
      matKernel(c1)=val;
  }

  return(matKernel);
}

// Create the Kernel matrix
// @param datMat  Matrix with the data
// @param function Kernel Function
// @param parms vector of parameters fot the kernel
// @return Kernel Matrix
Eigen::MatrixXd KernelMatrix(Eigen::MatrixXd datMat,const std::function<double (arma::vec, arma::vec, arma::vec)> kernel, arma::vec parms){
  //Get the number of rows
  int rows=datMat.rows();
  //Typecasting
  arma::mat datMat2 = convertEigenToArma(datMat);
  //Initialize the matriz
  Eigen::MatrixXd matKernel = Eigen::MatrixXd::Zero(rows,rows);
  for(unsigned int c1=0;c1<rows;c1++){
    for(unsigned int c2=c1;c2<rows;c2++){
      //First column with variables
      arma::vec vec1 = datMat2.row(c1).t();
      //Second column with variables
      arma::vec vec2 = datMat2.row(c2).t();
      //Calculate the kernel value
      double val= kernel(vec1,vec2,parms);
      //Store the kernel value
      matKernel(c1,c2)=matKernel(c2,c1)=val;
    }
  }

  return(matKernel);
}

//Define the Kernel functions
double CauchyKernel(arma::vec x,arma::vec y,arma::vec parms);
double PolynomialKernel(arma::vec x,arma::vec y,arma::vec parms);
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
double HistogramIntersectionKernel(arma::vec x,arma::vec y,arma::vec parms);
double WaveletKernel(arma::vec x,arma::vec y,arma::vec parms);
double MexicanHatKernel(arma::vec x,arma::vec y,arma::vec parms);
double MorletKernel(arma::vec x,arma::vec y,arma::vec parms);
double GeneralizedHistogramIntersectionKernel(arma::vec x,arma::vec y,arma::vec parms);
double CircularKernel(arma::vec x,arma::vec y,arma::vec parms);
double SphericalKernel(arma::vec x,arma::vec y,arma::vec parms);
double LogKernel(arma::vec x,arma::vec y,arma::vec parms);
double WaveKernel(arma::vec x,arma::vec y,arma::vec parms);
double HellingerKernel(arma::vec x,arma::vec y,arma::vec parms);
double DirichletKernel(arma::vec x,arma::vec y,arma::vec parms);
double PearsonKernel(arma::vec x,arma::vec y,arma::vec parms);
double SigmoidKernel(arma::vec x,arma::vec y,arma::vec parms);
double SquaredSincKernel(arma::vec x,arma::vec y,arma::vec parms);
double SymmetricTriangleKernel(arma::vec x,arma::vec y,arma::vec parms);
double ThinSplinePlateKernel(arma::vec x,arma::vec y,arma::vec parms);
double ANOVAKernel(arma::vec x,arma::vec y,arma::vec parms);
double SplineKernel(arma::vec x,arma::vec y,arma::vec parms);
double BesselKernel(arma::vec x,arma::vec y,arma::vec parms);
double ArccosKernel(arma::vec x,arma::vec y,arma::vec parms);

//TODO: Testar se o Kernel construido eh semi-positivo definido e caso nao seja, trocar os autovalores negativos por zero.
//TODO: http://scicomp.stackexchange.com/questions/10450/how-to-implement-the-spectral-decomposition-of-a-symmetric-dense-matrix-via-eige
//TOOD: https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html
// [[Rcpp::export]]
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,std::string stringValue, arma::vec parms){
  //Get the number of columns
  int rows=datMat.rows();
  //Initialize the matrix
  Eigen::MatrixXd Kernel = Eigen::MatrixXd::Zero(rows,rows);

  if(stringValue=="Cauchy"){
    Kernel = KernelMatrix(datMat, CauchyKernel ,parms);
  }
  else if(stringValue=="Chi-Square"){
    Kernel = KernelMatrix(datMat, ChiSquareKernel ,parms);
  }
  else if(stringValue=="Exponential"){
    Kernel = KernelMatrix(datMat, ExponentialKernel ,parms);
  }
  else if(stringValue=="Gaussian"){
    Kernel = KernelMatrix(datMat, GaussianKernel ,parms);
  }
  else if(stringValue=="Generalized T-Student"){
    Kernel = KernelMatrix(datMat, GeneralizedTStudentKernel ,parms);
  }
  else if(stringValue=="Hyperbolic Tangent"){
    Kernel = KernelMatrix(datMat, HyperbolicTangentKernel ,parms);
  }
  else if(stringValue=="Inverse Multiquadratic"){
    Kernel = KernelMatrix(datMat, InverseMultiquadraticKernel ,parms);
  }
  else if(stringValue=="Laplacian"){
    Kernel = KernelMatrix(datMat, LaplacianoKernel ,parms);
  }
  else if(stringValue=="Linear"){
    Kernel = KernelMatrix(datMat, LinearKernel ,parms);
  }
  else if(stringValue=="Log-Linear"){
    Kernel = KernelMatrix(datMat, LogLinearKernel ,parms);
  }
  else if(stringValue=="Polynomial"){
    Kernel = KernelMatrix(datMat, PolynomialKernel ,parms);
  }
  else if(stringValue=="Multiquadratic"){
    Kernel = KernelMatrix(datMat, MultiquadraticKernel ,parms);
  }
  else if(stringValue=="Power"){
    Kernel = KernelMatrix(datMat, PowerKernel ,parms);
  }
  else if(stringValue=="Rational Quadratic"){
    Kernel = KernelMatrix(datMat, RationalQuadraticKernel ,parms);
  }
  else if(stringValue=="Wavelet"){
    Kernel = KernelMatrix(datMat, WaveletKernel ,parms);
  }
  else if(stringValue=="Histogram Intersection"){
    Kernel = KernelMatrix(datMat, HistogramIntersectionKernel ,parms);
  }
  else if(stringValue=="Mexican-Hat"){
    Kernel = KernelMatrix(datMat, MexicanHatKernel ,parms);
  }
  else if(stringValue=="Morlet"){
    Kernel = KernelMatrix(datMat, MorletKernel ,parms);
  }
  else if(stringValue=="Generalized Histogram Intersection"){
    Kernel = KernelMatrix(datMat, GeneralizedHistogramIntersectionKernel ,parms);
  }
  else if(stringValue=="Circular"){
    Kernel = KernelMatrix(datMat, CircularKernel ,parms);
  }
  else if(stringValue=="Spherical"){
    Kernel = KernelMatrix(datMat, SphericalKernel ,parms);
  }
  else if(stringValue=="Log-Kernel"){
    Kernel = KernelMatrix(datMat, LogKernel ,parms);
  }
  else if(stringValue=="Wave"){
    Kernel = KernelMatrix(datMat, WaveKernel ,parms);
  }
  else if(stringValue=="Hellinger"){
    Kernel = KernelMatrix(datMat, HellingerKernel ,parms);
  }
  else if(stringValue=="Dirichlet"){
    Kernel = KernelMatrix(datMat, DirichletKernel ,parms);
  }
  else if(stringValue=="Pearson"){
    Kernel = KernelMatrix(datMat, PearsonKernel ,parms);
  }
  else if(stringValue=="Sigmoid"){
    Kernel = KernelMatrix(datMat, SigmoidKernel ,parms);
  }
  else if(stringValue=="Symmetric Triangle"){
    Kernel = KernelMatrix(datMat, SymmetricTriangleKernel ,parms);
  }
  else if(stringValue=="Thin Spline Plate"){
    Kernel = KernelMatrix(datMat, ThinSplinePlateKernel ,parms);
  }
  else if(stringValue=="ANOVA"){
    Kernel = KernelMatrix(datMat, ANOVAKernel ,parms);
  }
  else if(stringValue=="Spline"){
    Kernel = KernelMatrix(datMat, SplineKernel ,parms);
  }
  else if(stringValue=="Bessel"){
    Kernel = KernelMatrix(datMat, BesselKernel ,parms);
  }
  else if(stringValue=="Arccos"){
    Kernel = KernelMatrix(datMat, ArccosKernel ,parms);
  }
  return(Kernel);
}


// [[Rcpp::export]]
Eigen::MatrixXd KernelMatrixComputationPred(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat,std::string stringValue, arma::vec parms){
  //Get the number of columns
  int rows=datMat.rows();
  //Initialize the matrix
  Eigen::VectorXd Kernel = Eigen::VectorXd::Zero(rows);


  if(stringValue=="Cauchy"){
    Kernel = KernelMatrix(datMat, predMat, CauchyKernel ,parms);
  }
  else if(stringValue=="Chi-Square"){
    Kernel = KernelMatrix(datMat, predMat, ChiSquareKernel ,parms);
  }
  else if(stringValue=="Exponential"){
    Kernel = KernelMatrix(datMat, predMat, ExponentialKernel ,parms);
  }
  else if(stringValue=="Gaussian"){
    Kernel = KernelMatrix(datMat, predMat, GaussianKernel ,parms);
  }
  else if(stringValue=="Generalized T-Student"){
    Kernel = KernelMatrix(datMat, predMat, GeneralizedTStudentKernel ,parms);
  }
  else if(stringValue=="Hyperbolic Tangent"){
    Kernel = KernelMatrix(datMat, predMat, HyperbolicTangentKernel ,parms);
  }
  else if(stringValue=="Inverse Multiquadratic"){
    Kernel = KernelMatrix(datMat, predMat, InverseMultiquadraticKernel ,parms);
  }
  else if(stringValue=="Laplacian"){
    Kernel = KernelMatrix(datMat, predMat, LaplacianoKernel ,parms);
  }
  else if(stringValue=="Linear"){
    Kernel = KernelMatrix(datMat, predMat, LinearKernel ,parms);
  }
  else if(stringValue=="Log-Linear"){
    Kernel = KernelMatrix(datMat, predMat, LogLinearKernel ,parms);
  }
  else if(stringValue=="Polynomial"){
    Kernel = KernelMatrix(datMat, predMat, PolynomialKernel ,parms);
  }
  else if(stringValue=="Multiquadratic"){
    Kernel = KernelMatrix(datMat, predMat, MultiquadraticKernel ,parms);
  }
  else if(stringValue=="Power"){
    Kernel = KernelMatrix(datMat, predMat, PowerKernel ,parms);
  }
  else if(stringValue=="Rational Quadratic"){
    Kernel = KernelMatrix(datMat, predMat, RationalQuadraticKernel ,parms);
  }
  else if(stringValue=="Wavelet"){
    Kernel = KernelMatrix(datMat, predMat, WaveletKernel ,parms);
  }
  else if(stringValue=="Histogram Intersection"){
    Kernel = KernelMatrix(datMat, predMat, HistogramIntersectionKernel ,parms);
  }
  else if(stringValue=="Mexican-Hat"){
    Kernel = KernelMatrix(datMat, predMat, MexicanHatKernel ,parms);
  }
  else if(stringValue=="Morlet"){
    Kernel = KernelMatrix(datMat, predMat, MorletKernel ,parms);
  }
  else if(stringValue=="Generalized Histogram Intersection"){
    Kernel = KernelMatrix(datMat, predMat, GeneralizedHistogramIntersectionKernel ,parms);
  }
  else if(stringValue=="Circular"){
    Kernel = KernelMatrix(datMat, predMat, CircularKernel ,parms);
  }
  else if(stringValue=="Spherical"){
    Kernel = KernelMatrix(datMat, predMat, SphericalKernel ,parms);
  }
  else if(stringValue=="Log-Kernel"){
    Kernel = KernelMatrix(datMat, predMat, LogKernel ,parms);
  }
  else if(stringValue=="Wave"){
    Kernel = KernelMatrix(datMat, predMat, WaveKernel ,parms);
  }
  else if(stringValue=="Hellinger"){
    Kernel = KernelMatrix(datMat, predMat, HellingerKernel ,parms);
  }
  else if(stringValue=="Dirichlet"){
    Kernel = KernelMatrix(datMat, predMat, DirichletKernel ,parms);
  }
  else if(stringValue=="Pearson"){
    Kernel = KernelMatrix(datMat, predMat, PearsonKernel ,parms);
  }
  else if(stringValue=="Sigmoid"){
    Kernel = KernelMatrix(datMat, predMat, SigmoidKernel ,parms);
  }
  else if(stringValue=="Symmetric Triangle"){
    Kernel = KernelMatrix(datMat, predMat, SymmetricTriangleKernel ,parms);
  }
  else if(stringValue=="Thin Spline Plate"){
    Kernel = KernelMatrix(datMat, predMat, ThinSplinePlateKernel ,parms);
  }
  else if(stringValue=="ANOVA"){
    Kernel = KernelMatrix(datMat, predMat, ANOVAKernel ,parms);
  }
  else if(stringValue=="Spline"){
    Kernel = KernelMatrix(datMat, predMat, SplineKernel ,parms);
  }
  else if(stringValue=="Bessel"){
    Kernel = KernelMatrix(datMat, predMat, BesselKernel ,parms);
  }
  else if(stringValue=="Arccos"){
    Kernel = KernelMatrix(datMat, predMat, ArccosKernel ,parms);
  }

  return(Kernel);
}

