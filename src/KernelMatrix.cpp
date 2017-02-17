#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "KernelComputation.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
#include <functional>

using namespace Rcpp;

// Create the Kernel matrix
// @param datMat  Matrix with the data
// @param function Kernel Function
// @param parms vector of parameters fot the kernel
// @return Kernel Matrix
Eigen::VectorXd KernelMatrix(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat, const std::function<double (Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::RowVectorXd)> kernel, Eigen::RowVectorXd parms){
  //Get the number of rows
  int rows = datMat.rows();
  //Initialize the matriz
  Eigen::VectorXd matKernel = Eigen::VectorXd::Zero(rows);
  for(unsigned int c1=0;c1<rows;c1++){
      //First column with variables
      Eigen::RowVectorXd vec1 = datMat.row(c1);
      //Calculate the kernel value
      double val= kernel(vec1,predMat,parms);
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
Eigen::MatrixXd KernelMatrix(Eigen::MatrixXd datMat,const std::function<double (Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::RowVectorXd )> kernel, Eigen::RowVectorXd parms){
  //Get the number of rows
  int rows=datMat.rows();
  //Initialize the matriz
  Eigen::MatrixXd matKernel = Eigen::MatrixXd::Zero(rows,rows);
  for(unsigned int c1=0;c1<rows;c1++){
    for(unsigned int c2=c1;c2<rows;c2++){
      //First column with variables
      Eigen::RowVectorXd vec1 = datMat.row(c1);
      //Second column with variables
      Eigen::RowVectorXd vec2 = datMat.row(c2);
      //Calculate the kernel value
      double val= kernel(vec1,vec2,parms);
      //Store the kernel value
      matKernel(c1,c2)=matKernel(c2,c1)=val;
    }
  }

  return(matKernel);
}

//Define the Kernel functions
double CauchyKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double PolynomialKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
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
double HistogramIntersectionKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double WaveletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double MexicanHatKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double MorletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double GeneralizedHistogramIntersectionKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double CircularKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double SphericalKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double LogKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double WaveKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double HellingerKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double DirichletKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double PearsonKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double SigmoidKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double SquaredSincKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double SymmetricTriangleKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double ThinSplinePlateKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double ANOVAKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double SplineKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double BesselKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);
double ArccosKernel(Eigen::RowVectorXd x,Eigen::RowVectorXd y,Eigen::RowVectorXd parms);

//TODO: Testar se o Kernel construido eh semi-positivo definido e caso nao seja, trocar os autovalores negativos por zero.
//TODO: http://scicomp.stackexchange.com/questions/10450/how-to-implement-the-spectral-decomposition-of-a-symmetric-dense-matrix-via-eige
//TOOD: https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html
// [[Rcpp::export]]
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,std::string stringValue, Eigen::RowVectorXd parms){
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
Eigen::MatrixXd KernelMatrixComputationPred(Eigen::MatrixXd datMat, Eigen::RowVectorXd predMat, std::string stringValue, Eigen::RowVectorXd parms){
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


// [[Rcpp::export]]
double KernelMatrixComputationValue(Eigen::RowVectorXd datMat, Eigen::RowVectorXd predMat, std::string stringValue, Eigen::RowVectorXd parms){

  //Initialize the matrix
  double Kernel = 0.0;

  if(stringValue=="Cauchy"){
    Kernel = CauchyKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Chi-Square"){
    Kernel = ChiSquareKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Exponential"){
    Kernel = ExponentialKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Gaussian"){
    Kernel = GaussianKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Generalized T-Student"){
    Kernel = GeneralizedTStudentKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Hyperbolic Tangent"){
    Kernel = HyperbolicTangentKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Inverse Multiquadratic"){
    Kernel = InverseMultiquadraticKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Laplacian"){
    Kernel = LaplacianoKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Linear"){
    Kernel = LinearKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Log-Linear"){
    Kernel = LogLinearKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Polynomial"){
    Kernel = PolynomialKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Multiquadratic"){
    Kernel = MultiquadraticKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Power"){
    Kernel = PowerKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Rational Quadratic"){
    Kernel = RationalQuadraticKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Wavelet"){
    Kernel = WaveletKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Histogram Intersection"){
    Kernel = HistogramIntersectionKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Mexican-Hat"){
    Kernel = MexicanHatKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Morlet"){
    Kernel = MorletKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Generalized Histogram Intersection"){
    Kernel = GeneralizedHistogramIntersectionKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Circular"){
    Kernel = CircularKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Spherical"){
    Kernel = SphericalKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Log-Kernel"){
    Kernel = LogKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Wave"){
    Kernel = WaveKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Hellinger"){
    Kernel = HellingerKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Dirichlet"){
    Kernel = DirichletKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Pearson"){
    Kernel = PearsonKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Sigmoid"){
    Kernel = SigmoidKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Symmetric Triangle"){
    Kernel = SymmetricTriangleKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Thin Spline Plate"){
    Kernel = ThinSplinePlateKernel(datMat, predMat, parms);
  }
  else if(stringValue=="ANOVA"){
    Kernel = ANOVAKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Spline"){
    Kernel = SplineKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Bessel"){
    Kernel = BesselKernel(datMat, predMat, parms);
  }
  else if(stringValue=="Arccos"){
    Kernel = ArccosKernel(datMat, predMat, parms);
  }

  return(Kernel);
}
