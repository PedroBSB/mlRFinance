#include <RcppEigen.h>
#include "SVR.h"
#include "Utils.h"
#include "ErrorMeasure.h"
#include "progress.hpp"
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/
//Print Object at Console
void PrintObject(Eigen::MatrixXd mat);
void PrintObject(Eigen::VectorXd vec);

//C-SVR L1
Rcpp::List CSVRL1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, double epsilon, std::string kernel, Eigen::RowVectorXd parms);
Eigen::VectorXd PredictedCSVRL1(Rcpp::List CSVRL1, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev);

//Error Measure
Rcpp::List ErrorMeasures(Eigen::VectorXd y, Eigen::VectorXd yPred);


/***********************************************************************************************/
/*********************************    VOLATILITY PROXY FUNCTIONS  ******************************/
/***********************************************************************************************/

//Regular htt (Brooks, 2001; Brooks e Persand, 2003; Chen et al., 2010)
Eigen::VectorXd httFunction(Eigen::VectorXd y){
  Eigen::VectorXd htt = (y.array()-y.mean()).pow(2);
  return(htt);
}


/***********************************************************************************************/
/*********************************     GARCH SVC FUNCTIONS    **********************************/
/***********************************************************************************************/

//  TODO: Implementar outras proxies para a volatilidade Andersen e Bollerslev (1998)
//' @name Garch C-SVR L1
//' @title Garch C-SVR L1 - Garch Support Vector Regression with C cost and L1 regularization.
//' @description Training and Forecasting volatility
//'
//' @param train Vector of returns (training set). Dimension equal Nx1.
//' @param valid Vector of returns (validation set). Dimension equal Mx1.
//' @param C Cost parameter. Should be C>0.
//' @param epsilon Insentitive band. Should be epsilon>0.
//' @param kernel Name of the kernel that will be used for the mean equation.
//' @param parms Parameters associated with chosen kenel for the mean equation.
//' @param kernel Name of the kernel that will be used for the mean equation.
//' @param parms Parameters associated with chosen kenel for the mean equation.
//' @return List Support Vectors Mean, Forecast mean, EAMmean,
//' Support Vectors Volatility, Forecast volatility, EAMvolat, .
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML1(d, A, 1, 0.1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List GARCHCSVRL1(Eigen::VectorXd train, Eigen::VectorXd valid, double Cmean, double epsilonMean,
                       double Cvola, double epsilonVola,
                       std::string kernelMean, Eigen::RowVectorXd parmsMean,
                       std::string kernelVolat, Eigen::RowVectorXd parmsVola){

  //Size of the Time Series
  int sizeTrain = train.size();
  int sizeValid = valid.size();

  //Create the AR form
  Eigen::VectorXd y = train.tail(sizeTrain-1);
  Eigen::MatrixXd X(sizeTrain-1,1);
  X.col(0)= train.head(sizeTrain-1).array();

  //Create the AR form
  Eigen::VectorXd yValid = valid.tail(sizeValid-1);
  Eigen::MatrixXd Xvalid(sizeValid-1,1);
  Xvalid.col(0)= valid.head(sizeValid-1).array();

  //Step 1: Training and validating MEAN
  Rcpp::List SVRmean = CSVRL1(y, X, Cmean, epsilonMean, kernelMean, parmsMean);

  //Forecasting the results
  Eigen::VectorXd yPred = PredictedCSVRL1(SVRmean, y, X, X);
  Eigen::VectorXd yValidPred = PredictedCSVRL1(SVRmean, y, X, Xvalid);

  //Calculate the error measure
  Rcpp::List yPredError = ErrorMeasures(y,yPred);
  Rcpp::List yValidPredError = ErrorMeasures(yValid,yValidPred);

  //Calculate the squared residual
  Eigen::VectorXd residTrain2 = (y-yPred).array().pow(2);
  Eigen::VectorXd residValid2 = (yValid-yValidPred).array().pow(2);

  //Calculate the volatilility proxy
  Eigen::VectorXd httTrain;
  Eigen::VectorXd httValid;
  if(true){
    //TODO: Implement others proxies to htt
    httTrain = httFunction(y);
    httValid = httFunction(yValid);
  }

  //Step 2: Training and validating Volatility
  //Create the GARCH form (Tarining)
  Eigen::VectorXd httY = httTrain.tail(y.size()-1);
  Eigen::MatrixXd XGarchTrain(y.size()-1,2);
  XGarchTrain.col(0)= residTrain2.head(y.size()-1).array();
  XGarchTrain.col(1)= httTrain.head(y.size()-1).array();

  //Create the GARCH form (Validation)
  Eigen::VectorXd httYvalid = httValid.tail(yValid.size()-1);
  Eigen::MatrixXd XGarchValid(yValid.size()-1,2);
  XGarchValid.col(0)= residValid2.head(yValid.size()-1).array();
  XGarchValid.col(1)= httValid.head(yValid.size()-1).array();

  //Training the Garch-SVR
  Rcpp::List SVRgarch = CSVRL1(httY, XGarchTrain, Cvola, epsilonVola, kernelVolat, parmsVola);

  //Forecasting the results
  Eigen::VectorXd garchPred = PredictedCSVRL1(SVRgarch, httY, XGarchTrain, XGarchTrain);
  Eigen::VectorXd garchValidPred = PredictedCSVRL1(SVRgarch, httY, XGarchTrain, XGarchValid);

  //Calculate the error measure
  Rcpp::List garchPredError = ErrorMeasures(httY,garchPred);
  Rcpp::List garchValidPredError = ErrorMeasures(httYvalid,garchValidPred);

  //Insert NA in the first position Mean Vectors
  Eigen::VectorXd NA(1);
  NA(0) = std::numeric_limits<double>::quiet_NaN();

  Eigen::VectorXd yPredNew(yPred.size()+1);
  yPredNew << NA, yPred;
  Eigen::VectorXd yValidPredNew(yValidPred.size()+1);
  yValidPredNew << NA, yValidPred;

  //Insert NA in the first position Volaitility Vectors
  Eigen::VectorXd garchPredNew(garchPred.size()+1);
  garchPredNew << NA, garchPred;
  Eigen::VectorXd garchValidPredNew(garchValidPred.size()+1);
  garchValidPredNew << NA, garchValidPred;

  //Return the results
  return Rcpp::List::create(Rcpp::Named("PredictedMeanTraining") = yPredNew,
                            Rcpp::Named("PredictedMeanValidation") = yValidPredNew,
                            Rcpp::Named("ErrorMeasureTraining") = yPredError,
                            Rcpp::Named("ErrorMeasureValidation") = yValidPredError,
                            Rcpp::Named("PredictedGarchTraining") = garchPredNew,
                            Rcpp::Named("PredictedGarchValidation") = garchValidPredNew,
                            Rcpp::Named("ErrorMeasureTrainingGarch") = garchPredError,
                            Rcpp::Named("ErrorMeasureValidationGarch") = garchValidPredError);
}
























