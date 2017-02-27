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

//C-SVM L1
Rcpp::List CSVML1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, Eigen::RowVectorXd parms, bool biasTerm);
Eigen::VectorXd PredictedCSVML1(Rcpp::List CSVML1, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev, int typePredict, bool biasTerm);

//C-SVR L1
Rcpp::List CSVRL1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, double epsilon, std::string kernel, Eigen::RowVectorXd parms);
Eigen::VectorXd PredictedCSVRL1(Rcpp::List CSVRL1, Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev);

//Error Measure
Rcpp::List ErrorMeasures(Eigen::VectorXd y, Eigen::VectorXd yPred);
Rcpp::List ErrorMeasuresBinary(Eigen::VectorXd y, Eigen::VectorXd yPred);




/***********************************************************************************************/
/*********************************     PORTFOLIO SELECTION SVM FUNCTIONS  **********************/
/***********************************************************************************************/

//' @name Portfolio Selection C-SVR L1
//' @title Portfolio Selection C-SVR L1 - Portfolio Selection Support Vector Regression
//' with C cost and L1 regularization.
//' @description Training and Forecasting the portfolio
//'
//' @param y_train Binary Vector of good (+1) and bad (-1) firms. Dimension equal Nx1.
//' @param X_train Fundamental matrix explaining the y_train, y_train_t=f(X_train_(t-1)). Dimension equal NxP.
//' @param y_valid Binary Vector of good (+1) and bad (-1) firms. Dimension equal Mx1.
//' @param X_valid Fundamental matrix explaining the y_train, y_train_t=f(X_train_(t-1)). Dimension equal MxP.
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used for the mean equation.
//' @param parms Parameters associated with chosen kenel for the mean equation.
//' @param typePredict 0-Binary(-1,+1), 1-Probability, 2- Raw result
//' @return List Support Vectors, Forecast , EAM
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
Rcpp::List PortfolioSelectionCSVML1(Eigen::VectorXd y_train, Eigen::MatrixXd X_train,
                                    Eigen::VectorXd y_valid, Eigen::MatrixXd X_valid,
                                    double C,
                                    std::string kernel, Eigen::RowVectorXd parms,int typePredict){

  //Step 1: Training and validating MEAN
  Rcpp::List SVMport = CSVML1(y_train, X_train, C, kernel, parms, true);

  //Forecasting the results
  Eigen::VectorXd yPred = PredictedCSVML1(SVMport,y_train, X_train, X_train, typePredict, true);
  Eigen::VectorXd yValidPred = PredictedCSVML1(SVMport,y_train, X_train, X_valid,typePredict, true);

  //Calculate the error measure
  Rcpp::List yPredError;
  Rcpp::List yValidPredError;
  if(typePredict==0){
   yPredError = ErrorMeasuresBinary(y_train,yPred);
   yValidPredError = ErrorMeasuresBinary(y_valid,yValidPred);
  }
  else if(typePredict==1){
    yPredError = ErrorMeasures(y_train,yPred);
    yValidPredError = ErrorMeasures(y_valid,yValidPred);

  }
  else{
    yPredError = ErrorMeasures(y_train,yPred);
    yValidPredError = ErrorMeasures(y_valid,yValidPred);
  }

  //Return the results
  return Rcpp::List::create(Rcpp::Named("PredictedTraining") = yPred,
                            Rcpp::Named("PredictedValidation") = yValidPred,
                            Rcpp::Named("ErrorMeasureTraining") = yPredError,
                            Rcpp::Named("ErrorMeasureValidation") = yValidPredError);
}



// [[Rcpp::export]]
Rcpp::List PortfolioSelectionCSVRL1(Eigen::VectorXd y_train, Eigen::MatrixXd X_train,
                                    Eigen::VectorXd y_valid, Eigen::MatrixXd X_valid,
                                    double C, double epsilon,
                                    std::string kernel, Eigen::RowVectorXd parms){

  //Step 1: Training and validating
  Rcpp::List SVRport = CSVRL1(y_train, X_train, C, epsilon, kernel, parms);

  //Forecasting the results
  Eigen::VectorXd yPred = PredictedCSVRL1(SVRport,y_train, X_train, X_train);
  Eigen::VectorXd yValidPred = PredictedCSVRL1(SVRport,y_train, X_train, X_valid);

  //Calculate the error measure
  Rcpp::List yPredError;
  Rcpp::List yValidPredError;
  yPredError = ErrorMeasures(y_train,yPred);
  yValidPredError = ErrorMeasures(y_valid,yValidPred);

  //Return the results
  return Rcpp::List::create(Rcpp::Named("PredictedTraining") = yPred,
                            Rcpp::Named("PredictedValidation") = yValidPred,
                            Rcpp::Named("ErrorMeasureTraining") = yPredError,
                            Rcpp::Named("ErrorMeasureValidation") = yValidPredError);
}



















