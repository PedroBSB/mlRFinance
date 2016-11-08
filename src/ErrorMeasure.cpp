#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "Utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
using namespace Rcpp;

//Mean Square Error
double MSEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  double res = (y-yPred).squaredNorm();
  return(res);
}

// [[Rcpp::export]]
Rcpp::List ErrorMeasures(Eigen::VectorXd y, Eigen::VectorXd yPred){
  //Calculate Mean Square Error
  double mse = MSEfunction(y, yPred);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("MSE") = mse);
}


// [[Rcpp::export]]
Rcpp::List ErrorMeasuresBinary(Eigen::VectorXd y, Eigen::VectorXd yPred){
  //Source: https://en.wikipedia.org/wiki/Confusion_matrix
  //Get the size
  double n = (double) y.size();
  //Confusion matrix
  Eigen::MatrixXd CM = Eigen::MatrixXd::Zero(2,2);
  for(int i=0;i<n;i++){
    if(y(i)==+1.0 & yPred(i)==+1.0){
      CM(0,0)=CM(0,0)+1.0;
    }
    else if(y(i)==+1.0 & yPred(i)==-1.0){
      CM(0,1)=CM(0,1)+1.0;
    }
    else if(y(i)==-1.0 & yPred(i)==+1.0){
      CM(1,0)=CM(1,0)+1.0;
    }
    else{
      CM(1,1)=CM(1,1)+1.0;
    }
  }
  //Calculate Mean Square Error
  double dblPrevalence = (CM(0,0)+CM(0,1))/n;
  double dblTPR = CM(0,0)/(CM(0,0)+CM(0,1));
  double dblFNR = CM(0,1)/(CM(0,0)+CM(0,1));
  double dblFPR = CM(1,0)/(CM(0,0)+CM(0,1));
  double dblTNR = CM(1,1)/(CM(0,0)+CM(0,1));
  double dblACC = (CM(0,0)+CM(1,1))/n;
  double dblPPV = CM(0,0)/(CM(0,0)+CM(1,0));
  double dblFOR = CM(0,1)/(CM(0,0)+CM(1,0));
  double dblFDR = CM(1,0)/(CM(0,0)+CM(1,0));
  double dblNPV = CM(1,1)/(CM(0,0)+CM(1,0));
  double dblLRp = dblTPR/dblFPR;
  double dblLRn = dblFNR/dblTNR;
  double dblDOR = dblLRp/dblLRn;

  //Return the results
  return Rcpp::List::create(Rcpp::Named("Prevalence") = dblPrevalence,
                            Rcpp::Named("TruePositiveRate") = dblTPR,
                            Rcpp::Named("FalseNegativeRate") = dblFNR,
                            Rcpp::Named("FalsePositiveRate") = dblFPR,
                            Rcpp::Named("TrueNegativeRate") = dblTNR,
                            Rcpp::Named("Accuracy") = dblACC,
                            Rcpp::Named("PositivePredictiveValue") = dblPPV,
                            Rcpp::Named("FalseOmissionRate") = dblFOR,
                            Rcpp::Named("FalseDiscoveryRate") = dblFDR,
                            Rcpp::Named("NegativePredictiveValue") = dblNPV,
                            Rcpp::Named("PositiveLikelihoodRatio") = dblLRp,
                            Rcpp::Named("NegativeLikelihoodRatio") = dblLRn,
                            Rcpp::Named("DiagnosticOddsRatio") = dblDOR );
}



