#include <RcppEigen.h>
#include "Utils.h"
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>

/********************************************************************************************************************/
/**********************************************        SVR      *****************************************************/
/********************************************************************************************************************/

//Mean Square Error
double MSEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  double res = (y-yPred).squaredNorm();
  return(res);
}

//Mean Absolute Percentage Error
double MAPEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  Eigen::VectorXd res = (y-yPred).cwiseAbs();
  res = res.array()/y.array();
  return(res.mean());
}

//Theil-U statistic.
// [[Rcpp::export]]
double TheilUfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  int n=y.size()-1;
  //Numerator
  double num = ((y.tail(n).array()-yPred.tail(n).array()).square()/y.tail(n).array()).sum();
  //Denominator
  double den = ((y.head(n).array()-y.tail(n).array()).square()/y.head(n).array()).sum();
  double res = std::sqrt(num/den);
  return(res);
}

//Mean Error
double MEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  double res = (y-yPred).mean();
  return(res);
}

//Mean Absolute Error (MAE)
double MAEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  double res = (y-yPred).cwiseAbs().mean();
  return(res);
}

//Mean Percentage Error
double MPEfunction(Eigen::VectorXd y, Eigen::VectorXd yPred){
  Eigen::VectorXd res = (y-yPred);
  res = res.array()/y.array();
  return(res.mean());
}

//Harvey’s R2 statistic that uses the random walk model for comparison
double RWR2(Eigen::VectorXd y, Eigen::VectorXd yPred){
  int n=y.size()-1;
  double mu = y.mean();
  double rwVec = (y.tail(n).array()-y.head(n).array()-mu).square().sum();
  double vec = (y.array()-yPred.array()).square().sum();
  double res = 1.0 - (n/(n+1))*(vec/rwVec);
  return(res);
}

// [[Rcpp::export]]
Rcpp::List ErrorMeasures(Eigen::VectorXd y, Eigen::VectorXd yPred){
  //Calculate Mean Square Error
  double mse = MSEfunction(y, yPred);
  //Calculate Mean Absolute Percentage Error
  double mape = MAPEfunction(y, yPred);
  //Calculate Theil-U statistic
  double theilU = TheilUfunction(y, yPred);
  //Calculate Theil-U statistic
  double me = MEfunction(y, yPred);
  //Calculate Theil-U statistic
  double mae = MAEfunction(y, yPred);
  //Calculate Theil-U statistic
  double mpe = MPEfunction(y, yPred);
  //RWR2
  double rwr2 = RWR2(y, yPred);

  //Return the results
  return Rcpp::List::create(Rcpp::Named("MSE") = mse,
                            Rcpp::Named("MAPE") = mape,
                            Rcpp::Named("TheilU") = theilU,
                            Rcpp::Named("ME") = me,
                            Rcpp::Named("MAE") = mae,
                            Rcpp::Named("MPE") = mpe,
                            Rcpp::Named("RWR2") = rwr2
  );
}

/********************************************************************************************************************/
/**********************************************        SVWQR    *****************************************************/
/********************************************************************************************************************/



/********************************************************************************************************************/
/**********************************************   SVM   *************************************************************/
/********************************************************************************************************************/

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



