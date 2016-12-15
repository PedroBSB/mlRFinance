#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
using namespace Rcpp;

/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/

//Define the KernelMatrix function
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,
                                        std::string stringValue,
                                        arma::vec parms);
//Define the Solver for Quadratic Programming
Eigen::VectorXd rcppeigen_quadratic_solve(Eigen::MatrixXd & G,
                                          Eigen::VectorXd & g0,
                                          const Eigen::MatrixXd & CE,
                                          const Eigen::VectorXd & ce0,
                                          const Eigen::MatrixXd & CI,
                                          const Eigen::VectorXd & ci0);
//Test if the matrix is Positive Definite
bool IsPositiveDefinite(Eigen::MatrixXd mat);
//nearest positive semidefinite matrix in terms of Frobenius norm
void nearPositiveDefinite(Eigen::MatrixXd &mat,double noise);
//Nearest positive semidefinite matrix (Matrix::nearPD)
Eigen::MatrixXd nearPDefinite(Eigen::MatrixXd mat, int maxit, double eigtol, double conv_tol, double posd_tol, bool keepDiagonal);
//Add some noise to the matrix diagonal
void addNoise(Eigen::MatrixXd &mat,double noise);
//Print Object at Console
void PrintObject(Eigen::MatrixXd mat);
void PrintObject(Eigen::VectorXd vec);


/***********************************************************************************************/
/*********************************     SVM FUNCTIONS    ****************************************/
/***********************************************************************************************/

//Platt, J.: 2000, ‘Probabilistic outputs for support vector machines and comparison
//to regularized likelihood methods’. In: A. Smola, P. Bartlett, B. Sch¨olkopf, and
//D. Schuurmans (eds.): Advances in Large Margin Classifiers. Cambridge, MA.
// http://www.csie.ntu.edu.tw/~htlin/paper/doc/plattprob.pdf
Eigen::VectorXd predictProbability(Eigen::VectorXd predVec,Eigen::VectorXd y, Eigen::VectorXd SV){
  Eigen::VectorXd parms(2);
  //Initial parameters
  double A=0;
  double B=0;
  //Parameter setting
  int maxiter=100; //Maximum number of iterations
  double minstep=1e-10; //Minimum step taken in line search
  double sigma=1e-12; //Set to any value > 0
  double prior0=0;
  double prior1=0;
  //Count the number of elements
  for(int i=0;i<y.size();i++){
    if(y(i)<0){
      prior0=prior0+1;
    }
    else{
      prior1=prior1+1;
    }
  }
  //Construct initial values: target support in array t,
  // initial function value in fval
  double hiTarget=(prior1+1.0)/(prior1+2.0);
  double loTarget=1.0/(prior0+2.0);
  int len=prior1+prior0; // Total number of data
  Eigen::VectorXd t(len);
  for(int i = 0;i<len;i++) {
    if (y(i) > 0)
    {
      t(i)=hiTarget;
    }
    else{
      t(i)=loTarget;
    }
  }

  A=0.0;
  B=log((prior0+1.0)/(prior1+1.0));
  double fval=0.0;
  double fApB=0.0;
  for(int i = 0;i<len;i++){
    fApB=predVec(i)*A+B;
    if (fApB >= 0.0){
      fval = fval+t(i)*fApB+std::log(1+std::exp(-fApB));
    }
    else{
      fval = fval+(t(i)-1)*fApB+std::log(1+std::exp(fApB));
    }
  }
  int it=0;
  for(it = 0;it<maxiter;it++){
    //Update Gradient and Hessian (use H’ = H + sigma I)
    double h11=sigma;
    double h22=sigma;
    double h21=0.0;
    double g1=0.0;
    double g2=0.0;
    double p=0.0;
    double q=0.0;
    for(int i = 0;i<len;i++){
        fApB=predVec(i)*A+B;
        if (fApB >= 0){
          p=std::exp(-fApB)/(1.0+std::exp(-fApB));
          q=1.0/(1.0+std::exp(-fApB));
        }
        else{
          p=1.0/(1.0+std::exp(fApB));
          q=std::exp(fApB)/(1.0+std::exp(fApB));
        }
          double d2=p*q;
          h11 = h11+ predVec(i)*predVec(i)*d2;
          h22 += h22+d2;
          h21 += predVec(i)*d2;
          double d1=t[i]-p;
          g1 =g1+predVec(i)*d1;
          g2 = g2+d1;
      }
      if (std::abs(g1)<1e-5 && std::abs(g2)<1e-5){//Stopping criteria
        break;
      }
      //Compute modified Newton directions
      double det=h11*h22-h21*h21;
      double dA=-(h22*g1-h21*g2)/det;
      double dB=-(-h21*g1+h11*g2)/det;
      double gd=g1*dA+g2*dB;
      double stepsize=1.0;
      double newA=0;
      double newB=0;
      double newf=0.0;
      while (stepsize >= minstep){ //Line search
        newA=A+stepsize*dA;
        newB=B+stepsize*dB;
        newf=0.0;
        for(int i = 0;i<len;i++){
          fApB=predVec(i)*newA+newB;
          if (fApB >= 0){
            newf = newf+t(i)*fApB+std::log(1+std::exp(-fApB));
          }
          else{
            newf = newf+(t(i)-1)*fApB+std::log(1+std::exp(fApB));
          }
        }
        if (newf<fval+0.0001*stepsize*gd){
          A=newA;
          B=newB;
          fval=newf;
          break; //Sufficient decrease satisfied
        }
        else{
          stepsize = stepsize/2.0;
        }
        if (stepsize < minstep){
          //Didn't work.
          //print ’Line search fails’
          break;
        }
      }
    }
    if (it >= maxiter){
      //print ’Reaching maximum iterations’
    }
    //Results
    parms(0)=A;
    parms(1)=B;
    return(parms);
}


/************************************ C-SVM L1 *************************************************/
//' @name CSVML1
//' @title C-SVM L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVM L1:
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' 0<=u<=C1
//'
//' where d is the vector of dependent variable,
//' and Q=K.*(d*t(d))=DKD. C is the Cost parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML1(d, A, 1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List CSVML1(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, arma::vec parms, bool biasTerm){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0;
  //LHS equality
  Eigen::MatrixXd CE;
  if(biasTerm==true){
    //RHS equality
    ce0 = Eigen::VectorXd::Zero(1);
    //LHS equality
    CE = Eigen::MatrixXd::Ones(1,y.size());
  }
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(y.size(),y.size());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(y.size());
  ci2.fill(C);
  //Append RHS
  Eigen::VectorXd ci0(2.0*y.size());
  ci0 << ci1, ci2;
  //Append LHS
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
  //Diagonal matrix
  Eigen::VectorXd me(y.size());
  me.fill(-1.0);
  Eigen::MatrixXd mI = me.asDiagonal();
  //Vertical concatenation
  CI << CI1,
        mI;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(Q,1e-10);
  Q = nearPDefinite(Q, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE.transpose(),ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//' @name Predicted CSVML1
//' @title C-SVM L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVR L1:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVML1 List of Results of the CSVML1
//' @param y Numeric vector with the response variable. Dimension equal Nx1
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param Xprev Numeric matrix with the explanatory variables (predicted). Dimension equal MxP
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @param typePredict 0-Binary(-1,+1), 1-Probability, 2- Raw result
//' @return Eigen::VectorXd with the predicted values for Xpred
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
Eigen::VectorXd PredictedCSVML1(Rcpp::List CSVML1,Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Xprev, int typePredict, bool biasTerm){

  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVML1["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVML1["Kernel"]);

  //Get the parameters
  arma::vec parms = as<arma::vec> (CSVML1["Parameters"]);

  //Bias Term
  double gamma=0.0;
  if(biasTerm==true){
    for(int i=0;i<X.rows();i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      gamma=gamma+(y(i)-res);
    }
  }
  gamma = gamma / X.rows();

  //Total number of observations
  int size = Xprev.rows();
  Eigen::VectorXd predVec(size);
  if(typePredict==1){
    for(int i=0;i<size;i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      predVec(i)=res-gamma;
    }
      //Return the probability
    Eigen::VectorXd parms = predictProbability(predVec, y, SV);
    double A=parms(0);
    double B=parms(1);
    //Normalize the predVec;
    for(int i=0;i<size;i++){
      predVec(i)=1.0/(1+std::exp(A*predVec(i)+B));
    }
  }
  else{
    for(int i=0;i<size;i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      res = res - gamma;
      if(typePredict==0){
        //Return the signal
        if(res<0){
          predVec(i)=-1.0;
        }
        else{
          predVec(i)=+1.0;
        }
      }
      else{
        //Return the raw forecast
        res = res - gamma;
        predVec(i)=res;
      }
    }
  }
  return(predVec);
}


//' @name Pseudo R2 - Predicted CSVRL1
//' @title C-SVR L1 - Support Vector Regression with C cost and L1 regularization.
//' @description Prediction for the C-SVR L1:
//'
//' f(x)=Sum_{i=1}^{N}(lambda*-lambda)K(x_{i},x)
//' @param CSVRL1 List of Results of the CSVRL1
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return Eigen::VectorXd with the Pseudo R2 for each variable.
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
Eigen::VectorXd R2PredictedCSVML1(Rcpp::List CSVML1,Eigen::VectorXd y, Eigen::MatrixXd X, int typePredict, bool biasTerm){
  //Results
  Eigen::VectorXd R2vec(X.cols());

  //Get the SV
  Eigen::VectorXd SV = as<Eigen::VectorXd> (CSVML1["SupportVectors"]);

  //Get the kernel
  std::string kernel = as<std::string> (CSVML1["Kernel"]);

  //Get the parameters
  arma::vec parms = as<arma::vec> (CSVML1["Parameters"]);

  //Bias Term
  double gamma=0.0;
  if(biasTerm==true){
    for(int i=0;i<X.rows();i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      gamma=gamma+(y(i)-res);
    }
  }
  gamma = gamma / X.rows();

  //Prediction for the full model
  //Total number of observations
  int size = X.rows();
  Eigen::VectorXd predVec(size);
  if(typePredict==1){
    for(int i=0;i<size;i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      predVec(i)=res - gamma;
    }
    //Return the probability
    Eigen::VectorXd parms = predictProbability(predVec, y, SV);
    double A=parms(0);
    double B=parms(1);
    //Normalize the predVec;
    for(int i=0;i<size;i++){
      predVec(i)=1.0/(1+std::exp(A*predVec(i)+B));
    }
  }
  else{
    for(int i=0;i<size;i++){
      //Create the Kernel Matrix
      Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
      Eigen::VectorXd F = y.array()*SV.array() *K.array();
      double res = F.sum();
      res = res - gamma;
      if(typePredict==0){
        //Return the signal
        if(res<0){
          predVec(i)=-1.0;
        }
        else{
          predVec(i)=+1.0;
        }
      }
      else{
        //Return the raw forecast
        predVec(i)=res-gamma;
      }
    }
  }

  //Sum of squared errors
  double SSE = predVec.squaredNorm();

  //For each variable:
  for(int v=0;v<X.cols();v++){
    //Zero columns
    Eigen::MatrixXd Xprev = X;

    //Zero the variable
    Xprev.col(v).fill(0.0);

    //Bias Term
    gamma=0.0;
    if(biasTerm==true){
      for(int i=0;i<X.rows();i++){
        //Create the Kernel Matrix
        Eigen::VectorXd K = KernelMatrixComputationPred(X,X.row(i),kernel,parms);
        Eigen::VectorXd F = y.array()*SV.array() *K.array();
        double res = F.sum();
        gamma=gamma+(y(i)-res);
      }
    }
    gamma = gamma / X.rows();

    //Total number of observations
    int size = X.rows();
    Eigen::VectorXd predVec2(size);
    if(typePredict==1){
      for(int i=0;i<size;i++){
        //Create the Kernel Matrix
        Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
        Eigen::VectorXd F = y.array()*SV.array() *K.array();
        double res = F.sum();
        predVec2(i)=res-gamma;
      }
      //Return the probability
      Eigen::VectorXd parms = predictProbability(predVec2, y, SV);
      double A=parms(0);
      double B=parms(1);
      //Normalize the predVec;
      for(int i=0;i<size;i++){
        predVec2(i)=1.0/(1+std::exp(A*predVec2(i)+B));
      }
    }
    else{
      for(int i=0;i<size;i++){
        //Create the Kernel Matrix
        Eigen::VectorXd K = KernelMatrixComputationPred(X,Xprev.row(i),kernel,parms);
        Eigen::VectorXd F = y.array()*SV.array() *K.array();
        double res = F.sum();
        res = res - gamma;
        if(typePredict==0){
          //Return the signal
          if(res<0){
            predVec2(i)=-1.0;
          }
          else{
            predVec2(i)=+1.0;
          }
        }
        else{
          //Return the raw forecast
          predVec2(i)=res-gamma;
        }
      }
    }

    double SSEvar = predVec2.squaredNorm();
    R2vec(v) = SSEvar/SSE;
  }
  return(R2vec);
}


/************************************ C-SVM L2 *************************************************/

//' @name CSVML2
//' @title C-SVM L2 - Support Vector Regression with C cost and L2 regularization.
//' @description Optimize the Lagrange multiplier for the C-SVM L2:
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' u>=0
//'
//' where d is the vector of dependent variable,
//' and Q=(K+I/C).*(t(d)*d). C is the Cost parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- CSVML2(d, A, 1, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List CSVML2(Eigen::VectorXd y, Eigen::MatrixXd X, double C, std::string kernel, arma::vec parms, bool biasTerm){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0;
  //LHS equality
  Eigen::MatrixXd CE;
  if(biasTerm==true){
    //RHS equality
    ce0 = Eigen::VectorXd::Zero(1);
    //LHS equality
    CE = Eigen::MatrixXd::Zero(1,y.size());
    CE.row(0) = y;
  }
  //RHS: Inequality 1
  Eigen::VectorXd ci0 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI = Eigen::MatrixXd::Identity(y.size(),y.size());
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  K=K+(CI/C);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(Q,1e-10);
  Q = nearPDefinite(Q, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE.transpose(),ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//' @name nuSVM
//' @title nu-SVM - Support Vector Regression with nu parameter.
//' @description The m-support vector classification (Scholkopf, Smola,
//' Williamson, & Bartlett, 2000) uses a new parameter nu which controls
//' the number of support vectors and training errors. The
//' parameter nu in (0, 1] is an upper bound on the fraction of training
//' errors and a lower bound of the fraction of support vectors.
//'
//' Min (1/2)u^{t}Qu-1^{t}u
//' s.t.
//' d^{t}*u=0
//' nu <=1^t*u
//' 0<=u<=1/l
//'
//' where d is the vector of dependent variable,
//' and Q=K.*(t(d)*d). nu is the parameter.
//'
//' @param y Vector with dependent variables should be -1 or +1. Dimension equal Nx1.
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>0.
//' @param kernel Name of the kernel that will be used.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used and parameters.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//' 2,4,1,2),nrow=4,ncol=2)
//' d<-c(-1,-1,+1,-1)
//' svm1<- nuSVM(d, A, 0.2, "Gaussian", c(0.5))
//'
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite soman2009machine @chang2001training
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List nuSVM(Eigen::VectorXd y, Eigen::MatrixXd X, double nu, std::string kernel, arma::vec parms){
  //Support Vectors
  Eigen::VectorXd SV(y.size());
  //Create the one vector Nx1
  Eigen::VectorXd e = Eigen::VectorXd::Ones(y.size());
  e=(-1.0)*e;
  //RHS equality
  Eigen::VectorXd ce0 = Eigen::VectorXd::Zero(1);
  //LHS equality
  Eigen::MatrixXd CE(1,y.size());
  CE = y;
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(y.size());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(y.size(),y.size());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(y.size());
  ci2.fill(1.0/y.size());
  //Append RHS
  Eigen::VectorXd ci0(2.0*y.size()+1);
  ci0 << ci1, ci2, -nu;
  //Append LHS
  Eigen::MatrixXd CI3(1,y.size());
  CI3.fill(1);
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows()+1, CI1.cols());
  CI <<        CI1,
        (-1.0)*CI1,
               CI3;
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Create matrix D
  Eigen::MatrixXd D = y*y.transpose();
  //Create matrix Q
  Eigen::MatrixXd Q = K.cwiseProduct(D);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //  nearPositiveDefinite(Q,1e-10);
  Q = nearPDefinite(Q, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,e, CE,ce0, CI.transpose(), ci0);
  //Return the results
  return Rcpp::List::create(Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


// [[Rcpp::export]]
Eigen::VectorXd solveTest(Eigen::MatrixXd Dmat, Eigen::VectorXd dvec, Eigen::MatrixXd Amat,Eigen::VectorXd bvec, Eigen::MatrixXd CE, Eigen::VectorXd ce){
  //Convert to the quadprog
  Eigen::VectorXd d0=-dvec;
  Eigen::VectorXd b0=-bvec;
  Eigen::VectorXd ce0=-ce;
  //Get the solution
  Eigen::VectorXd x = rcppeigen_quadratic_solve(Dmat,d0, CE,ce0, Amat, b0);
  return(x);
}


// [[Rcpp::export]]
Eigen::VectorXd solveTest2(Eigen::MatrixXd Dmat, Eigen::VectorXd dvec, Eigen::MatrixXd Amat,Eigen::VectorXd bvec){
  Eigen::MatrixXd CE;
  Eigen::VectorXd ce;
  //Convert to the quadprog
  Eigen::VectorXd d0=-dvec;
  Eigen::VectorXd b0=-bvec;
  Eigen::VectorXd ce0=-ce;
  //Get the solution
  Eigen::VectorXd x = rcppeigen_quadratic_solve(Dmat,d0, CE,ce0, Amat, b0);
  return(x);
}
