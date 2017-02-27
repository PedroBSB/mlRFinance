#include <RcppEigen.h>
#include "eiquadprog.h"
#include "KernelMatrix.h"
#include "Utils.h"
#include "progress.hpp"
#include <cmath>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


/***********************************************************************************************/
/*********************************    HEADER FUNCTIONS  ****************************************/
/***********************************************************************************************/

//Define the KernelMatrix function
Eigen::MatrixXd KernelMatrixComputation(Eigen::MatrixXd datMat,
                                        std::string stringValue,
                                        Eigen::RowVectorXd parms);
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
/*********************************     SVC FUNCTIONS    ****************************************/
/***********************************************************************************************/

//TODO: Deixar mais eficiente a normalização zMat=zMat.array().rowwise()/vecSum.transpose().array();
//' @name WOC-SCM
//' @title WOC-SCM - Support Vector Clustering
//' @description Optimize the Lagrange multiplier for the WOC-SCM:
//'
//' Min (1/2)u^{t}Qu+g^{t}u
//' s.t.
//' 0<=u<=wi*C
//' sum ui=1
//' where g=diag(K) and Q=-2K
//' C is the Cost parameter, wi weights for each observation
//'
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param C Cost parameter. Should be C>=0.
//' @param k Total number of clusters.
//' @param sigma Similarity parameter.
//' @param inter Total number of interations.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used, parameters and similarity matrix.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//'             5,5,2,1,
//'             8,1,1,7),nrow=4,ncol=3)
//' svc<-WOCSCM(A, 1, 2, 1, 100, "Gaussian", c(0.5))
//' svc
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite Bicego, Manuele, and Mario AT Figueiredo.
//       "Soft clustering using weighted one-class support vector machines."
//        Pattern Recognition 42.1 (2009): 27-32.
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List WOCSCM(Eigen::MatrixXd X, double C, int k,double sigma,int inter, std::string kernel, Eigen::RowVectorXd parms){
  //Cluster weights
  Eigen::VectorXd gammaWeight(k);
  gammaWeight.fill(1.0/(double)k);
  //Support Vectors
  Eigen::VectorXd SV(X.cols());
  //Z matrix
  Eigen::MatrixXd zMat = Eigen::MatrixXd::Random(X.rows(),k);
  zMat = zMat.cwiseAbs();
  //Get the row sum
  Eigen::VectorXd vecSum = zMat.array().rowwise().sum().eval();
  //Normalize zMat
  for(int l=0;l<zMat.rows();l++){
    zMat.row(l)=zMat.row(l)/vecSum(l);
  }
  Eigen::MatrixXd simVec(X.rows(),k);
  //Initialize logLikelihood
  Eigen::VectorXd llVec(inter);
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(K,1e-10);
  K = nearPDefinite(K, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Training the WOC-SCM
  Eigen::VectorXd g(X.rows());
  g = (-1.0)*K.diagonal();
  //Quadratic programming matrix
  Eigen::MatrixXd Q = (+2.0)*K;
  //RHS equality
  Eigen::VectorXd ce0(1);
  ce0.fill(-1.0);
  //LHS equality
  Eigen::MatrixXd CE = Eigen::MatrixXd::Ones(1,X.rows());
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(X.rows());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(X.rows(),X.rows());
  //Intialize the progressbar
  Progress p(inter, true);
  for(int it=0;it<inter;it++){
    //Verify if everything is ok
    if (Progress::check_abort()) return -1.0;
    //Initialize the loop
    for(int c=0;c<k;c++){
      //RHS: Inequality 2
      Eigen::VectorXd ci2(X.rows());
      ci2.fill(C);
      //Weighted Cost parameter
      ci2=ci2.array()*zMat.col(c).array();
      //Append RHS
      Eigen::VectorXd ci0(2.0*X.rows());
      ci0 << ci1, ci2;
      //Append LHS
      Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
      //Diagonal matrix
      Eigen::VectorXd me(X.rows());
      me.fill(-1.0);
      Eigen::MatrixXd mI = me.asDiagonal();
      //Vertical concatenation
      CI << CI1,
            mI;
      //Get the solution Support Vectors
      SV = rcppeigen_quadratic_solve(Q,g, CE.transpose(),ce0, CI.transpose(), ci0);
      //Get the center of the Hypersphere
      Eigen::RowVectorXd centerA(X.rows());
      centerA = SV.transpose()*X;

      //For each line
      Eigen::VectorXd dist0 = K.diagonal();
      Eigen::VectorXd dist1 = (-2.0)*(SV.transpose()*K).array();
      double dist2 = SV.transpose()*K*SV;
      Eigen::VectorXd zVec = dist0.array()+dist1.array()+dist2;

      zVec=(-1.0)*zVec/sigma;
      //For each cluster
      zVec = zVec.array().exp()*gammaWeight(c);
      //Store the similarity
      simVec.col(c) = zVec;
      //Store the column
      zMat.col(c) = zVec;
    }

    //M-Step
    //Get the row sum
    vecSum = zMat.array().rowwise().sum().eval();
    //Normalize zMat
    for(int l=0;l<zMat.rows();l++){
      zMat.row(l)=zMat.row(l)/vecSum(l);
    }

    //Update gammaWeight
    gammaWeight = (1.0/zMat.rows())*zMat.colwise().sum();

    //Compute the log-likelihood
    double ll=simVec.colwise().sum().array().log().sum();
    llVec(it)=ll;

    //Increment the progress bar
    p.increment();
  }

  //Return the results
  return Rcpp::List::create(Rcpp::Named("LogLikelihood") = llVec,
                            Rcpp::Named("Zmat") = zMat,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


double R2Distance(Eigen::MatrixXd K, Eigen::MatrixXd X, Eigen::VectorXd SV, Eigen::RowVectorXd x, std::string kernel, Eigen::RowVectorXd parms){
  //Calculate the R2 distance.
  double res=0;
  double K1 = KernelMatrixComputationValue(x, x, kernel, parms);
  Eigen::VectorXd Krow = KernelMatrixComputationPred(X,x,kernel,parms);
  double K2 = SV.transpose()*Krow;
  double K3 = SV.transpose()*K*SV;
  res = K1-2*K2+K3;
  return(res);
}


// [[Rcpp::export]]
Rcpp::List CSVC(Eigen::MatrixXd X, double C, std::string kernel, Eigen::RowVectorXd parms){
  //Support Vectors
  Eigen::VectorXd SV(X.cols());
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);
  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(K,1e-10);
  K = nearPDefinite(K, 1e+6, 1e-06, 1e-07, 1e-08, true);
  //Training the WOC-SCM
  Eigen::VectorXd g(X.rows());
  g = (-1.0)*K.diagonal();
  //Quadratic programming matrix
  Eigen::MatrixXd Q = (+2.0)*K;
  //RHS equality
  Eigen::VectorXd ce0(1);
  ce0.fill(-1.0);
  //LHS equality
  Eigen::MatrixXd CE = Eigen::MatrixXd::Ones(1,X.rows());
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(X.rows());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(X.rows(),X.rows());
  //RHS: Inequality 2
  Eigen::VectorXd ci2(X.rows());
  ci2.fill(C);
  //Append RHS
  Eigen::VectorXd ci0(2.0*X.rows());
  ci0 << ci1, ci2;
  //Append LHS
  Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
  //Diagonal matrix
  Eigen::VectorXd me(X.rows());
  me.fill(-1.0);
  Eigen::MatrixXd mI = me.asDiagonal();
  //Vertical concatenation
  CI << CI1,
        mI;
  //Get the solution Support Vectors
  SV = rcppeigen_quadratic_solve(Q,g, CE.transpose(),ce0, CI.transpose(), ci0);

  //Get the center of the Hypersphere
  Eigen::RowVectorXd centerA(X.rows());
  centerA = SV.transpose()*X;
  //Find the Support Vectors
  double R2 = 0.0;
  double cont = 0 ;
  for(int i=0;i<X.rows();i++){
    //Support Vector 0 < lambda < C
    if(SV(i)>1e-5 & SV(i)< C-1e-5){
      R2 = R2 + R2Distance(K, X, SV, X.row(i), kernel, parms);
      cont = cont+1.0;
    }
  }
  R2 = R2/cont;

  //Caculating adjacency matrix
  R2 = R2+1e-7;
  int inter = (X.rows()*(X.rows()-1)/2);
  //Intialize the progressbar
  Progress p(inter, true);
  //Create the Adjancet Matrix
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(X.rows(),X.rows());
  for(int i=0;i<X.rows();i++){
    Eigen::RowVectorXd x = X.row(i);
    for(int j=(i+1);j<X.rows();j++){
      Eigen::RowVectorXd y = X.row(i);
      //Number of segments
      int iSegments=20;
      double k = 1.0/(double) iSegments;
      //Boolean connected
      int conect = 0;
      //Check every point between them
      for(int s=0;s<iSegments;s++){
        //Calculate the distance
        Eigen::RowVectorXd z = x+k*(y-x);
        double dist = R2Distance(K, X, SV, z, kernel, parms);
        //Increment k
        k = k+(1.0/(double)iSegments);
        if(dist <= R2){
          conect=conect+1;
        }
      }
      if(conect==iSegments){
        A(i,j)=1;
        A(j,i)=1;
      }
      //Increment the progress bar
      p.increment();
    }
  }


  //Return the results
  return Rcpp::List::create(Rcpp::Named("AdjacencyMatrix") = A,
                            Rcpp::Named("SupportVectors") = SV,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}


//TODO: Deixar mais eficiente a normalização zMat=zMat.array().rowwise()/vecSum.transpose().array();
//' @name WOC-SCM
//' @title WOC-SCM - Support Vector Clustering
//' @description Optimize the Lagrange multiplier for the WOC-SCM:
//'
//' Min (1/2)u^{t}Qu+g^{t}u
//' s.t.
//' 0<=u<=wi*C
//' sum ui=1
//' where g=diag(K) and Q=-2K
//' C is the Cost parameter, wi weights for each observation
//'
//' @param X Numeric matrix with the explanatory variables. Dimension equal NxP
//' @param wMat Weight Numeric matrix. Dimension equal NxN
//' @param C Cost parameter. Should be C>=0.
//' @param k Total number of clusters.
//' @param sigma Similarity parameter (between 0 and 1).
//' @param inter Total number of interations.
//' @param parms Parameters associated with chosen kenel.
//' @return List Support Vectors, Kernel used, parameters and similarity matrix.
//' If the results for the Support Vectors are NaN it means that
//' there is no Support Vector and the Quadratic Programming Problem
//' is unfeasible.
//' @examples
//'
//' A<-matrix(c(1,2,5,6,
//'             5,5,2,1,
//'             8,1,1,7),nrow=4,ncol=3)
//' svc<-WOCSCM(A, 1, 2, 1, 100, "Gaussian", c(0.5))
//' svc
//' @seealso See \code{\link{.CallOctave}}, \code{\link{o_source}}, \code{\link{o_help}}
// @cite Bicego, Manuele, and Mario AT Figueiredo.
//       "Soft clustering using weighted one-class support vector machines."
//        Pattern Recognition 42.1 (2009): 27-32.
// @bibliography ~/vignettes/bibliography.bib
// [[Rcpp::export]]
Rcpp::List SpatialWOCSCM(Eigen::MatrixXd X,Eigen::MatrixXd wMat, double C, int k,double gamma1, double gamma2,int inter, std::string kernel, Eigen::RowVectorXd parms){
  //Cluster weights
  Eigen::VectorXd gammaWeight(k);
  gammaWeight.fill(1.0/(double)k);
  //Support Vectors
  Eigen::VectorXd SV(X.cols());
  //Z matrix
  Eigen::MatrixXd zMat = Eigen::MatrixXd::Random(X.rows(),k);
  zMat = zMat.cwiseAbs();
  //Get the row sum
  Eigen::VectorXd vecSum = zMat.array().rowwise().sum().eval();
  //Normalize zMat
  for(int l=0;l<zMat.rows();l++){
    zMat.row(l)=zMat.row(l)/vecSum(l);
  }
  Eigen::MatrixXd simVec(X.rows(),k);
  //Initialize logLikelihood
  Eigen::VectorXd llVec(inter);
  //Create the Kernel Matrix
  Eigen::MatrixXd K = KernelMatrixComputation(X,kernel,parms);

  //Nearest positive semidefinite matrix in terms of Frobenius norm
  //nearPositiveDefinite(K,1e-10);
  K = nearPDefinite(K, 1e+6, 1e-06, 1e-07, 1e-08, true);

    //Training the WOC-SCM
  Eigen::VectorXd g(X.rows());
  g = (-1.0)*K.diagonal();
  //Quadratic programming matrix
  Eigen::MatrixXd Q = (+1.0)*K;
  //RHS equality
  Eigen::VectorXd ce0(1);
  ce0.fill(-1.0);
  //LHS equality
  Eigen::MatrixXd CE = Eigen::MatrixXd::Ones(1,X.rows());
  //RHS: Inequality 1
  Eigen::VectorXd ci1 = Eigen::VectorXd::Zero(X.rows());
  //LHS: Inequality 1
  Eigen::MatrixXd CI1 = Eigen::MatrixXd::Identity(X.rows(),X.rows());
  //Intialize the progressbar
  Progress p(inter, true);
  for(int it=0;it<inter;it++){
    //Verify if everything is ok
    if (Progress::check_abort()) return -1.0;
    //Initialize the loop
    for(int c=0;c<k;c++){
      //RHS: Inequality 2
      Eigen::VectorXd ci2(X.rows());
      ci2.fill(C);
      //Weighted Cost parameter
      ci2=ci2.array()*zMat.col(c).array();
      //Append RHS
      Eigen::VectorXd ci0(2.0*X.rows());
      ci0 << ci1, ci2;
      //Append LHS
      Eigen::MatrixXd CI(CI1.rows()+CI1.rows(), CI1.cols());
      //Diagonal matrix
      Eigen::VectorXd me(X.rows());
      me.fill(-1.0);
      Eigen::MatrixXd mI = me.asDiagonal();
      //Vertical concatenation
      CI << CI1,
            mI;
      //Get the solution Support Vectors
      SV = rcppeigen_quadratic_solve(Q,g, CE.transpose(),ce0, CI.transpose(), ci0);

      //Get the center of the Hypersphere
      Eigen::RowVectorXd centerA(X.rows());
      centerA = SV.transpose()*X;

      //For each line
      Eigen::VectorXd dist0 = K.diagonal();
      Eigen::VectorXd dist1 = (-2.0)*(SV.transpose()*K).array();
      double dist2 = SV.transpose()*K*SV;
      Eigen::VectorXd zVec1 = dist0.array()+dist1.array()+dist2;

      //Weight based on distance
      Eigen::VectorXd zVec2 = zVec1.transpose()*wMat;

      //Smooth factor
      zVec1=(-1.0)*zVec1*(1/gamma1);

      //Smooth factor
      zVec2=(-1.0)*zVec2*(1/gamma2);

      //Join the two similarities
      Eigen::VectorXd zVec = zVec1+zVec2;

      //For each cluster
      zVec = zVec.array().exp()*gammaWeight(c);

      //Store the similarity
      simVec.col(c) = zVec;

      //Store the column
      zMat.col(c) = zVec;
    }

    //M-Step
    //Get the row sum
    vecSum = zMat.array().rowwise().sum().eval();
    //Normalize zMat
    for(int l=0;l<zMat.rows();l++){
      zMat.row(l)=zMat.row(l)/vecSum(l);
    }

    //Update gammaWeight
    gammaWeight = (1.0/zMat.rows())*zMat.colwise().sum();

    //Compute the log-likelihood
    double ll=simVec.colwise().sum().array().log().sum();
    llVec(it)=ll;

    //Increment the progress bar
    p.increment();
  }

  //Return the results
  return Rcpp::List::create(Rcpp::Named("LogLikelihood") = llVec,
                            Rcpp::Named("Zmat") = zMat,
                            Rcpp::Named("Kernel") = kernel,
                            Rcpp::Named("Parameters") = parms);
}




