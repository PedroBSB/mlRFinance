#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
#include "conversion.h"
using namespace Rcpp;
/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/
void PrintObjectLine(arma::vec mat){
  int nelem=mat.n_elem;
  for(int i=0;i<nelem;i++){
    std::cout<<mat(i)<<" ";
  }
  std::cout<<std::endl;
}

void PrintObjectLine(arma::uvec mat){
  for(int i=0;i<mat.n_elem;i++){
    std::cout<<mat(i)<<" ";
  }
  std::cout<<std::endl;
}


void PrintObject(Eigen::VectorXd mat){
  for(int i=0;i<mat.size();i++){
    std::cout<<mat(i)<<" ";
  }
  std::cout<<std::endl;
}

//Print Object at Console
void PrintObject(Eigen::MatrixXd mat){
  for(int i=0;i<mat.cols();i++){
    for(int j=0;j<mat.rows();j++){
      std::cout<<mat(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }
}

//Print Object at Console
void PrintObject(arma::mat mat){
  for(int i=0;i<mat.n_rows;i++){
    for(int j=0;j<mat.n_cols;j++){
      std::cout<<mat(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }
}


//Print time
void PrintTime(){
  // current date/time based on current system
  time_t now = time(0);
  // convert now to string form
  char* dt = ctime(&now);
  std::cout << "The local date and time is: " << dt << std::endl;
}

//Identify if the matrix is Positive Definite
bool IsPositiveDefinite(Eigen::MatrixXd mat){
  //Set the result
  bool res=false;
  //Setting the tolerance
  double tol=1e-08;
  //Getting the eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXd> es(mat,true);
  for(int i=0;i<mat.rows();i++){
    //Casting to complex
    std::complex<double> lambda = es.eigenvalues()[i];
    if (std::abs(lambda.real()) < tol) {
      lambda.real(0);
    }
    else{
      if(lambda.real()<0){
        return(true);
      }
    }
  }
return(res);
}


//Adding some noise to diagonal
void addNoise(Eigen::MatrixXd &mat,double noise){
  //Create some perturbation
  Eigen::VectorXd me1(mat.cols());
  me1.fill(noise);
  //Adding the noise
  Eigen::MatrixXd mI1 = me1.asDiagonal();
  mat=mat+mI1;
}


//Construct the the nearest positive semidefinite matrix in terms of Frobenius norm
//Higham, Nicholas J. "Computing a nearest symmetric positive semidefinite matrix."
//Linear algebra and its applications 103 (1988): 103-118.
void nearPositiveDefinite(Eigen::MatrixXd &mat,double noise){
  //Set the result
  bool res=false;
  //Setting the tolerance
  double tol=1e-08;
  //Getting the eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXd> es(mat,true);
  Eigen::MatrixXcd Lambda = es.eigenvalues().asDiagonal();
  //Real part
  Eigen::MatrixXd rLambda = Lambda.real();
  for(int i=0;i<mat.cols();i++){
   if(rLambda(i,i)<tol){
     rLambda(i,i)=noise;
    }
  }
  //Spectral Decomposition
  mat = es.eigenvectors().real()*rLambda*es.eigenvectors().real().inverse();
}

//Remove columns
Eigen::MatrixXd removeColumns(Eigen::MatrixXd Q, Eigen::VectorXd d, double e1){
  //Find the number of columns
  int ncols = (d.array()>e1).count();
  //Create the new matrix
  Eigen::MatrixXd Q0(Q.rows(),ncols);
  //Initialize the counter
  int cont=0;
  for(int c=0;c<Q.cols();c++){
    if(d(c)>e1){
      for(int r=0;r<Q.rows();r++){
        Q0(r,cont) = Q(r,c);
      }
      cont = cont+1;
    }
  }
  return(Q0);
}

//Repeated Vector
Eigen::MatrixXd repetVector(Eigen::VectorXd d,double e1, int nrows){
  //Find the number of columns
  int ncols = (d.array()>e1).count();
  //Initialize the new matrix
  Eigen::MatrixXd vecMat(nrows,ncols);
  //Intialize the vector column
  Eigen::VectorXd vec(nrows);
  //Initialize the counter
  int cont=0;
  for(int i=0;i<d.size();i++){
    if(d(i)>e1){
      double dbl = d(i);
      vec.fill(dbl);
      vecMat.col(cont) = vec;
      cont=cont+1;
    }
  }
  return(vecMat);
}


//Near Positive Definite matrix
// [[Rcpp::export]]
Eigen::MatrixXd nearPDefiniteOld(Eigen::MatrixXd mat, int maxit=1e+6, double eigtol = 1e-06, double conv_tol = 1e-07, double posd_tol = 1e-08, bool keepDiagonal=false){
  //Number of columns
  int n = mat.cols();
  //Store the diagonal
  Eigen::ArrayXd diagX0 = mat.diagonal();
  //Dykstra matrix
  Eigen::MatrixXd D_S = Eigen::MatrixXd::Zero(mat.rows(),mat.cols());
  //Store the original matrix
  Eigen::MatrixXd X = mat;
  int iter = 0;
  bool converged = false;
  double conv = std::numeric_limits<double>::infinity();
  while (iter < maxit && !converged) {
    Eigen::MatrixXd Y = X;
    Eigen::MatrixXd R = Y - D_S;
    //Getting the eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> e(R,true);
    Eigen::MatrixXd Q = e.eigenvectors().real();
    Eigen::VectorXd d = e.eigenvalues().real();
    double e1 = eigtol*d(0);
    //Test if the matrix seems be negative definite
    if (!(d.array()>e1).any()) stop("Matrix seems negative semi-definite");
    //Remove columns with d <= eig.tol * d[0]
    Eigen::MatrixXd Q0 = removeColumns(Q,d,e1);
    //Create repeated vector
    Eigen::MatrixXd repVec = repetVector(d,e1,Q0.rows());
    //Elementwise multiplication
    Eigen::MatrixXd Q0temp = Q0.cwiseProduct(repVec);
    //Calculate the tcrossprod
    X = Q0temp*Q0.transpose();
    //Dykstra's correction
    D_S = X - R;
    //Keep the same diagonal
    if(keepDiagonal) X.diagonal() = mat.diagonal();
    //Get the infinity norm
    double convNum = (Y-X).cwiseAbs().rowwise().sum().maxCoeff();
    double convDem = Y.cwiseAbs().rowwise().sum().maxCoeff();
    conv = convNum/convDem;
    //Update the interaction and convergence;
    iter = iter + 1;
    //Update the convergence criteria
    converged = (conv <= conv_tol);
  }

  //Eigen step applied to the result of the Higham (2002) algorithm.
  //Get the eigenvalues
  Eigen::EigenSolver<Eigen::MatrixXd> es(X,true);
  Eigen::MatrixXcd Lambda = es.eigenvalues().asDiagonal();
  //Real part
  Eigen::MatrixXd rLambda = Lambda.real();
  //Get the maximum
  double d = rLambda.maxCoeff();
  double Eps = posd_tol * std::abs(d);
  if(rLambda.minCoeff()<Eps){
    //Replace all eigenvalues lower than Eps
    for(int i=0;i<rLambda.rows();i++){
      if(rLambda(i,i)<Eps){
        rLambda(i,i) = Eps;
      }
    }
    //Eigen vectors
    Eigen::MatrixXcd Q = es.eigenvectors();
    Eigen::MatrixXd rQ = Q.real();
    //Get the diagonal as an array
    Eigen::ArrayXd o_diag = X.diagonal();
    //Calculate the product
    Eigen::MatrixXd rQtemp = rQ.transpose();
    for(int i=0;i<rQtemp.rows();i++){
      rQtemp.row(i) = rQtemp.row(i).array()*rLambda(i,i);
    }
    X = rQ*rQtemp;
    //Reconstruct the Diagonal
    Eigen::MatrixXd D = X.diagonal().asDiagonal();
    //Parallel Maximum
    Eigen::VectorXd vecD(o_diag.size());
    for(int i=0;i< o_diag.size();i++){
      if(o_diag(i)<Eps){
        o_diag(i) = Eps ;
      }
      vecD(i) = std::sqrt(o_diag(i)/X(i,i));
    }

    for(int r=0;r<X.rows();r++){
      for(int c=0;c<X.cols();c++){
        X(r,c) = vecD(r)*X(r,c)*vecD(c);
      }
    }
  }

  if(keepDiagonal){
    X.diagonal()=diagX0;
  }

  //Return X
return(X);
}



// [[Rcpp::export]]
Eigen::MatrixXd nearPDefinite(Eigen::MatrixXd X, int maxit=1e+6, double eigtol = 1e-06, double conv_tol = 1e-07, double posd_tol = 1e-08, bool keepDiagonal=false){
  Rcpp::Environment G = Rcpp::Environment::global_env();
  Rcpp::Environment Matrix("package:Matrix");
  Function f = Matrix["nearPD"];
  Rcpp::List res = f(X, false, keepDiagonal, true, false, true, false, false,  eigtol, conv_tol, posd_tol, 100, "I", false);
  Rcpp::NumericMatrix mat = internal::convert_using_rfunction(wrap(res[0]), "as.matrix");
  Eigen::MatrixXd Xpd=convertMatrix<Eigen::MatrixXd,Rcpp::NumericMatrix>(mat);
  return(Xpd);
}
