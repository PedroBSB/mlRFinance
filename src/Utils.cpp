#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <cmath>
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
Eigen::VectorXd repetVector(Eigen::VectorXd d,double e1, int nrows){
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


