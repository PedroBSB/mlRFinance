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


/*
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();

  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

  matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;

  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

  matrix.conservativeResize(numRows,numCols);
}


Eigen::VectorXi compareGreat(Eigen::VectorXd d,double e1){
  Eigen::VectorXi v(d.size());
  for(int i=0;i<d.size();i++){
    if(d(i)>e1){
      v(i)=1; //True
    }
    else{
      v(i)=0; //False
    }
  }
  return(v);
}

Eigen::MatrixXd removeColumns(Eigen::MatrixXd matrix, Eigen::VectorXi colToRemove)
{
  int col = colToRemove.sum();
  int cont = 0;
  Eigen::MatrixXd matrix0(matrix.rows(),col);
  for(int c=0;c<matrix.cols();c++){
    if(colToRemove(c)==1){
      for(int r=0;r<matrix.rows();r++){
            matrix0(r,cont)=matrix(r,c);
      }
      cont=cont+1;
    }
  }
  return(matrix0);
}

void nearPDefinite(Eigen::MatrixXd &mat, int maxit, double eigtol = 1e-06){
  int n = mat.cols();
  Eigen::MatrixXd D_S = Eigen::MatrixXd::Zero(mat.rows(),mat.cols());
  Eigen::MatrixXd X =mat;
  int iter = 0;
  bool converged = false;
  double conv = std::numeric_limits<double>::infinity();
  while (iter < maxit && !converged) {
    Eigen::MatrixXd Y = X;
    Eigen::MatrixXd R = Y - D_S;
    //Getting the eigenvalues
    Eigen::EigenSolver<Eigen::MatrixXd> e(R,true);
    Eigen::MatrixXd Q = e.eigenvectors().real();
    Eigen::VectorXd d = e.eigenvalues();
    double e1 = eigtol*d(1);
    Eigen::VectorXi p = compareGreat(d,e1);
    bool test = (p.array() > 0).any();
    if (!test){
      stop("Matrix seems negative semi-definite")
    }
    Eigen::MatrixXd Q0 = removeColumns(Q,p);
    Eigen::VectorXd d0 = removeElements(d,p);
    Eigen::MatrixXd Qtemp
    X <- tcrossprod(Q0 * rep(d[p], each = nrow(Q0)), Q0)
      if (doDykstra)
        D_S <- X - R
        if (doSym)
          X <- (X + t(X))/2
        if (corr)
          diag(X) <- 1
        else if (keepDiag)
          diag(X) <- diagX0
        conv <- norm(Y - X, conv.norm.type)/norm(Y, conv.norm.type)
        iter <- iter + 1
        if (trace)
          cat(sprintf("iter %3d : #{p}=%d, ||Y-X|| / ||Y||= %11g\n",
                      iter, sum(p), conv))
        converged <- (conv <= conv.tol)
  }
  if (do2eigen || only.values) {
    e <- eigen(X, symmetric = TRUE)
    d <- e$values
    Eps <- posd.tol * abs(d[1])
    if (d[n] < Eps) {
      d[d < Eps] <- Eps
      if (!only.values) {
        Q <- e$vectors
        o.diag <- diag(X)
        X <- Q %*% (d * t(Q))
        D <- sqrt(pmax(Eps, o.diag)/diag(X))
        X[] <- D * X * rep(D, each = n)
      }
    }
  }
  return(X);
}
*/
