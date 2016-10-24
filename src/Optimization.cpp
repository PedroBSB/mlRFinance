#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      NUMERICAL              *********************************************/
/***************************************************************************************************************************/
//Nocedal, J., & Wright, S. (2006).
//Numerical optimization. Springer Science & Business Media. Page 453
//Solve the quadratic equality program
// Min (1/2) t(x)*G*x + t(x)*c
// Subject to:
// A*x = b
//@param G arma::mat (n x n) symetric and positive semidefinite nonlinear constants objective function
//@param c arma::vec (n x 1) linear constants objective function
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
// [[Rcpp::export]]
arma::vec quadprogEquality(arma::mat G,arma::vec c, arma::mat A, arma::vec b) {
  //Step 1.1: Create the full matrix
  arma::mat F = arma::join_rows(G, -A.t());
  //Step 1.2: Create the zero matrix
  arma::mat zero(A.n_rows,A.n_rows);
  //Step 1.2: Create A temp
  arma::mat Atemp = arma::join_rows(A, zero);
  F = join_cols(F, Atemp);
  //Step 2: Create the constant vector
  arma::vec f = arma::join_cols(-c, b);
  //Step 3.1: Solve the Linear System of Equations
  arma::vec sol = arma::solve(F, f);
  //Step 3.2: Get the solution
  int from = 0;
  int to = G.n_rows - 1;
  return sol.subvec(from, to);
}


//Identify the Working Set, i.e., Ax<=b
//@param solIni arma::vec (n x 1) Initial Feasible Solution
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
//Identify the Working Set, i.e., Ax<=b
//@param solIni arma::vec (n x 1) Initial Feasible Solution
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
arma::uvec findWorkSet(arma::vec solIni, arma::mat A,arma::vec b){
  //Initialize the vector with Active Constraints
  int iCount=0;
  //Count the number of ActiveConstraints
  //Realize the product
  arma::vec Asol=A*solIni;
  //Initialize the vector
  for(int i=0;i<Asol.n_elem;i++){
    if(Asol(i)==b(i)){
      iCount=iCount+1;
    }
  }
  //Initialize the object
  arma::uvec actCons(iCount);
  if(iCount>0){
    iCount=0;
    for(int i=0;i<Asol.n_elem;i++){
      if(Asol(i)==b(i)){
        actCons(iCount)=i;
        iCount=iCount+1;
      }
    }
  }
  //Pass only the positions
  return(actCons);
}


//@param sol current solution
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
bool testingFeasibilitySolution(arma::vec sol, arma::mat A, arma::vec b){
  bool feasible=true;
  //Calculate the product
  arma::vec Ab = A*sol;
  for(int i=0;i<Ab.n_elem;i++){
    if(Ab(i)>b(i)){
      feasible=false;
      return(feasible);
    }
  }
  return(true);
}


//@param sol current solution
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
arma::uvec findConstrNFeasibilitySolution(arma::vec sol, arma::mat A, arma::vec b){
  //Calculate the product
  arma::vec Ab = A*sol;
  //Create the count
  int cont=0;
  for(int i=0;i<Ab.n_elem;i++){
    if(Ab(i)>b(i)){
      cont=cont+1;
    }
  }
  //Intialize the new workSet
  arma::uvec newWorkSet(cont);
  cont=0;
  for(int i=0;i<Ab.n_elem;i++){
    if(Ab(i)>b(i)){
      newWorkSet(cont)=i;
      cont=cont+1;
    }
  }
  return(newWorkSet);
}

//Find all constraints with lambda greater than the negative minimum
arma::uvec verifyLambda(arma::vec lam){
  //Initialize the vector
  arma::uvec iKeep;
  double min_val = lam.min();
  if(min_val<0){
    //How many are negative ?
    int size=0;
    for(int i=0;i<lam.n_elem;i++){
      if(lam(i)>min_val){
        size=size+1;
      }
    }
    //Recreate the iKeep
    arma::uvec iKeep(size);
    int iCont=0;
    for(int i=0;i<lam.n_elem;i++){
      if(lam(i)>min_val){
        iKeep(iCont)=i;
        iCont=iCont+1;
      }
    }
    return(iKeep);
  }
  else{
    return(iKeep);
  }
}

double findMax(arma::vec ca){
  //Number of elements
  int nelem=ca.n_elem;
  //Initialize the counter
  int cont=0;
  //Initialize Maximum
  double max=0;
  //Count the number of problem elements
  for(int i=0;i<nelem;i++){
    if(!std::isnan(ca(i)) && !std::isinf(ca(i))){
      if(ca(i)>max){
        max=ca(i);
      }
    }
  }
  return(max);
}

//Solve the quadratic inequality program
// Min (1/2) t(x)*G*x + t(x)*c
// Subject to:
// A*x <= b
//@param solIni arma::vec (n x 1) Initial Feasible Solution
//@param G arma::mat (n x n) symetric and positive semidefinite nonlinear constants objective function
//@param c arma::vec (n x 1) linear constants objective function
//@param A arma::mat (m x n) jacobian of constraints (m<n)
//@param b arma::vec (m x 1) constraints constants
// [[Rcpp::export]]
arma::vec quadprogInequality(arma::vec sol, arma::mat G,arma::vec c, arma::mat A, arma::vec b) {
  //Define the maximum number of interactions
  double MAXIT=1e+5;
  double niter=0.0 ;
  //Total number of decision variables
  int nvars=G.n_cols;
  //Find the workingSet
  arma::uvec workSet = findWorkSet(sol,A,b);
  //Initialize the current solution
  arma::vec xSol;
  //Initialize lambda
  arma::vec lam;
  do{
    //workingSet with zero length
    if(workSet.n_elem>0){
      //Subset working set
      arma::mat Awork=A.rows(workSet);
      Awork=Awork.t();
      arma::vec bwork=b.elem(workSet);
      //Create KKT matrix
      // QR decomposition
      arma::mat KKT0=arma::join_rows(G,Awork);
      arma::mat zeros(Awork.n_cols,Awork.n_cols);
      zeros.fill(0.0);
      arma::mat KKT1=arma::join_rows(Awork.t(),zeros);
      arma::mat KKT= arma::join_cols(KKT0, KKT1);
      //Right Hand Side
      arma::vec rhs=arma::join_cols(-c, bwork);
      //Solve the linear system
      arma::vec solLinear=arma::solve(KKT,rhs);
      //Store the results
      xSol=solLinear.subvec(0, (nvars-1));
      lam =solLinear.subvec(nvars,(solLinear.n_elem-1));
    }
    else{
      //EQP without constraints
      xSol=arma::solve(G,-c);
    }
    //Is the solution feasible ?
    bool feasible = testingFeasibilitySolution(xSol, A, b);
    if(feasible){
      double minLambda=lam.min();
      if(minLambda<0){
        //Remove the constraint with negative lambda
        arma::uvec iKeep = verifyLambda(lam);
        if(iKeep.n_elem>0){
          workSet=workSet.elem(iKeep);
          //Update the solution
          sol=xSol;
        }
        else{
          arma::uvec null;
          workSet=null;
          sol=xSol;
        }
      }
      else{
        sol=xSol;
        return(sol);
        break;
      }
    }
    else{
      //Add the constraint that violates in workSet
      workSet=findConstrNFeasibilitySolution(xSol, A, b);
      //Fin the maximum
      //Max t, subject to:
      //Amat%*%(sol+t*(xSol-sol))<=b
      arma::vec a=A*(xSol-sol);
      arma::vec cc=b-(A*sol);
      arma::vec ca=cc/a;
      //Remove INF and NA and Find Maximum
      double t= findMax(ca);
      //Increment the solution
      sol=sol+t*(xSol-sol);
    }
    //Update the niter
    niter=niter+1;
  } while (niter<MAXIT);
  return(sol);
}

