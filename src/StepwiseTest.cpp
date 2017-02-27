#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//https://github.com/cran/StepwiseTest/tree/master/src
//A Generalized Stepwise Procedure with Improved Power for Multiple Inequalities Testing



arma::vec HF(arma::vec test_stat, double cv) {
  // Decision function: 1(test_stat>cv)
  arma::vec y = arma::zeros(test_stat.n_rows);
  y.rows(arma::find(test_stat > cv)).ones();
  return y;
}


// [[Rcpp::export]]
Rcpp::List FWERkControl(arma::vec test_stat,arma::mat boot_stat,int k, double alpha) {
  // test_stat: m x 1 column vector of test statistics
  // boot_stat: m x B matrix of bootstrap statistics
  // k: Number of false rejections
  // alpha: FWER(k) level
  int m = test_stat.n_rows;
  int B = boot_stat.n_cols;
  int num_reject=0;
  int num_reject1=-m;
  int q = floor(alpha*B);

  arma::vec reject=arma::zeros(m);
  double CV;

  arma::uvec test_index = arma::stable_sort_index(test_stat,"descend"); // ranking index based on test statistics
  arma::mat ranked_boot_stat = boot_stat.rows(test_index);  // each column of boot_stat is ranked based on the test_index order
  while(num_reject>num_reject1){
    num_reject1=num_reject;
    if(num_reject < k){
      arma::mat sim_CV = arma::sort(ranked_boot_stat,"descend");
      arma::vec k_sim_max = sim_CV.row( k-1 ).t();
      arma::vec sort_k_sim_max = arma::sort(k_sim_max,"descend");
      CV=arma::as_scalar(sort_k_sim_max.row(q));
      CV=(std::abs(CV)+CV)/2;
    }
    else{
      // (k-1) least significant models + not rejected models
      arma::uvec retained_m = arma::linspace<arma::uvec>(num_reject-k+1, m-1, m-num_reject+k-1);
      arma::mat sim_CV=sort(ranked_boot_stat.rows(retained_m),"descend");
      arma::vec k_sim_max = sim_CV.row(k-1).t();
      arma::vec sort_k_sim_max=sort(k_sim_max,"descend");
      CV=arma::as_scalar(sort_k_sim_max.row(q));
    }

    reject = HF(test_stat,CV); // this gives the models that are rejected after this step
    num_reject=arma::sum(reject);
  } // end of critical value calculation (while loop)

  return Rcpp::List::create(
    Rcpp::Named("Reject") = reject.t(),
    Rcpp::Named("CV") = CV
  ) ;
}


// [[Rcpp::export]]
Rcpp::List FDPControl(arma::vec test_stat,arma::mat boot_stat, double gamma, double alpha) {
  // test_stat: m x 1 column vector of test statistics
  // boot_stat: m x B matrix of bootstrap statistics
  // gamma: The false discovery proportion (FDP) parameter
  // alpha: The bound for the probability of FDP
  int k = 1;
  Rcpp::List res = FWERkControl(test_stat,boot_stat,k,alpha);
  Rcpp::NumericVector Rejected = res[0];
  double Nk = Rcpp::sum(Rejected);
  while( k/(Nk+1) <= gamma ) {
    k = k + 1;
    res = FWERkControl(test_stat,boot_stat,k,alpha);
    Rejected = res[0];
    Nk = sum(Rejected);
  }
  Rcpp::NumericVector CV = res[1];
  return Rcpp::List::create(
    Rcpp::Named("Reject") = Rejected,
    Rcpp::Named("k_stopped") = k,
    Rcpp::Named("CV") = CV
  ) ;
}

