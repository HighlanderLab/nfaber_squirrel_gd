// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat fastDist(const arma::mat& X){
  arma::colvec Xn = sum(square(X),1);
  arma::mat D = -2*(X*X.t());
  D.each_col() += Xn;
  D.each_row() += Xn.t();
  D = sqrt(D);
  D.diag().zeros(); //Removes NaN values
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0);
  }
  return D;
}

// [[Rcpp::export]]
arma::mat fastPairDist(const arma::mat& X, const arma::mat& Y){
  arma::colvec Xn = sum(square(X),1);
  arma::colvec Yn = sum(square(Y),1);
  arma::mat D = -2*(X*Y.t());
  D.each_col() += Xn;
  D.each_row() += Yn.t();
  D = sqrt(D);
  if(D.has_nan()){
    D.elem(find_nonfinite(D)).fill(0.0);
  }
  return D;
}


