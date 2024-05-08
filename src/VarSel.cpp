#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

// Functions in Rcpp for variational mixture model with variable selection

// [[Rcpp::export(name = ".cmatrixCalc")]]
NumericMatrix cmatrixCalc(NumericMatrix nullphi, NumericMatrix X, NumericVector c, double N, double D) {
  NumericMatrix cmatrix(N, D );
  for (int n = 0; n < N; n++){
    for(int d = 0; d < D; d++){
      cmatrix(n, d) = (1 - c(d)) * log(nullphi(d, X(n, d) - 1));
    }
  }
  return cmatrix;
}

// [[Rcpp::export(name = ".ElogthetaCalcCat")]]
NumericMatrix ElogthetaCalcCat(NumericMatrix beta, double K, double J) {
  NumericMatrix v(K, J);
  for (int k = 0; k < K; k++){
    double sum = 0; // Sum value
    for (int j = 0; j < J; j++){
      sum += beta(k, j);
    }
    for (int j = 0; j < J; j++){
      v(k, j) = boost::math::digamma(beta(k, j)) - boost::math::digamma(sum);
    }
  }
  return v;
}

// [[Rcpp::export(name = ".betaCalc")]]
NumericMatrix betaCalc(NumericVector priorbeta, NumericVector y, double K, double J, double N, NumericMatrix rnk){
  NumericMatrix v(K, J);
  for (int k = 0; k < K; k++){
    for (int j = 0; j < J; j++){
      double sum = 0; // Sum value
      for(int n = 0; n < N; n++){
        if(y(n) == j+1){
          sum += rnk(n, k);
        } 
      }
      v(k, j) = priorbeta(j) + sum;
    }
  }
  return v;
}


// [[Rcpp::export(name = ".CpriorbetaCalc")]]
NumericVector CpriorbetaCalc(NumericVector priorbeta, double K, double J){
  NumericVector v = NumericVector(Dimension(K));
  for (int k = 0; k < K; k++){
    double sum1 = 0; // Sum value
    double sum2 = 0;
    for(int j = 0; j < J; j++){
      sum1 += priorbeta(j);
      sum2 += lgamma(priorbeta(j));
    }
    v(k) = lgamma(sum1) - sum2;
  }
  return v;
}

// [[Rcpp::export(name = ".CpostbetaCalc")]]
NumericVector CpostbetaCalc(NumericMatrix beta, double K, double J){
  NumericVector v = NumericVector(Dimension(K));
  for (int k = 0; k < K; k++){
    double sum1 = 0; // Sum value
    double sum2 = 0;
    for(int j = 0; j < J; j++){
      sum1 += beta(k, j);
      sum2 += lgamma(beta(k, j));
    }
    v(k) = lgamma(sum1) - sum2;
  }
  return v;
}

// [[Rcpp::export(name = ".respthetaCalc")]]
NumericMatrix respthetaCalc(NumericMatrix Elogtheta, NumericMatrix rnk, NumericVector y, double N, double K){
  NumericMatrix v(N, K);
  for (int k = 0; k < K; k++){
    for (int n = 0; n < N; n++){
      v(n, k) = rnk(n, k) * Elogtheta(k, y(n) - 1);
    }
  }
  return v;
}

