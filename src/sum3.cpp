#include <Rcpp.h>
using namespace Rcpp;
//' @title A sum function using Rcpp
//' @description A function to calculate the sum of a vector
//' @param x the numeric vector
//' @return sum
//' @examples
//' \dontrun{
//' rnC <- sum3(1:10)
//' }
//' @export
// [[Rcpp::export]]
double sum3(NumericVector x) {
  double total = 0;
  
  NumericVector::iterator it;
  for(it = x.begin(); it != x.end(); ++it) {
    total += *it;
  }
  return total;
}