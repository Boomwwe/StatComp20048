#include <Rcpp.h>
using namespace Rcpp;

//' @title A function of Fibonacci sequence using Rcpp
//' @description A function to calculate the nth term of the Fibonacci sequence
//' @param n the nth term of the Fibonacci sequence
//' @return a int number
//' @examples
//' \dontrun{
//' fib1 <- fib_cpp_1(10)
//' }
//' @export
//[[Rcpp::export]]
int fib_cpp_1(int n)
{
  if(n==1||n==2) return 1;
  return fib_cpp_1(n-1)+fib_cpp_1(n-2);
}