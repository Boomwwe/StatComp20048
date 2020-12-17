#' @title A plot function of naive kernel estimation
#' @description By the function, we can estimate a mixed normal distribution by the naive kernel estimation and draw the plot.
#' @param N the number of random number
#' @param h the bandwidth of the naive kernel estimation
#' @return a plot of a naive kernel estimation function
#' @examples
#' \dontrun{
#' nk <- naivekernelplot(1000,0.1)
#' par(mfrow=c(1,1));
#' nk
#' }
#' @export
#' @importFrom graphics plot
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20048

naivekernelplot<-function(N,h){
    rf<-1:N
    index<-sample(1:6,N,replace=TRUE,prob=c(0.5,0.1,0.1,0.1,0.1,0.1))
    for (i in 1:N) {
       if(index[i]==1){rf[i]<-rnorm(1,0,1)}
       else rf[i]<-rnorm(1,(index[i]-2)/2-1,0.1)
}
    naive<-function(x,X,h){
    n<-length(X)
    K<-function(l){
    if(-1<l&&l<1) return(1)
    else return(0)
    }
    sum=0
    for (i in 1:n) {
    sum=sum+(1/n)*(1/h)*K((X[i]-x)/h)
  }
    return(sum)
}
    xl<-seq(-3,3,0.001)
    m<-1:length(xl)
    for(i in 1:length(xl)){
    m[i]<-naive(xl[i],rf,h)
  }
    plot(xl,m,type='l',xlab ='x',ylab='naive estimate f(x)')
}
