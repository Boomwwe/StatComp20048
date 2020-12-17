#' @title A function to generate random number
#' @description A function to generate random number for my research. This function can control the number and the dimensions of of the random numbers. And it can generate three kinds of random numbers.
#' @param k the number of random number
#' @param l the dimensions of linear random number changed by time
#' @param fl the dimensions of log random number changed by time
#' @param e the dimensions of exponential random number changed by time
#' @return a dataframe of random numbers
#' @examples
#' \dontrun{
#' redata(10,1,1,1)
#' }
#' @export
#' @importFrom stats rnorm runif
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20048
redata<-function(k,l,fl,e){
  id<-c(1:k,1:k)
  time<-c(rep(1,k),rep(2,k))
  data<-cbind(id,time)
    linear<-matrix(rep(0,2*k*l),ncol=l)
    for(j in 1:l){
      linear1<-1:k
      linear2<-1:k
      for(i in 1:k){
        a1<-runif(1,1,2)
        linear1[i]<-a1*time[i]+rnorm(1,0,1)
        linear2[i]<-a1*time[i+k]+rnorm(1,0,1)
      }
      li<-c(linear1,linear2)
      linear[,j]<-li  
    }
    data<-cbind(data,linear)
    flog<-matrix(rep(0,2*k*fl),ncol=fl)
    for(j in 1:fl){
      flog1<-1:k
      flog2<-1:k
      for(i in 1:k){
        a1<-runif(1,1,2)
        flog1[i]<-a1*log(time[i])+rnorm(1,0,1)
        flog2[i]<-a1*log(time[i+k])+rnorm(1,0,1)
      }
      flo<-c(flog1,flog2)
      flog[,j]<-flo 
    }
    data<-cbind(data,flog)
    fexp<-matrix(rep(0,2*k*e),ncol=e)
    for(j in 1:e){
      fexp1<-1:k
      fexp2<-1:k
      for(i in 1:k){
        a1<-runif(1,1,2)
        fexp1[i]<-a1*exp(time[i])+rnorm(1,0,1)
        fexp2[i]<-a1*exp(time[i+k])+rnorm(1,0,1)
      }
      fe<-c(fexp1,fexp2)
      fexp[,j]<-fe  
    }
    data<-cbind(data,fexp)
  colnames(data)<-c("id","time",rep("linear",l),rep("log",fl),rep("exp",e))
  return(data)
}