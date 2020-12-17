## -----------------------------------------------------------------------------
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
 ##Example
 redata(10,2,2,2)

## -----------------------------------------------------------------------------
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
##Example:
  set.seed(1998)
  naivekernelplot(1000,0.1)

