## ----echo=FALSE---------------------------------------------------------------
a <- letters[1:3]
table(a, sample(a), deparse.level = 0) # dnn is c("", "")

## -----------------------------------------------------------------------------
set.seed(1998)
n <- 1000
u <- runif(n)
x <- 2/sqrt(1-u) 
hist(x, breaks=500,xlim=c(2,10),prob = TRUE, main = expression(f(x)==(8/x^3)))
y <- seq(2, 10, .01)
lines(y, 8*y^(-3),col="red")

## -----------------------------------------------------------------------------
set.seed(1998)
n<-50000
u1 <- runif(n,-1,1)
u2 <- runif(n,-1,1)
u3 <- runif(n,-1,1)
u4 <-rep(0,n)
for (i in 1:n) {
  if(abs(u3[i])>=abs(u1[i])&&abs(u3[i])>=abs(u2[i]))
    u4[i]<-u2[i]
  else u4[i]<-u3[i]
}
hist(u4, breaks=500,xlim=c(-1,1),xlab="x",prob = TRUE,main="rescaled Epanechnikov kernel")
y <- seq(-1, 1, .0002)
lines(y,0.75*(1-y^2),col="red")

## -----------------------------------------------------------------------------
set.seed(1998)
n <- 1000
u <- runif(n)
x <- 2/(1-u)^(0.25)-2 #反解
hist(x, breaks=500,xlim=c(0,6),prob = TRUE, main = expression(f(x)==64(x+2)^(-5)))
y <- seq(0, 6, .01)
lines(y,64*(y+2)^(-5),col="red" )

## -----------------------------------------------------------------------------
set.seed(1998)
k<-10000
X<-runif(k,0,pi/3)
Y<-mean(pi/3*sin(X))
Y#蒙特卡洛积分估计值#

## -----------------------------------------------------------------------------
truevalue<--cos(pi/3)+cos(0)
truevalue

## -----------------------------------------------------------------------------
set.seed(1998)
m<-100000
X<-runif(m,0,1)
Y<-mean(exp(X))
theta<-exp(X)
FANGCHA<-1/m*var(theta)
FANGCHA#使用传统MC方法的方差#
FANGCHA2<-1/(2*m)*var(exp(X[1:(k/2)])+exp((1-X)[1:(k/2)]))
FANGCHA2#使用Inverse transformation method后的方差#
(FANGCHA-FANGCHA2)/FANGCHA#降低方差的效率#
a1<-10*exp(1)-3*exp(1)^2-5
a2<-2*exp(1)-0.5*exp(1)^2-1.5
(a2-a1/2)/a2#理论效率#

## -----------------------------------------------------------------------------
x<- seq(1.01,4,.01)
y<-((x^2)*exp(-0.5*x^2))/((2*pi)^0.5)
g<-function(x){
   ((x^2)*exp(-0.5*x^2))/((2*pi)^0.5)*(x>1)
}
plot(x,y,type = "l")

## -----------------------------------------------------------------------------
set.seed(1998)
m<-10000
X1<-rnorm(m,1.5,1)
Y1<-mean(g(X1)/dnorm(X1,1.5,1))
var1<-var(g(X1)/dnorm(X1,1.5,1))
Y1
var1#这里var1是g(X)/f1(X)的方差

## -----------------------------------------------------------------------------
set.seed(1998)
m<-10000
X2<-rexp(m,1)
Y2<-mean(g(X2)/dexp(X2,1))
var2<-var(g(X2)/dexp(X2,1))
Y2
var2#这里var2是g(X)/f2(X)的方差

## -----------------------------------------------------------------------------
set.seed(1998)
m<-2000
g<-function(x){
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}
f<-function(x){
  exp(-x)/(1-exp(-1))
}
a<-rep(0,6)
a[1]<-0
for(i in 1:5){
a[i+1]<-qexp((1-exp(-1))*i/5,1)}#生成区间段点
result<-matrix(rep(0,10),ncol=2)
for (i in 1:5) {
  rn<- runif(m,0,1)
  q<- -log(1-((1-exp(-1))*(rn + i -1)/5))
  re<- g(q)/((5*(exp(-q)))/(1-exp(-1)))
  result[i,1]<- mean(re)
  result[i,2]<- var(re)
}
cat("theta的估计是:",sum(result[,1]))

cat("标准差是:",sqrt(sum(result[,2]))) 


## -----------------------------------------------------------------------------
set.seed(1998)
n<- 20
alpha<- 0.05
t<-0
s<- matrix(rep(0,20000),nrow =10000)
for (i in 1:10000) {
  x<- rlnorm(n,0,1)
  s[i,1]<- mean(log(x))+sd(log(x))*qt(alpha/2,n-1)/(n-1)^0.5
  s[i,2]<- mean(log(x))-sd(log(x))*qt(alpha/2,n-1)/(n-1)^0.5
  if(s[i,2]>0 & s[i,1]< 0) 
    t<-t+1
}
t/10000

## -----------------------------------------------------------------------------
set.seed(1998)
n<- 20
alpha<- 0.05
t<-0
re<- matrix(rep(0,20000),nrow =10000)
for (i in 1:10000) {
  x<- rchisq(n,2)
  re[i,1]<- mean(x)+sd(x)*qt(alpha/2,n-1)/sqrt(n-1)
  re[i,2]<- mean(x)-sd(x)*qt(alpha/2,n-1)/sqrt(n-1)
  if(re[i,2]>2 & re[i,1]< 2){
    t<-t+1}
}
t/10000

## ----eval=FALSE, include=FALSE------------------------------------------------
#  par(mfcol=c(2,2))
#  curve(dbeta(x, 2, 2), from = 0, to = 1)
#  curve(dbeta(x, 10, 10), from = 0, to = 1)
#  curve(dbeta(x, 30, 30), from = 0, to = 1)
#  curve(dnorm(x,0.5,0.1),from =0,to =1,col="red")

## -----------------------------------------------------------------------------
set.seed(1998)
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
canshu <- c(seq(0.5,30,0.5))
N <- length(canshu)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  e <- canshu[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rbeta(n,e,e)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(canshu, pwr, type = "l",
     xlab = bquote(canshu), ylim = c(0,0.2),col='red')
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(canshu, pwr+se, lty = 3)
lines(canshu, pwr-se, lty = 3)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  set.seed(1998)
#  par(mfcol=c(2,2))
#  curve(dt(x,2), from = -3, to = 3)
#  curve(dt(x,10), from = -3, to = 3)
#  curve(dt(x,30), from = -3, to = 3)
#  curve(dnorm(x,0,1),from =-3,to =3,col="red")

## -----------------------------------------------------------------------------
set.seed(1998)
sk <- function(x) {
#computes the sample skewness coeff.
xbar <- mean(x)
m3 <- mean((x - xbar)^3)
m2 <- mean((x - xbar)^2)
return( m3 / m2^1.5 )
}
alpha <- .1
n <- 30
m <- 2500
canshu <- c(seq(1,30,1))
N <- length(canshu)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  e <- canshu[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    x <- rt(n,e)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(canshu, pwr, type = "l",
     xlab = bquote(canshu), ylim = c(0,1),col='red')
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(canshu, pwr+se, lty = 3)
lines(canshu, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}
Ftest<- function(x,y){
  xl<- length(x)
  yl<- length(y)
  sx<- var(x)
  sy<- var(y)
  T<- sy/sx
  return(as.integer(T<qf(0.055/2,xl-1,yl-2)|T>qf(1-0.055/2,xl-1,yl-1)))
}

## -----------------------------------------------------------------------------
set.seed(1998)
compare<-function(n){
m <- 1000
sigma1 <- 1
sigma2 <- 1.5
power1 <- mean(replicate(m, expr={
x <- rnorm(n, 0, sigma1)
y <- rnorm(n, 0, sigma2)
count5test(x, y)
}))
power2 <- mean(replicate(m, expr={
x <- rnorm(n, 0, sigma1)
y <- rnorm(n, 0, sigma2)
Ftest(x, y)
}))
cat("这是样本量为n的时候count5检验的功效",power1)
cat("\n这是样本量为n的时候F检验的功效",power2)}
compare(10)#小样本情况
compare(50)#中样本情况
compare(500)#大样本情况

## ----eval=FALSE, include=FALSE------------------------------------------------
#  set.seed(1998)
#  d=1
#  m=1000
#  n <- c(10, 20, 30, 50, 100, 500)
#  pM<- 1:6
#  cv<-matrix(rep(0,2*6),ncol=2)
#  for(i in 1:6){
#     cv[i,1]<-6/n[i]*qchisq(.025,d*(d+1)*(d+2)/6)
#     cv[i,2]<-6/n[i]*qchisq(.975,d*(d+1)*(d+2)/6)
#  }
#  geneb<-function(x){
#    n<-length(x)
#    vcov<- var(x)
#    vrev<- 1/vcov
#    xbar<- mean(x)
#    b<- 0
#    for (i in 1:n) {
#      for (j in 1:n) {
#        b<-b+((x[i]-xbar)*vrev*(x[j]-xbar))^3
#      }
#    }
#    return(b/(n^2))
#  }
#  set.seed(1234)
#  for (i in 1:6) {
#    sktests<- 1:1000
#    for (j in 1:1000) {
#      x<- rnorm(n[i])
#      v<- geneb(x)
#      if(v<cv[i,1]|v>cv[i,2])
#        sktests[j]<-1 else
#          sktests[j]<-0
#    }
#    pM[i]<- mean(sktests)
#  }
#  cat("对应的模拟置信系数",pM)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  set.seed(1998)
#  alpha <- .05
#  cv<-rep(0,2)
#  cv[1]<-6/30*qchisq(.025,d*(d+1)*(d+2)/6)
#  cv[2]<-6/30*qchisq(.975,d*(d+1)*(d+2)/6)
#  n <- 30
#  m <- 1000
#  epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  set.seed(1235)
#  #critical value for the skewness test
#  for (j in 1:N) { #for each epsilon
#    e <- epsilon[j]
#    matests <- numeric(m)
#    for (i in 1:m) { #for each replicate
#      sigma <- sample(c(1, 10), replace = TRUE,
#                      size = n, prob = c(1-e, e))
#      x <- rnorm(n,0,sigma)
#      if(geneb(x)<cv[1]|geneb(x)>cv[2])
#        matests[i]<-1 else
#          matests[i]<-0
#    }
#    pwr[j] <- mean(matests)
#  }
#  #plot power vs epsilon
#  plot(epsilon, pwr, type = "l",
#       xlab = bquote(epsilon), ylim = c(0,1),col='red')
#  abline(h = .05, lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m) #add standard errors
#  lines(epsilon, pwr+se, lty = 3)
#  lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
law<-law82
n <- nrow(law)
LSAT<- law$LSAT
GPA<- law$GPA
cor<-cor(LSAT,GPA)
print (cor)
#compute the jackknife replicates
cor.jack <- numeric(n)
for (i in 1:n)
cor.jack[i] <- cor(LSAT[-i],GPA[-i])
bias <- (n - 1) * (mean(cor.jack) -cor)
print(bias) #jackknife estimate of bias

## -----------------------------------------------------------------------------
se <- sqrt((n-1) *mean((cor.jack - mean(cor.jack))^2))
print(se)#jackknife estimate of standard error

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
set.seed(1237)
air<-aircondit
meMLE<-function(x,i){
  xi<-x[i,1]
  lamda<-length(xi)/sum(xi)
  return(1/lamda)
}
obj <- boot(data =air, statistic=meMLE, R = 2000)
print(boot.ci(obj,type = c("basic", "norm", "perc","bca")))

## -----------------------------------------------------------------------------
library(boot)
thett<-function(x){
eg<-eigen(cov(x))$values
theta<-eg[1]/sum(eg)
return(theta)}
theta<-thett(scor)
n<-nrow(scor)
theta.jack <- numeric(n)
for (i in 1:n){
theta.jack[i] <-thett(scor[-i,])
}
bias <- (n - 1)* (mean(theta.jack) -theta)
bias

## -----------------------------------------------------------------------------
se <- sqrt((n-1) *mean((theta.jack - mean(theta))^2))
print(se)#jackknife estimate of standard error

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits
L1 <- lm(magnetic ~ chemical)
yhat1 <- L1$coef[1] + L1$coef[2] * a
L2 <- lm(magnetic ~ chemical + I(chemical^2))
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
L3 <- lm(log(magnetic) ~ chemical)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
L4 <- lm(log(magnetic) ~ log(chemical))
logyhat4 <- L4$coef[1] + L4$coef[2] * log(a)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- matrix(rep(0,n^2),nrow=n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:(n-1)) {
  for (m in (k+1):n) {
    y <- magnetic[-k][-(m-1)]
    x <- chemical[-k][-(m-1)]
J1 <- lm(y ~ x)
yhat11 <- J1$coef[1] + J1$coef[2] * chemical[k]
yhat12 <- J1$coef[1] + J1$coef[2] * chemical[m]
e1[k,k] <- magnetic[k] - yhat11
e1[k,m] <- magnetic[m] - yhat12
J2 <- lm(y ~ x + I(x^2))
yhat21 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
yhat22 <- J2$coef[1] + J2$coef[2] * chemical[m] +J2$coef[3] * chemical[m]^2
e2[k,k] <- magnetic[k] - yhat21
e2[k,m] <- magnetic[m] - yhat22
J3 <- lm(log(y) ~ x)
logyhat31 <- J3$coef[1] + J3$coef[2] * chemical[k]
logyhat32 <- J3$coef[1] + J3$coef[2] * chemical[m]
yhat31 <- exp(logyhat31)
yhat32 <- exp(logyhat32)
e3[k,k] <- magnetic[k] - yhat31
e3[k,m] <- magnetic[m] - yhat32
J4 <- lm(log(y) ~ log(x))
logyhat41 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
logyhat42 <- J4$coef[1] + J4$coef[2] * log(chemical[m])
yhat41 <- exp(logyhat41)
yhat42 <- exp(logyhat42)
e4[k,k] <- magnetic[k] - yhat41
e4[k,m] <- magnetic[k] - yhat41
}
}
ee1<-e1[e1!=0]
ee2<-e2[e2!=0]
ee3<-e3[e3!=0]
ee4<-e4[e4!=0]
print(c(mean(ee1^2),mean(ee2^2),mean(ee3^2),mean(ee4^2)))



## -----------------------------------------------------------------------------
set.seed(2443)
countmtest <- function(x,y,m) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > m))
}
n1 <- 200
n2 <- 300
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 10000
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:(length(z))
reps <- numeric(R) #storage for replicates
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size =length(x), replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
reps[i] <- countmtest(x1, y1,7)
}
mean(reps)



## ----eval=FALSE, include=FALSE------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  m <- 300; k<-3; p<-2; mu <- 0.3; set.seed(12346)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
#    y <- matrix(rnorm(n2*p,0,2.04),ncol=p);
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12346)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE, include=FALSE------------------------------------------------
#  m <- 300; k<-3; p<-2; mu <- 0.3; set.seed(12347)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
#    y <- matrix(rnorm(n2*p,0.5,1.8),ncol=p);
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12347)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE, include=FALSE------------------------------------------------
#  m <- 300; k<-3; p<-2; mu <- 0.3; set.seed(12348)
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rt(n1*p,1),ncol=p);
#    y1 <- matrix(rnorm(n2*p,0,1),ncol=p)
#    y2 <- matrix(rnorm(n2*p,0,10),ncol=p)
#    index<-sample(0:1,n2,replace=TRUE,c(0.5,0.5))
#    y<-matrix(rep(0,p*n2),ncol=p)
#    for (j in 1:(n2)) {
#    y[j,]<-index[j]*y1[j,]+(1-index[j])*y2[j,]
#    }
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12348)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## ----eval=FALSE, include=FALSE------------------------------------------------
#  m <- 300; k<-3; p<-2; mu <- 0.3; set.seed(12349)
#  n1 <-50;n2 <- 100; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1),ncol=p);
#    y <- matrix(rnorm(n2*p,0.1,1.3),ncol=p);
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12349)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow

## -----------------------------------------------------------------------------
fl<-function(x) 0.5*exp(-abs(x))
x<-seq(-5,5,0.01)
y<-fl(x)
plot(x,y,type='l',main='The Standard Laplace Distribution')
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (fl(y) / fl(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 20
set.seed(1278)
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)
#number of candidate points rejected
print(c(1-rw1$k/2000, 1-rw2$k/2000, 1-rw3$k/2000, 1-rw4$k/2000))#the acceptance rates
par(mfrow=c(1,1))
plot(1:length(rw2$x),rw2$x,type = 'l',xlab = 'σ=0.5',ylab = 'X')

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

## ----include=FALSE------------------------------------------------------------
normal.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for Normal(0,1)
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
xt <- x[i-1]
y <- rnorm(1, xt, sigma) #candidate point
r1 <- fl(y) * dnorm(xt, y, sigma)
r2 <- fl(x[i-1]) * dnorm(y, xt, sigma)
r <- r1 / r2
if (u[i] <= r) x[i] <- y else
x[i] <- xt
}
return(x)
}

## ----eval=FALSE, include=FALSE------------------------------------------------
#  sigma <-0.05 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 200000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#  X[i, ] <- normal.chain(sigma, n, x0[i])
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  #plot psi for the four chains
#  par(mfrow=c(1,1)) #restore default
#  #plot the sequence of R-hat statistics
#  rhat1 <- rep(0, n)
#  for (j in (b+1):n)
#  rhat1[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat1[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  sigma <-0.5 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 30000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#  X[i, ] <- normal.chain(sigma, n, x0[i])
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  #plot psi for the four chains
#  par(mfrow=c(1,1)) #restore default
#  #plot the sequence of R-hat statistics
#  rhat2 <- rep(0, n)
#  for (j in (b+1):n)
#  rhat2[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat2[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  sigma <-2 #parameter of proposal distribution
#  k <- 4 #number of chains to generate
#  n <- 10000 #length of chains
#  b <- 1000 #burn-in length
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#  X[i, ] <- normal.chain(sigma, n, x0[i])
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#  psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  #plot psi for the four chains
#  par(mfrow=c(1,1)) #restore default
#  #plot the sequence of R-hat statistics
#  rhat3 <- rep(0, n)
#  for (j in (b+1):n)
#  rhat3[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat3[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
soroot<-function(K){
N<-length(K)
root<-1:N
for (i in 1:length(K)) {
  k<-K[i]
  S1<-function(a){
   re1<-1-pt(sqrt(a^2*(k-1)/(k-a^2)),df=(k-1))
   re2<-1-pt(sqrt(a^2*k/((k+1)-a^2)),df=k)
   return(re1-re2)
   }
  root[i]<-uniroot(S1,c(0.001,2))$root
}
print(root)
}
soroot(4:25)#k=4:25
soroot(4:100)#k=4:100

## -----------------------------------------------------------------------------
nA<-444
nB<-132
nOO<-361
nAB<-63
p0<-0.3
q0<-0.2
LL<-function(theta) {
  p <- theta[1]
  q <- theta[2]
  loglik<-2*nA*((p0^2)/(p0^2+2*p0*(1-p0-q0)))* log(p)+2*nB*((q0^2)/(q0^2+2*q0*(1-p0-q0)))* log(q)+2*nOO*log(1-p-q)+nA*((2*p0*(1-p0-q0))/(p0^2+2*p0*(1-p0-q0)))*log(p*(1-p-q))+nB*((2*q0*(1-p0-q0))/(q0^2+2*q0*(1-p0-q0)))*log(q*(1-p-q))+nAB*log(p*q)
  return(-loglik)
}
BEM<-function(n){
pars<-matrix(rep(0,2*(n+1)),ncol=2)
values<-1:(n+1)
opt<-optim(c(0.3,0.2), LL)
theta0<-opt$par
pars[1,]<-theta0
values[1]<-opt$value
for(i in 1:n){
  p1<-pars[i,1]
  q1<-pars[i,2]
  LL1<-function(theta) {
    p <- theta[1]
    q <- theta[2]
    loglik<-2*nA*((p1^2)/(p1^2+2*p1*(1-p1-q1)))* log(p)+2*nB*((q1^2)/(q1^2+2*q1*(1-p1-q1)))* log(q)+2*nOO*log(1-p-q)+nA*((2*p1*(1-p1-q1))/(p1^2+2*p1*(1-p1-q1)))*log(p*(1-p-q))+nB*((2*q1*(1-p1-q1))/(q1^2+2*q1*(1-p1-q1)))*log(q*(1-p-q))+nAB*log(p*q)
    return(-loglik)
 }
opt<-optim(c(pars[i,1],pars[i,2]), LL1)
theta1<-opt$par
pars[(i+1),]<-theta1
values[i+1]<-opt$value
}
result<-cbind(pars,-values)
colnames(result)<-c('p','q','ml values')
return(result)
}
BEM(5)

## -----------------------------------------------------------------------------
mpg<-mtcars$mpg
disp<-mtcars$disp
wt<-mtcars$wt
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
for (i in 1:4) {
  print(lm(formulas[[i]]))
}

## -----------------------------------------------------------------------------
mpg<-mtcars$mpg
disp<-mtcars$disp
wt<-mtcars$wt
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
lapply(formulas,lm)


## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
sapply(trials,function(x) x$p.value)

## -----------------------------------------------------------------------------
library(parallel)
summary <- function(x) {
funs <- c(mean, median, sd, mad, IQR)
sapply(funs, function(f) f(x, na.rm = TRUE))
}
df<-data.frame(replicate(6,sample(c(1:10,NA),10,rep=T)))
round(vapply(df,summary,FUN.VALUE=c(mean=0,median=0,sd=0,mad=0,IQR=0)),3)
L1<-matrix(as.vector(unlist(Map(summary,df))),nrow=5,ncol=6)
colnames(L1)<-c('x1','x2','x3','x4','x5','x6')
rownames(L1)<-c('mean','median','sd','mad','IQR')
L1

## -----------------------------------------------------------------------------
library(Rcpp)
library(RcppArmadillo)
sourceCpp(
code = '
#include<Rmath.h>
#include<RcppCommon.h>
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;
// [[Rcpp::export]]
extern "C" SEXP rw_cpp(
double sigma,
double x0,
int N){
vec x(N, fill::zeros);
x[0] = x0;
vec u = randu<vec>(N);
double k = 0.0;
for (int i = 1; i < N; i++){
double y = ::Rf_rnorm(x(i - 1), sigma);
if ( u(i) <= exp(abs(x(i - 1))) / exp(abs(y)) ){
x(i) = y;
} else{
x(i) = x(i - 1);
k++;
}
}
double accept_rate = 1 - (k / N);
return Rcpp::List::create(
Rcpp::Named("x") = x,
Rcpp::Named("accept.rate") = accept_rate
);
}
')

## -----------------------------------------------------------------------------
fl<-function(x) 0.5*exp(-abs(x))
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (fl(y) / fl(x[i-1])))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 20
set.seed(1278)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rwc3<- rw_cpp(sigma[3], x0, N)
#number of candidate points rejected
par(mfrow=c(1,2))
plot(1:length(rw3$x),rw3$x,type = 'l',xlab = 'σ=2',ylab = 'X',main='R function')
plot(1:length(rwc3$x),rwc3$x,type = 'l',xlab = 'σ=2',ylab = 'X',main='Rcpp function')
qqnorm(rw3$x,main='R function',)
qqline(rw3$x)
qqnorm(rwc3$x,main='Rcpp function')
qqline(rwc3$x)

## -----------------------------------------------------------------------------
library('microbenchmark')
print(c(1-rw3$k/2000, rwc3$accept.rate))
compare<- microbenchmark(rw3=rw.Metropolis(sigma[3], x0, N), rwc3=rw_cpp(sigma[3], x0, N))
compare

