## USEFUL FUNCTIONS
## 1. CER function
cer <- function(P,Q){
  if(length(P)!=length(Q)){stop('le due partizioni devono avere la stessa lunghezza')}
  cer.comp <- 0
  for(i in 1:(length(P)-1)){for(j in (i+1):length(P)){cer.comp <- cer.comp + abs((P[i]==P[j])-(Q[i]==Q[j]))}}
  cer.comp <- cer.comp/choose(length(P),2)
  return(cer.comp)
}



## DATA GENERATION 1: NO-MISALIGNMENT CASE
generate.data.FV17 <- function(n, x, paramC=0.5, plots=FALSE){
  # this function generates a set of simulated functional data in 2 clusters
  # that reproduce the examples in Simulations 2A and 2B in Floriello & Vitelli (2017)
  
  # n = number of curves
  # x = curves' domain
  # paramC = proportion of cluster overlap (default 0.5, as in Simulation 2A)
  # plots = boolean; should plots be drawn (FALSE default)
  
  # returns a list including:
  # $data = matrix (n x length(x)) with the simulated data
  # $true.partition = vector (length = n) with the true cluster assignments

  a <- 3
  b <- 0
  bpert <- .5
  c <- 2
  sd1 <- .5
  sd2 <- .25

  # means
  temp <- a-4*(1-x)*paramC/(1-paramC)
  temp[which(x<=paramC)]<-(a-4*x)[which(x<=paramC)]
  media1 <- (c*sin(c*pi*x)+a)*(a-4*x)+b
  temp2 <- rep(bpert,length(x))
  temp2[which(x>paramC)]<-(bpert*(1-x)/(1-paramC))[which(x>paramC)]
  media2 <- (c*sin(c*pi*x)+a)*temp+temp2

  # group 1
  fx <- NULL
  for(i in 1:n){
    a1 <- rnorm(1,mean=a,sd=sd1)
    b1 <- rnorm(1,mean=b,sd=sd1)
    c1 <- rnorm(1,mean=c,sd=sd2)
    fx <- cbind(fx,(c1*sin(c1*pi*x)+a1)*(a1-4*x)+b1)
  }

  # group 2
  fx2 <- NULL
  for(i in 1:n){
    a1 <- rnorm(1,mean=a,sd=sd1)
    b1 <- rnorm(1,mean=b+bpert,sd=sd1)
    c1 <- rnorm(1,mean=c,sd=sd2)
    temp <- a1-4*(1-x)*paramC/(1-paramC)
    temp[which(x<=paramC)]<-(a1-4*x)[which(x<=paramC)]
    temp2 <- rep(bpert,length(x))
    temp2[which(x>paramC)]<-(bpert*(1-x)/(1-paramC))[which(x>paramC)]
    fx2 <- cbind(fx2,(c1*sin(c1*pi*x)+a1)*temp+temp2)
  }
  data <- t(cbind(fx,fx2))
  true.partition <- c(rep(1,n),rep(2,n))

  # plots
  if(plots){

    x11()
    par(mfrow = c(1,2))
    matplot(x,cbind(media1,media2),
            type='l',lty=1,col=2:3,ylab='',main='True cluster means')

    matplot(x,t(data),type='l',lty=1,col=true.partition+1, main='Set of synthetic data')
    lines(x,media1,lwd=2)
    lines(x,media2,lwd=2)

  }

  return(list(data=data, true.partition=true.partition))
}
