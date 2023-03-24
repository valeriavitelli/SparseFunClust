GetOptimalW <- function(b, c_star){
  # b is the function b(x)
  # c_star is the parameter c* to decide what component is to be set zero
  # returns the optimal w(x)
  b_star <- b
  b_star[which(b <= c_star)] <- 0
  norm_b_star <- sqrt(sum((b_star)^2), na.rm = TRUE)
  w <- (1/norm_b_star)*b_star
  return(w)
}

GetOptimalClusters <- function(data, K, w, method){
  # data is the nxp matrix of functions
  # K is the number of clusters 
  # w is the function w(x) - vector of length p
  # method is a string defining the clustering method ('kmea','pam','hier')
  # returns the cluster assignments
  qualim=c('kmea','pam','hier')
  qualem <- pmatch(method,qualim)
  weighted_data <- scale(data,center=FALSE,scale=1/w)
  switch(qualem,{
    km <- kmeans(weighted_data, K)$cluster
  },{
    km <- pam(weighted_data, K)$cluster
  },{
    km <- cutree(hclust(dist(weighted_data)),K)
  })
  return(km)
}

GetTemplates <- function(data, clusters, w){
  # data is the n x p matrix of functions
  # clusters is the vector of cluster assignments to K clusters - vector of length n
  # w is the function w(x) - vector of length p
  # returns the K x p matrix of cluster templates
  mytmp <- colMeans(data, na.rm = TRUE)
  ind.mod <- which(w!=0)
  template <- NULL
  K <- max(unique(clusters))
  for(k in 1:K){
    ksel <- which(clusters==k)
    mytmp2 <- colMeans(data[ksel,,drop=FALSE], na.rm = TRUE)
    mytmp[ind.mod] <- mytmp2[ind.mod]
    template <- rbind(template, mytmp)
  }
  
  return(template)  
}

## function to use with L2 distance:
GetWCSSalign <- function(Y, Cs, warp, xreg, Xall){
  # Y    is the nxp matrix of functions
  # Cs   is the cluster assignment (vector length n)
  # warp is the nx2 matrix of warping functions coefficients (intercept, slope)
  # xreg is the nxp matrix of the registered abscissas
  # Xall is the vector giving the union of the functions' registered abscissas
  
  # returns the Within and Between Cluster Sum of Squares, for each domain point
  
  #manage the domain first
  n <- dim(xreg)[1]
  p <- length(Xall)
  invh.domain <- array(0, dim=c(n, n, p))
  invh.matrix <- (xreg - warp[,1])/warp[,2]
  ai <- apply(invh.matrix, 1, min, na.rm=TRUE)
  bi <- apply(invh.matrix, 1, max, na.rm=TRUE)
  for(i in 1:n)for(j in 1:n)invh.domain[i,j,which(Xall <= min(bi[i], bi[j]) & Xall >= max(ai[i],ai[j]))] <- 1
  invh.Dmeasure <- apply(invh.domain, c(1,2), integral, x=Xall)
  
  # now compute the wcss & bcss
  baseline <- numeric(ncol(Y))
  wcss.perfeature <- numeric(ncol(Y))
  for(k in unique(Cs)){
    whichers <- which(Cs==k)
    if(length(whichers)>1){
      tmp <- numeric(ncol(Y))
      for(i in 1:length(whichers)){for(j in 1:length(whichers)){
        ab <- range(Xall[which(invh.domain[whichers[i],whichers[j],]==1)],na.rm = TRUE)
        indi <- indj <- baseline
        indi[which(xreg[whichers[i],] <= ab[2] & xreg[whichers[i],] >= ab[1])] <- 1
        indj[which(xreg[whichers[j],] <= ab[2] & xreg[whichers[j],] >= ab[1])] <- 1
        tmp <- tmp + (Y[whichers[i],] - Y[whichers[j],])^2*(!is.na(invh.matrix[whichers[i],])*indi*!is.na(invh.matrix[whichers[j],]*indj)/sqrt(invh.Dmeasure[whichers[i],whichers[j]]))
      }}
      wcss.perfeature <- wcss.perfeature + tmp/(length(whichers))
    } 
  }
  tot.perfeature <- numeric(ncol(Y))
  for(i in 1:n){for(j in 1:n){
    ab <- range(Xall[which(invh.domain[i,j,]==1)],na.rm = TRUE)
    indi <- indj <- baseline
    indi[which(xreg[i,] <= ab[2] & xreg[i,] >= ab[1])] <- 1
    indj[which(xreg[j,] <= ab[2] & xreg[j,] >= ab[1])] <- 1
    tot.perfeature <- tot.perfeature + (Y[i,] - Y[j,])^2*(!is.na(invh.matrix[i,])*indi*!is.na(invh.matrix[j,]*indj)/sqrt(invh.Dmeasure[i,j]))
  }}
  tot.perfeature <- tot.perfeature/(n)
  bcss.perfeature <- tot.perfeature - wcss.perfeature
  
  # return wcss and bcss
  return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}


## function to use with H1 similarity:
GetWCSSalignRho <- function(Y, Cs, warp, xreg, Xall){
  # Y    is the nxp matrix of functions
  # Cs   is the cluster assignment (vector length n)
  # warp is the nx2 matrix of warping functions coefficients (intercept, slope)
  # xreg is the nxp matrix of the registered abscissas
  # Xall is the vector giving the union of the functions' registered abscissas
  
  # returns the Within Cluster Sum of Squares, for each domain point
  
  #manage the domain first
  n <- dim(xreg)[1]
  p <- length(Xall)
  mydim <- dim(Y)[3]
  invh.domain <- array(0, dim=c(n, n, p))
  invh.matrix <- (xreg - warp[,1])/warp[,2]
  ai <- apply(invh.matrix, 1, min, na.rm=TRUE)
  bi <- apply(invh.matrix, 1, max, na.rm=TRUE)
  for(i in 1:n)for(j in 1:n)invh.domain[i,j,which(Xall <= min(bi[i], bi[j]) & Xall >= max(ai[i],ai[j]))] <- 1
  invh.Dmeasure <- apply(invh.domain, c(1,2), integral, x=Xall)
  
  # now compute the wcss
  baseline <- numeric(p)
  wcss.perfeature <- NULL
  for(k in unique(Cs)){
    whichers <- which(Cs==k)
    if(length(whichers)>1){
      dim.tmp <- NULL
      for(d in 1:mydim){
        tmp <- NULL
        for(i in 1:length(whichers)){for(j in 1:length(whichers)){
          ab <- range(Xall[which(invh.domain[whichers[i],whichers[j],]==1)],na.rm = TRUE)
          indi <- indj <- baseline
          indi[which(xreg[whichers[i],] <= ab[2] & xreg[whichers[i],] >= ab[1])] <- 1
          indj[which(xreg[whichers[j],] <= ab[2] & xreg[whichers[j],] >= ab[1])] <- 1
          tmp <- cbind(tmp, (Y[whichers[i],,d]/L2norm(Xall,Y[whichers[i],,d])*Y[whichers[j],,d]/L2norm(Xall,Y[whichers[j],,d]))*(!is.na(invh.matrix[whichers[i],])*indi*!is.na(invh.matrix[whichers[j],])*indj))
        }}
        dim.tmp <- cbind(dim.tmp, apply(tmp, 1, sum, na.rm = TRUE)/length(whichers))
      }
      wcss.perfeature <- cbind(wcss.perfeature, apply(dim.tmp, 1, mean, na.rm = TRUE))
    } 
  }
  wcss.perfeature <- apply(wcss.perfeature, 1, sum, na.rm=TRUE)
  
  # return wcss
  return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature)))
}



### 4 functions to re-approx the registered data
matapprox.x <- function(X, xout){# abscissa
  Xnew <- matrix(NA, dim(X)[1], length(xout))
  for(i in 1:dim(X)[1])
    Xnew[i,] <- approx(X[i,], X[i,],xout,ties='ordered')$y
  return(Xnew)
}
matapprox.Y <- function(X, Y, xout){# "all-to-all" function
  Ynew <- matrix(NA, dim(X)[1], length(xout))
  for(i in 1:dim(X)[1])
    Ynew[i,] <- approx(X[i,], Y[i,],xout,ties='ordered')$y
  return(Ynew)
}
matapprox.multiY <- function(X, Y, xout){# "all-to-all" function for multidim functions
  Ynew <- array(NA, dim=c(dim(X)[1], length(xout), dim(Y)[3]))
  for(d in 1:(dim(Y)[3]))for(i in 1:(dim(X)[1]))Ynew[i,,d] <- approx(X[i,], Y[i,,d], xout, ties='ordered')$y
  return(Ynew)
}
matapprox.y <- function(X, y, xout){# "all-to-one" function
  Ynew <- matrix(NA, dim(X)[1], length(xout))
  for(i in 1:dim(X)[1])Ynew[i,] <- approx(X[i,], y,xout,ties='ordered')$y
  return(Ynew)
}

### integral approximation (trapezoid method)
integral <- function(x,y){sum((y[-1]+y[-length(y)])*diff(x)/2, na.rm=TRUE)}
### L2 norm
L2norm <- function(x,y){sqrt(integral(x,y^2))}
