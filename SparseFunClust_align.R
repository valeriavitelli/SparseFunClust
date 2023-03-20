### MAIN sparseKMA FUNCTION (with the normalized L^2 distance)
sparseKMA <- function(data, x, K, m.prop = .3, perc=0.03, tol=0.01, iter.max=50, 
                      n.out = 500, vignette=TRUE){
  
  # data is the matrix representing the functions (n x p)
  #   [I assume to have a vectorized version of the data AFTER smoothing]
  # x is a (n x p) matrix giving the domain of each function,
  #   or a p-dimensional vector giving the common domain 
  # K is the number of clusters
  # m.prop is the sparsity parameter (proportion of unrelevant domain 
  #               where w(x) = 0) --> (default 30%)
  # perc is the alignment parameter (max proportion of shift / dilation at each iter
  #               of the warping procedure) --> (3% default)
  #             WARNING: 5% is already extreme; don't set this above 8-10%
  # tol is the tolerance criterion on the weighting function to exit the loop (1% default)
  # iter.max is the maximum number of iteration (50 default)
  # n.out is the number of abscissa points on which w(x) is estimated
  # vignette is a boolean (should the algorithm progress be reported?)
  
  # preliminaries
  n.obs <- dim(data)[1]
  n.abs <- dim(data)[2]
  if(length(dim(x))==0){
    if(length(x)!=n.abs)print('Warning: dimension mismatch between data and corresponding domain')
    x.reg <- matrix(x,n.obs,length(x),byrow=TRUE)
  }
  if(length(dim(x))==2){
    if(sum(dim(x)!=(dim(data)[1:2])))print('Warning: dimension mismatch between data and corresponding domain')
    x.reg <- x
  }
  index <- rep(1000,n.obs)
  index.old <- rep(5000,n.obs)
  labels <- sample.int(K,n.obs,replace=TRUE)
  labels.old <- sample.int(K,n.obs,replace=TRUE)
  grid.shift <- seq(-perc,perc,len=9)
  grid.dil <- 1 + seq(-perc,perc,len=9)
  grid.val <- t(expand.grid(grid.shift, grid.dil))
  xout <- seq(min(x.reg), max(x.reg), len=n.out)
  w <- rep(1, length(xout))
  w.old <- rep(10, length(xout))
  b.old <- b <- rep(1, length(x))
  iter <- 0
  
  # initialize the template
  mytmp <- colMeans(matapprox.Y(x.reg, data, xout), na.rm = TRUE)
  template <- NULL
  for(k in 1:K)template <- rbind(template, mytmp)
  
  while( (mean(abs(index - index.old)) > tol | sum(abs(labels - labels.old)) > 0 ) & iter < iter.max ){
    
    iter <- iter + 1
    index.old <- index
    labels.old <- labels
    b.old <- b
    xout.old <- xout
    xout <- seq(min(x.reg), max(x.reg), len=n.out)
    if(iter>1)w <- as.vector(matapprox.y(rbind(xout.old), w, xout))
    w.old <- w
    #this is not necessary for the algorithm:
    #D <- diff(range(xout))
    #m <- D*m.prop
    
    # alignment step
    warp.def <- matrix(NA,n.obs,2)
    for(i in (1:n.obs)){
      cand <- scale(t(x.reg[i,,drop=FALSE])%*%grid.val[2,,drop=FALSE], center=-grid.val[1,], scale=FALSE)
      datir <- matapprox.y(t(cand),data[i,],xout=xout)# rows are possible warpings, columns are time instances
      index.temp <- NULL
      for(k in 1:K){# normalized L^2 norm between functions, with sparsity weighting function
        matrice.norma <- scale(scale(datir, center=template[k,], scale=FALSE)^2, center = FALSE, scale = 1/w)
        normak <- apply(matrice.norma,1,integral,x=xout)/(diff(range(xout))*rowSums(!is.na(matrice.norma))/dim(matrice.norma)[2])
        # index.temp: matrix with normalized L2 norm to the template, for every possible warping (rows) and cluster (columns)
        index.temp <- cbind(index.temp,  normak)
      }
      quale <- which(index.temp==min(index.temp), arr.ind=TRUE)
      quale <- quale[sample(dim(quale)[1],1),]
      index[i] <- index.temp[quale[1],quale[2]]
      warp.def[i,] <- t(grid.val[,quale[1]])
      labels[i] <- quale[2]
    }
    
    # within-cluster normalization step
    for(k in 1:K){
      ind <- which(labels == k)
      warp.def[ind,1] <- (warp.def[ind,1] - colMeans(warp.def[ind,])[1])/colMeans(warp.def[ind,])[2]
      warp.def[ind,2] <- warp.def[ind,2]/colMeans(warp.def[ind,])[2]
    }
    x.reg <- x.reg*warp.def[,2]+warp.def[,1]
    
    # update the weighting function
    b <- GetWCSSalign(matapprox.Y(x.reg, data, xout), labels, warp.def, matapprox.x(x.reg, xout), xout)$bcss.perfeature
    b.ord <- sort(b)
    c.star <- b.ord[ceiling(length(b.ord)*m.prop)]
    w <- GetOptimalW(b, c.star)
    
    # update the templates
    mytmp <- colMeans(matapprox.Y(x.reg, data, xout), na.rm = TRUE)
    ind.mod <- which(w!=0)
    template <- NULL
    for(k in 1:K){
      ksel <- which(labels==k)
      mytmp2 <- colMeans(matapprox.Y(x.reg[ksel,,drop=FALSE], data[ksel,,drop=FALSE], xout), na.rm = TRUE)
      mytmp[ind.mod] <- mytmp2[ind.mod]
      #mytmp2[which(colSums(!is.na(mytmp)) < 2)] <- NA
      template <- rbind(template, mytmp)
    }
    
    if(vignette)print(paste('Iteration: ',iter,', mean distance: ',mean(index),sep=''))
    
  }
  
  return(list(template=template, temp.abscissa=xout, labels=labels, warping=warp.def, reg.abscissa=x.reg, distance=index, w=w, x.bcss=b))
}


### MAIN sparseKMA FUNCTION (with the \rho(f1^\prime, f2^\prime) similarity of Sangalli et al., 2010)
sparseKMArho <- function(data, x, K, m.prop=0.3, perc=0.03, tol=0.01, template.est = 'raw',  
                         n.out=500, iter.max=50, vignette=TRUE){
  
  # data is the matrix representing the functions (n x p)
  # NOTE: data can be an array (n x p x d) in case of multidimensional functions R -> R^d 
  #   [I assume to have a vectorized version of the data AFTER smoothing]
  #   [here I also assume that these are the functions FIRST DERIVATIVES]
  #   [because the similarity index is ONLY based on derivatives, so no reasons to treat the original functions]
  # x is a (n x p) matrix giving the domain of each function,
  #   or a p-dimensional vector giving the common domain 
  # K is the number of clusters
  # m.prop is the sparsity parameter (proportion of unrelevant domain 
  #               where w(x) = 0) --> (default 30%)
  # perc is the alignment parameter (max proportion of shift / dilation at each iter
  #               of the warping procedure) --> (3% default)
  #             WARNING: 5% is already extreme; don't set this above 8-10%
  # tol is the tolerance criterion on the weighting function to exit the loop (1% default)
  # template.est is a text string giving choices for the template estimation method
  #   [currently 2 choices are supported:'raw' or 'loess']
  #   ['raw' just computes the vector means across functions]
  #   ['loess' estimates the template via the R loess function]
  # n.out is the number of abscissa points on which w(x) is estimated
  # iter.max is the maximum number of iteration (50 default)
  # vignette is a boolean (should the algorithm progress be reported?)
  
  # preliminaries
  n.obs <- dim(data)[1]
  n.abs <- dim(data)[2]
  if(length(dim(data))<3){
    data <- array(data, dim=c(n.obs,n.abs,1))
  }
  n.dim <- dim(data)[3]
  if(length(dim(x))==0){
    if(length(x)!=n.abs)print('Warning: dimension mismatch between data and corresponding domain')
    x.reg <- matrix(x,n.obs,length(x),byrow=TRUE)
  }
  if(length(dim(x))==2){
    if(sum(dim(x)!=(dim(data)[1:2])))print('Warning: dimension mismatch between data and corresponding domain')
    x.reg <- x
  }
  index <- rep(1000,n.obs)
  index.old <- rep(5000,n.obs)
  n.opt <- 9 # keep value small for computational reasons
  labels <- sample.int(K,n.obs,replace=TRUE)
  labels.old <- sample.int(K,n.obs,replace=TRUE)
  grid.shift <- seq(-perc,perc,len=n.opt)
  grid.dil <- 1 + seq(-perc,perc,len=n.opt)
  grid.val <- t(expand.grid(grid.shift, grid.dil))
  xout <- seq(min(x.reg, na.rm = TRUE), max(x.reg, na.rm = TRUE), len=n.out)
  w <- rep(1, length(xout))
  w.old <- rep(10, length(xout))
  b.old <- b <- rep(1, length(x))
  iter <- 0
  if(template.est == 'loess')loess.span <- 0.15
  
  # initialize the template
  template <- array(NA, dim=c(K, n.out, n.dim))
  for(d in 1:n.dim){
    if(template.est == 'raw'){
      mytmp <- colMeans(matapprox.Y(x.reg, data[,,d], xout), na.rm = TRUE)
    }else if(template.est == 'loess'){
      loess.temp <- loess(as.numeric(data[,,d])~as.numeric(x.reg), span=loess.span)
      mytmp <- predict(loess.temp,xout)
    }
    for(k in 1:K)template[k,,d] <- mytmp
  }
  
  while( mean(abs(index - index.old)/index.old) > tol & sum(abs(labels - labels.old)) > 0 & iter < iter.max ){
    
    iter <- iter + 1
    index.old <- index
    labels.old <- labels
    b.old <- b
    xout.old <- xout
    xout <- seq(min(x.reg, na.rm = TRUE), max(x.reg, na.rm = TRUE), len=n.out)
    if(iter>1)w <- as.vector(matapprox.y(rbind(xout.old), w, xout))
    w.old <- w
    #this is not necessary for the algorithm:
    #D <- diff(range(xout))
    #m <- D*m.prop
    
    # alignment step
    warp.def <- matrix(NA,n.obs,2)
    for(i in (1:n.obs)){# loop over samples
      cand <- scale(t(x.reg[i,,drop=FALSE])%*%grid.val[2,,drop=FALSE], center=-grid.val[1,], scale=FALSE)
      index.temp <- array(NA, dim=c(n.opt^2, K, n.dim))
      for(d in 1:n.dim){# loop over dimensions
        datir <- matapprox.y(t(cand),data[i,,d],xout=xout)# rows are possible warpings, columns are abscissa instances
        datirnorm <- apply(datir, 1, L2norm, x=xout)# compute also the normalization term of each warped function
        for(k in 1:K){# loop over clusters
          matrice.norma <- scale(datir, center=FALSE, scale=L2norm(xout, template[k,,d])/(template[k,,d]*w))/datirnorm # similarity index rho between functions, with sparsity weighting function
          normak <- apply(matrice.norma,1,integral,x=xout)
          index.temp[,k,d] <- normak
          # index.temp: array with rho similarity to the template, for every possible warping (rows) and cluster (columns)
          #             the third dimension is the functions' original dimension
        }
      }
      # average the similarity index over dimensions
      index.temp <- apply(index.temp, c(1,2), mean)
      # maximize the similarity to the template
      quale <- which(index.temp==max(index.temp), arr.ind=TRUE)
      quale <- quale[sample(dim(quale)[1],1),]
      index[i] <- index.temp[quale[1],quale[2]]
      warp.def[i,] <- t(grid.val[,quale[1]])
      labels[i] <- quale[2]
    }
    
    # normalization step:
    # warping functions must have identity mean within each cluster
    for(k in 1:K){
      ind <- which(labels == k)
      warp.def[ind,1] <- (warp.def[ind,1] - colMeans(warp.def[ind,,drop=FALSE])[1])/colMeans(warp.def[ind,,drop=FALSE])[2]
      warp.def[ind,2] <- warp.def[ind,2]/colMeans(warp.def[ind,,drop=FALSE])[2]
    }
    x.reg <- x.reg*warp.def[,2]+warp.def[,1]
    
    # update the weighting function:
    # w(x) is non-zero on the portion of the domain with maximal within-cluster similarity
    b <- GetWCSSalignRho(matapprox.multiY(x.reg, data, xout), labels, warp.def, matapprox.x(x.reg, xout), xout)$wcss.perfeature
    b.ord <- sort(b)
    c.star <- b.ord[ceiling(length(b.ord)*m.prop)]
    w <- GetOptimalW(b, c.star)
    
    # update the templates
    ind.mod <- which(w!=0)
    template <- array(NA, dim=c(K, n.out, n.dim))
    for(k in 1:K){
      ksel <- which(labels==k)
      if(length(ksel) > 0){
        if(template.est == 'raw'){
          mytmp2 <- apply(matapprox.multiY(x.reg[ksel,,drop=FALSE], data[ksel,,,drop=FALSE], xout), c(2,3), mean, na.rm = TRUE)
        }else if(template.est == 'loess'){
          mytmp2 <- matrix(NA, n.out, n.dim)
          for(d in 1:n.dim){
            loess.temp <- loess(as.numeric(data[ksel,,d,drop=FALSE])~as.numeric(x.reg[ksel,,drop=FALSE]), span=loess.span)
            mytmp2[,d] <- predict(loess.temp,xout)
          }
        }
      }else{
        pickatrandom <- sample(1:n.obs,1)
        mytmp2 <- matrix(NA, n.out, n.dim)
        for(d in 1:n.dim)mytmp2[,d] <- approx(x.reg[pickatrandom,], data[pickatrandom,,d], xout,ties='ordered')$y
        labels[pickatrandom] <- k
      }
      template[k,,] <- mytmp2
    }
    
    if(vignette)print(paste('Iteration: ',iter,', mean distance: ',mean(index),sep=''))
    
  }
  
  return(list(template=template, temp.abscissa=xout, labels=labels, warping=warp.def, reg.abscissa=x.reg, 
              distance=index, w=w, x.bcss=b))
}

