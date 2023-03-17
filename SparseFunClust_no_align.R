FKMSparseClustering <- function(data, x, K, m, method=c('kmea','pam','hier'), maxiter = 50){
  # data is the matrix representing the functions (dim nxp)
  # K is the number of clusters
  # x the common domain of the functions (length p)
  # m is the sparsity parameter (0 < m < \mu(D), where \mu(D) is the domain length)
  # method is the chosen clustering method ('kmea','pam','hier'), no default
  # maxiter is the maximum number of iteration (50 default)
  mu <- x[length(x)] - x[1]
  if(m > mu){
    stop("m has to be less than the measure of the domain")
  }
  qualim=c('kmea','pam','hier')
  qualem <- pmatch(method,qualim)
  w_old <- rep(1, length(x)) 
  switch(qualem,{
    k_old <- kmeans(data, K)$cluster
  },{
    k_old <- pam(data, K)$cluster
  },{
    k_old <- cutree(hclust(dist(data)),K)
  })
  b_old <- GetWCSS(data, k_old)$bcss.perfeature
  perc <- m/mu
  b_ord <- sort(b_old)
  c_star <- b_ord[ceiling(length(b_ord)*perc)]
  niter <- 1
  w <- rep(0, length(x)) 
  k <- rep(0, length(k_old))
  b <- rep(0, length(x))
  cluster_difference <- sum(abs(k_old - k))
  epsilon <- 1e-6
  w_difference <- sqrt(sum((w - w_old)^2))
  out <- list()
  while(w_difference >= epsilon && cluster_difference > 0 && niter < maxiter){
    niter <- niter + 1
    w_old <- w
    k_old <- k
    w <- GetOptimalW(b_old, c_star) 
    k <- GetOptimalClusters(data, K, w, method)
    b <- GetWCSS(data, k)$bcss.perfeature
    b_old <- b
    b_ord <- sort(b_old)
    c_star <- b_ord[ceiling(length(b_ord)*perc)]
  }
  obj <- sum(w*b)
  out <- list(W = w, CLUSTER = k, OBJ=obj, ITERATION = niter)
  return(out)
}

FKMSparseClustering.permute <- function(data, x, K, mbound = NULL, method=c('kmea','pam','hier'), nperm = 20, maxiter = 50){
  # data is the matrix representing the functions (dim nxp)
  # K is the number of clusters
  # x the common domain of the functions (length p)
  # mbound is the sparsity parameter bound, i.e., the maximal value of m to be tried;
  #       recall: 0 < m < \mu(D), where \mu(D) is the domain length
  #       if set to NULL, mboud = 60%\mu(D) (defaults to NULL)
  # method is the chosen clustering method ('kmea','pam','hier')
  # nperm is the number of permutations for the gap statistics (20 default)
  # maxiter is the maximum number of iteration (50 default)
  mu <- x[length(x)] - x[1]
  n <- dim(data)[1]
  p <- dim(data)[2]
  if(length(mbound)>0){if(mbound > mu){
    stop("m has to be less than the measure of the domain")
  }
  }else{
    mbound <- .6*mu
  }
  qualim <- seq(2*min(diff(x)),mbound,len=.1*length(x))
  GAP <- numeric(length(qualim))
  iter <- 1
  for(m in qualim){
    resTRUE <- FKMSparseClustering(data, x, K, m, method)$OBJ
    resPERM <- NULL
    for(perm in 1:nperm){
      dataperm <- data
      for(j in 1:p)dataperm[,j] <- data[sample(1:n),j]
      resPERM <- c(resPERM, FKMSparseClustering(dataperm, x, K, m, method)$OBJ)
    }
    GAP[iter] <- log(resTRUE) - mean(log(resPERM))
    iter <- iter + 1
  }
  return(list(GAP = max(GAP), m = qualim[which.max(GAP)]))
}

