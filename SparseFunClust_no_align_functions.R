GetWCSS <- function(x, Cs, ws=NULL){
  # x is the nxp matrix of functions
  # Cs is the cluster assignment (vector length n)
  # ws is the current w(x)
  # returns the Within Cluster Sum of Squares, for each domain point
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}

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
