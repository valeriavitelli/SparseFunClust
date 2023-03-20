## main function computing Sparse Functional Clustering & Alignment
 
# data = matrix representing the functions (n x p)
# NOTES: 
#   [1. assumed to be a vectorized version of the functional data AFTER smoothing]
#   [2. when using the H1 functional measure, assumed to include the functions FIRST DERIVATIVES]
#   [3. when using the H1 functional measure, it supports multidimensional functions R -> R^d;
#       then data can be an array (n x p x d)]

# x = matrix giving the domain of each function (n x p),
#   or a p-dimensional vector giving the common domain 

# K = number of clusters

# do.alignment = boolean (should alignment be performed?)

# funct.measure = the functional measure to be used to compare the functions in
#                 both the clustering and alignment procedures;
#                 can be 'L2' or 'H1' (default 'L2'); see Vitelli (2019) for details
#                 [NOTE: 'H1' only supported with alignment]

# clust.method = the clustering method to be used; can be:
#             'kmea', for k-means clustering 'pam','hier')
#             'pam', for PAM (see ..) 
#             'hier', for hierarchical clustering clustering
# [NOTE: 'pam' and 'hier' only supported for the case of NO ALIGNMENT]

# m.prop = the sparsity parameter (proportion of unrelevant domain where w(x) = 0)
#          (default 30%)
# [NOTE: needs to be a proportion for compatibility with alignment, values >1 not supported]

# tuning.m = boolean (should the sparsity parameter be tuned via a permutation-based approach?)
# [NOTE: tuning only supported for the case of NO ALIGNMENT]

# tuning.par = list of settings for the tuning of the sparsity parameter
#              list(mbound = NULL, nperm = 20)
#              mbound = max value of the sparsity parameter to be tested, default 60%
#              [NOTE: mbound must be lower than 1; the minimal value tested is 0]  
#              nperm = number of permutations to be performed in the tuning, default 20
#              [NOTE: nperm > 50 is unadvisable for computational reasons]

# perc = alignment parameter (max proportion of shift / dilation at each iter
#               of the warping procedure) --> (default 3%)
# [NOTE: 5% is already extreme; don't set this above 8-10%]

# tol = tolerance criterion on the weighting function to exit the loop (default 1%)

# template.est = text string giving choices for the template estimation method
# NOTES:
#   [1. only supported with H1 measure + ALIGNMENT]
#   [2. currently 2 choices are supported:'raw' or 'loess']
#   [   'raw' just computes the vector means across functions (default choice)]
#   [   'loess' estimates the template via the R loess function]

# n.out = number of abscissa points on which w(x) is estimated (default 500)
# iter.max = maximum number of iterations of the clustering loop (default 50)
# vignette = boolean (should the algorithm progress be reported?)

SparseFunClust <- function(data, x, K, do.alignment, funct.measure = 'L2', clust.method = 'kma',
                           m.prop=0.3, tuning.m = FALSE, tuning.par = list(mbound = NULL, nperm = 20),
                           perc=0.03, tol=0.01, template.est = 'raw',  
                           n.out=500, iter.max=50, vignette=TRUE){
  
  # preliminary checks
  n.obs <- dim(data)[1]
  n.abs <- dim(data)[2]
  
  if(length(dim(x))==0){
    if(length(x)!=n.abs)stop('Error: dimension mismatch between data and corresponding domain')
    x.reg <- matrix(x,n.obs,length(x),byrow=TRUE)
  }
  if(length(dim(x))==2){
    if(sum(dim(x)!=(dim(data)[1:2])))stop('Error: dimension mismatch between data and corresponding domain')
    x.reg <- x
  }
  
  if(!is.integer(K) | length(K)!=1){stop('Error: number of clusters K must be a single integer number')}
  
  if(do.alignment){ # joint sparse clustering WITH alignment
    
    # other initial checks
    
    if(funct.measure == 'L2' & length(dim(data))>2){
      stop('Error: array of data provided; note that multivariate functional data
            are currently  not supported when L2 functional measure is used.')
    }
    
    if(clust.method == 'pam' | clust.method == 'hier'){
      print('WARNING: Only k-means clustering is supported when alignment is performed.')
      print('Functional k-means clustering will be used instead.')
    } 
    
    if(tuning.m){
      print('WARNING: Tuning of the sparsity parameter not supported when alignment is performed. 
            Default value will be used instead.')
    }
    
    # run sparse functional clustering with alignment
    if(funct.measure == 'L2'){
      out <- sparseKMA(data, x, K, m.prop = m.prop, perc = perc, tol = tol, iter.max = iter.max, 
                                   n.out = n.out, vignette = vignette)
    }else{
      out <- sparseKMArho(data, x, K, m.prop = m.prop, perc = perc, tol = tol, template.est = template.est,  
                                   n.out = n.out, iter.max = iter.max, vignette = vignette)
    }
    
  }else{ # joint sparse clustering WITHOUT alignment
    
    if(length(dim(data))>2){
      stop('Error: array of data provided; note that multivariate functional data 
            are currently  not supported when alignment is not performed.')
    }  
    
    if(funct.measure == 'H1'){
      stop('Error: functional H1 measure not supported when alignment is not performed')
    }  
    

    if(tuning.m){
      out <- FKMSparseClustering.permute(data, x, K, mbound = tuning.par$mbound, nperm = tuning.par$nperm,
                                         method = clust.method, maxiter = iter.max)
      m.rescaled <- out$m[1]
    }else{
      mu <- diff(range(x))
      m.rescaled <- m.prop*mu
    }
    out.no.align <- FKMSparseClustering(data, x, K, m = m.rescaled, method = clust.method, maxiter = iter.max)
    # create the same output as for the case with alignment
    final.templates <- GetTemplates(data, clusters = out.no.align$CLUSTER, w = out.no.align$W)
    final.b <- GetWCSS(data, out.no.align$CLUSTER)$bcss.perfeature
    final.index <- apply((data - final.templates[out.no.align$CLUSTER,]), 1, function(y){L2norm(x,y)})
    out <- list(template=final.templates, temp.abscissa=x, labels=out.no.align$CLUSTER, 
                warping=NULL, reg.abscissa=NULL, distance=final.index, w=out.no.align$W, x.bcss=final.b)
  }
 
  return(out)
   
}