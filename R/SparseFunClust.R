#' @title Compute Sparse Functional Clustering & Alignment
#' @param data matrix representing the functions (n x p)
#' @param x matrix giving the domain of each function (n x p), or a
#'  p-dimensional vector giving the common domain
#' @param K number of clusters
#' @param template.est text string giving choices for the template estimation method
#' @param n.out number of abscissa points on which w(x) is estimated (default 500)
#' @param iter.max maximum number of iterations of the clustering loop (default 50)
#' @param vignette boolean (should the algorithm progress be reported?)
#' @param tol tolerance criterion on the weighting function to exit the loop (default 1%)
#' @param perc alignment parameter (max proportion of shift / dilation at each iter
#'   of the warping procedure) --> (default 3%)
#' @param tuning.par list of settings for the tuning of the sparsity parameter
#'  (defaults to \code{list(mbound = NULL, nperm = 20)}: mbound = max value of
#'  the sparsity parameter to be tested, default 60%; nperm = number of
#'  permutations to be performed in the tuning, default 20
#' @param tuning.m boolean (should the sparsity parameter be tuned via a
#'  permutation-based approach?)
#' @param m.prop the sparsity parameter (proportion of unrelevant domain
#'  where w(x) = 0); default 30%
#' @param clust.method the clustering method to be used; can be:
#'  'kmea' for k-means clustering,'pam','hier' for hierarchical clustering
#' @param funct.measure the functional measure to be used to compare the
#'  functions in both the clustering and alignment procedures;
#'  can be 'L2' or 'H1' (default 'L2'); see Vitelli (2019) for details
#' @param do.alignment boolean (should alignment be performed?)

#' @note
#' \code{data}:
#' \enumerate{
#'    \item assumed to be a vectorized version of the functional data AFTER smoothing
#'    \item when using the H1 functional measure, assumed to include the
#'      functions FIRST DERIVATIVES
#'    \item when using the H1 functional measure, it supports multidimensional
#'      functions R -> R^d, then data can be an array (n x p x d)]
#' }
#'
#' \code{funct.measure}: 'H1' only supported with alignment
#'
#' \code{clust.method}: 'pam' and 'hier' only supported for the case of
#'  NO ALIGNMENT
#'
#' \code{m.prop}: needs to be a proportion for compatibility with alignment,
#'  values > 1 not supported
#'
#' \code{tuning.m}: tuning only supported for the case of NO ALIGNMENT
#'
#' \code{tuning.par}:
#' \itemize{
#'    \item \code{mbound} must be lower than 1; the minimal value tested is 0
#'    \item \code{nperm} > 50 is unadvisable for computational reasons
#' }
#'
#' \code{perc}: 5% is already extreme; don't set this above 8-10%
#'
#' \code{template.est}:
#' \enumerate{
#'   \item only supported with H1 measure + ALIGNMENT
#'  \item currently 2 choices are supported:'raw' or 'loess'. 'raw' just
#'    computes the vector means across functions (default choice); 'loess'
#'    estimates the template via the R loess function
#' }
#' @return A list, with elements:
#' \describe{
#'  \item{template}{matrix (dim=K x n.out) with the final cluster templates}
#'  \item{temp.abscissa}{vector (length=n.out) of the abscissa values on which the template is defined}
#'  \item{labels}{vector (length=n) of the cluster assignments}
#'  \item{warping}{matrix (dim=n x 2) with the intercept (1st column) and slope (2nd column)
#'    of the estimated warping function for each of the n curves}
#'  \item{reg.abscissa}{matrix (dim=n x n.out) of each of the n curves registered abscissa}
#'  \item{distance}{vector (length=n) of each curve's final distance to the assigned cluster template}
#'  \item{w}{vector (length=n.out) of the estimated weighting function w(x)}
#'  \item{x.bcss}{vector (length=n.out) of the final point-wise between-cluster sum-of-squares}
#' }
#' @examples
#' set.seed(8988327)
#' x <- seq(0, 1, len = 500)
#' out <- generate.data.FV17(50, x)
#' result <- SparseFunClust(out$data, x, K = 2, do.alignment = FALSE)
#' str(result)
#' @export
SparseFunClust <- function(data, x, K, do.alignment, funct.measure = 'L2', clust.method = 'kmea',
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

  if(!is.double(K) | length(K)!=1){stop('Error: number of clusters K must be a single integer number')}

  if(do.alignment){ # joint sparse clustering WITH alignment

    # other initial checks

    if(funct.measure == 'L2' & length(dim(data))>2){
      stop('Error: array of data provided; note that multivariate functional data
            are currently  not supported when L2 functional measure is used.')
    }

    if(clust.method == 'pam' | clust.method == 'hier'){
      warning(
        'Only k-means clustering is supported when alignment is performed.',
        'Functional k-means clustering will be used instead.'
      )
    }

    if(tuning.m){
      warning(
        'Tuning of the sparsity parameter not supported when alignment is performed.
        Default value will be used instead.'
      )
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
    final.index <- apply((data - final.templates[out.no.align$CLUSTER,]), 1, function(y){L2norm(x=x,y)})
    out <- list(template=final.templates, temp.abscissa=x, labels=out.no.align$CLUSTER,
                warping=NULL, reg.abscissa=NULL, distance=final.index, w=out.no.align$W, x.bcss=final.b)
  }

  return(out)
}
