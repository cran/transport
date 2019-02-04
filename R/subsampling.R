#' Approximate Computation of Wasserstein Distances via Subsampling.
#'
#' Samples \code{S} elements each of a source and a target measure and
#' computes the Wasserstein distance between the samples.
#' The mean distance out of \code{K} tries is returned.
#' 
#' @param source The source measure has to be either a weight vector or an object
#'        of one of the classes \code{"pgrid"}, \code{"wpp"} or \code{"pp"}.
#' @param target The target measure needs to be of the same type as the source measure.
#' @param S The sample size.
#' @param K The number of tries. Defaults to 1.
#' @param costM The cost matrix between the source and target measures. Ignored unless source
#'        and target are weight vectors.
#' @param p The order of the Wasserstein metric (i.e. the power of the distances). Defaults to 1.
#' @param prob logical. Should the objects a, b be interpreted as probability measures, i.e. their
#'        total mass be normalized to 1?        
#' @param method A string with the name of the method used for optimal transport distance computation.
#'        Options are "revsimplex", "shortsimplex" and "primaldual". Defaults to "revsimplex".
#'    
#' @details
#' For the classes \code{"pgrid"}, \code{"wpp"} and \code{"pp"}, currently the  whole distance matrix
#' is precomputed. Then for each of the \code{K} samples the function \code{\link{wasserstein}} is
#' called using the appropriate submatrix. The computation of the distance matrix makes this function
#' slow on large measures. This will be fixed in a later version.
#'                
#' @return The mean of the K values of the Wasserstein distances between
#'         the subsampled measures.
#'         
#' @author Jörn Schrieber \email{joern.schrieber-1@mathematik.uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{dominic.schuhmacher@mathematik.uni-goettingen.de}
#' 
#' @references M. Sommerfeld, J. Schrieber, Y. Zemel and A. Munk (2018)
#'             Optimal Transport: Fast Probabilistic Approximation with Exact Solvers
#'             preprint: \href{https://arxiv.org/abs/1802.05570}{arXiv:1802.05570}
#' 
#' @examples 
#' \dontrun{
#' subwasserstein(random64a, random64b, S=1000)
#' wasserstein(random64a, random64b)
#' }
#' 
#' @export

subwasserstein <- function(source, target, S, K = 1, p = 1, costM = NULL, prob=TRUE, method = "revsimplex") {
  
  # Get mass vectors and the cost matrix.
  # This depends on the object class of source and target.
  if (is(source, "pgrid") && is(target, "pgrid")) {
    
    # generate the grid points in matrix form from the pgrid generators
    x <- expand.grid(source$generator)
    y <- expand.grid(target$generator)
    
    # compute the cost matrix
    n <- source$N
    m <- target$N
    fulldist <- as.matrix(dist(rbind(x,y)))
    costM <- fulldist[1:n,((n+1):(n+m))]
    
    # fetch weight vectors
    massA <- as.vector(source$mass)
    massB <- as.vector(target$mass)
    
  } else if (is(source, "pp") && is(target, "pp")) {
    
    # compute the cost matrix
    n <- source$N
    m <- target$N
    fulldist <- as.matrix(dist(rbind(source$coordinates,target$coordinates)))
    costM <- fulldist[1:n,((n+1):(n+m))]
    
    # set weight vectors to ones
    massA <- rep(1,source$N)
    massB <- rep(1,target$N)
    
  } else if (is(source, "wpp") && is(target, "wpp")) {
    
    # compute the cost matrix
    n <- source$N
    m <- target$N
    fulldist <- as.matrix(dist(rbind(source$coordinates,target$coordinates)))
    costM <- fulldist[1:n,((n+1):(n+m))]
    
    # fetch weight vectors
    massA <- source$mass
    massB <- target$mass
    
  } else if (is(source, "numeric") && is(target, "numeric")) {
    
    stopifnot(all.equal(dim(costM), c(length(source),length(target))))
    massA <- source
    massB <- target
    
  } else {
    stop("source measure is type ", class(source), ", target measure is type ", class(target))
  }
  
  # Compute the total mass for later scaling and stop if measures have different total masses.
  totalMass <- sum(massA)
  stopifnot(totalMass == sum(massB))
  
  res <- rep(0,K)
  # Subsample and compute the optimal transport K times.
  for (k in 1:K) {
    
    # sampling
    sourcesub <- rmultinom(1, size = S, prob = massA)
    targetsub <- rmultinom(1, size = S, prob = massB)
    
    # compute Wasserstein distance by the appropriate method;
    # actually this is typically a case for the auction algorithm -> TO DO
    res[k] <- wasserstein(sourcesub[sourcesub != 0], targetsub[targetsub != 0], p=p, tplan=NULL,
                          costm = costM[sourcesub != 0, targetsub != 0], prob=prob, method = method)
  }
  # return the mean of the results, scaled by the appropriate factor
  if (!prob) {
    res <- res * (totalMass/S)^(1/p)
  }
  return(mean(res))
}