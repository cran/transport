\name{transport}
\alias{transport}
\alias{transport.default}
\alias{transport.pgrid}
\alias{transport.pp}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find Optimal Transference Plan (a.k.a. Optimal Transport Map) Between Two Objects.
}
\description{
  Given two objects \code{a} and \code{b} that specify distributions of mass and an object that specifies (a way to compute) costs,
  find the transference plan for going from \code{a} to \code{b} that minimizes the total cost.
}
\usage{
transport(a, b, ...)
\method{transport}{default}(a, b, costm, method = c("shortsimplex", "revsimplex", "primaldual"),
  control = list(), ...)
\method{transport}{pgrid}(a, b, p = NULL, method = c("shortsimplex", "revsimplex",
  "primaldual", "aha"), control = list(), ...)
\method{transport}{pp}(a, b, p = 1, method = c("auction", "auctionbf", "shortsimplex",
  "revsimplex", "primaldual"), control = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, b}{
    two objects that describe mass distributions, between which the optimal transport map is to be computed. For the default
    method these are vectors of non-negative values. For the other two methods these are objects of the respective classes.
}
  \item{costm}{
    for the default method a \code{length(a)} by \code{length(b)} matrix specifying the cost of transporting single
    units of mass between the corresponding source and destination points. 
}
  \item{p}{
    for the two specialized methods the power \eqn{>0} to which the Euclidean distance between points is taken in order
    to compute costs.
}
  \item{method}{
    the name of the algorithm to use. See details below.
}
  \item{control}{
  	a named list of parameters for the chosen method or the result of a call to \code{\link{trcontrol}}. Any parameters
  	that are not set by the control argument will get reasonable (sometimes problem specific) defaults.
  }
  \item{...}{
  	currently without effect.
  }
}
\details{
  There is a number of algorithms that are currently implemented and more will be added in future versions of the package. 
  The following is a brief description of each abbreviation used. Much more details can be found in the cited references 
  and in a forthcoming package vignette.
  
  \code{shortsimplex}: The shortlist method based an a revised simplex algorithm, as described in Gottschlich and Schuhmacher (2014).
  
  \code{revsimplex}: The revised simplex algorithm as described in Luenberger and Ye (2008, Section 6.4) with various speed improvements, including
  a multiscale approach.
  
  \code{primaldual}: The primal-dual algorithm as described in Luenberger (2003, Section 5.9).
  
  \code{aha}: The Aurenhammer--Hoffmann--Aronov (1998) method with the multiscale approach presented in Mérigot (2011). 
  
  \code{auction}: The auction algorithm by Bertsekas (1988) with epsilon-scaling, see Bertsekas (1992).
  
  \code{auctionbf}: A refined auction algorithm that combines forward and revers auction, see Bertsekas (1992).
  
  The following table gives information about the applicability of the various algorithms (or sometimes rather
  their current implementations).
  
  \tabular{llll}{
  	             \tab default \tab pgrid \tab pp    \cr
    shortsimplex \tab   +     \tab   +   \tab   *   \cr  
    revsimplex   \tab   +     \tab   +   \tab   *   \cr  
    primaldual   \tab   *     \tab   *   \tab   *   \cr
    aha (p=2!)   \tab   -     \tab   +   \tab   @   \cr
    auction      \tab   -     \tab   -   \tab   +   \cr
    auctionbf    \tab   -     \tab   -   \tab   +   \cr
  } 
  where: + recommended, * applicable (usually slow), @ planned, - no implementation planned / combination does not make sense.

  Each algorithm has certain parameters supplied by the \code{control} argument. The following table gives an overview of parameter names and 
  their applicability.
  
  \tabular{llll}{
  	                        \tab \code{start} \tab multiscale \tab individual parameters    \cr
    \code{shortsimplex}     \tab   -     \tab   -   \tab   \code{slength}, \code{kfound}, \code{psearched}   \cr  
    \code{revsimplex}       \tab   +     \tab   +   \tab      \cr  
    \code{primaldual}       \tab   -     \tab   -   \tab      \cr
    \code{aha} (\eqn{p=2}!) \tab   -     \tab   +   \tab   \code{factr}, \code{maxit}  \cr
    \code{auction}          \tab   -     \tab   -   \tab   \code{lasteps}, \code{epsfac}   \cr
    \code{auctionbf}        \tab   -     \tab   -   \tab   \code{lasteps}, \code{epsfac}   \cr
  } 
 
  \code{start} specifies an algorithm for computing a starting solution. Currently only the North-West Corner Rule (\code{start="nwcorner"}) and the method by Russell (1969) (\code{start="russell"}) are implemented.
  
  The word \dQuote{multiscale} in the above table summarizes the three parameter \code{nscales}, \code{scmult}, \code{returncoarse}
  that control a possible multiscale version of the corresponding algorithm. The default value of \code{nscales=1} suppresses
  the multiscale version. For larger problems it is advisable to use the multiscale version, which currently is only implemented for
  square pgrids in two dimensions. The algorithm proceeds then by coarsening the pgrid \code{nscales-1} times, summarizing
  each time \code{scmult^2} pixels into one larger pixels, and then solving the various transport problems starting from the coarsest and
  using each previous problem to compute a starting solution to the next finer problem. If \code{returncoarse} is \code{TRUE}, the coarser
  problems and their solutions are returned as well (\code{revsimplex} only).
  
  \code{slength}, \code{kfound}, \code{psearched} are the shortlist length, the number of pivot candidates needed, and the percentage of
  shortlists searched, respectively.
  
  \code{factr}, \code{maxit} are the corresponding components of the \code{control} argument in the \code{\link{optim}} L-BFGS-B method.
  
  \code{lasteps}, \code{epsfac} are parameters used for epsilon scaling in the auction algortihm. The algorithm starts with a \dQuote{transaction cost} per bid of \code{epsfac^k * lasteps} for some reasonable \code{k} generating finer and finer approximate solutions as the \code{k} counts down to zero. Note that in order for the procedure to make sense, \code{epsfac} should be larger than one (typically two- to three-digit) and in order for the final solution to be exact \code{lasteps} should be smaller than \code{1/n}, where \code{n} is the total number of points in either of the point patterns.
}
\value{
A data frame with columns \code{from}, \code{to} and \code{mass} that specifies from which knot to which other knot what amount of mass is sent in the optimal transference plan. Knots are given as indices in terms of the usual column major enumeration of the matrices \code{a$mass} and \code{b$mass}. There are \code{\link[=plot.pgrid]{plot}} methods for the classes \code{pgrid} and \code{pp}, which can plot this solution.

If \code{returncoarse} is \code{TRUE} for the \code{revsimplex} method, a list with components \code{sol} and \code{prob} giving the solutions and problems on the various scales considered. The solution on the finest scale (i.e. the output we obtain when setting \code{returncoarse} to \code{FALSE}) is in \code{sol[[1]]}.
}
\references{
F. Aurenhammer, F. Hoffmann and B. Aronov (1998). Minkowski-type theorems and least-squares clustering. Algorithmica 20(1), 61--76.

D. P. Bertsekas (1988). The auction algorithm: a distributed relaxation method for the assignment problem. Annals of Operations Research 14(1), 105--123.

D. P. Bertsekas (1992). Auction algorithms for network flow problems: a tutorial introduction. Computational Optimization and Applications 1, 7-66. 

C. Gottschlich and D. Schuhmacher (2014). The shortlist method for finding an optimal solution to the transportation problem. Preprint, available on \url{arxiv.org}.

D.G. Luenberger (2003). Linear and nonlinear programming, 2nd ed. Kluwer.

D.G. Luenberger and Y. Ye (2008). Linear and nonlinear programming, 3rd ed. Springer.

Q. Mérigot (2011). A multiscale approach to optimal transport. Computer Graphics Forum 30(5), 1583--1592.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de} \cr  
  Björn Bähre \email{bjobae@gmail.com} \cr (code for \code{\link{aha}}-method) \cr  
  Carsten Gottschlich \email{gottschlich@math.uni-goettingen.de} \cr (original java code for \code{shortlist}- and \code{revsimplex}-methods)
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=plot.pgrid]{plot}}, \code{\link{wasserstein}}
}
\examples{
#
# example for class 'pgrid'
#
dev.new(width=9, height=4.5)
par(mfrow=c(1,2), mai=rep(0.1,4))
image(random32a$mass, col = grey(0:200/200), axes=FALSE)
image(random32b$mass, col = grey(0:200/200), axes=FALSE)
res <- transport(random32a,random32b)
dev.new()
par(mai=rep(0,4))
plot(random32a,random32b,res,lwd=1)

#
# example for class 'pp'
#
set.seed(27)
x <- pp(matrix(runif(1000),500,2))
y <- pp(matrix(runif(1000),500,2))
res <- transport(x,y)
dev.new()
par(mai=rep(0.02,4))
plot(x,y,res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{main function}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
