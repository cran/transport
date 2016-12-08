\name{semidiscrete}
\alias{semidiscrete}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find Optimal Transport Partition Between pgrid and wpp.
}
\description{
  Given an object \code{a} of class \code{\link[=pgrid-object]{pgrid}} specifying an image and an object \code{b}
  of class \code{\link[=wpp-object]{wpp}} specifiying a more flexible mass distribution at finitely many points,
  find the partition of the image (and hence the optimal transport map) that minimizes the total transport cost
  for going from \code{a} to \code{b}.
}
\usage{
  semidiscrete(a, b, p = 2, method = c("aha"), control = list(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
  	an object of class \code{\link[=pgrid-object]{pgrid}} usually representing an image or the discretization of a measure.
}
  \item{b}{
    an object of class \code{\link[=wpp-object]{wpp}} usually having the same total mass as \code{a}.   	
}
  \item{p}{
    the power \eqn{\geq 1}{>=1} to which the Euclidean distance between points is taken in order
    to compute costs.
}
  \item{method}{
    the name of the algorithm to use. Currently only \code{aha} is supported.
}
  \item{control}{
  	a named list of parameters for the chosen method or the result of a call to \code{\link{trcontrol}}. Currently only
  	the parameters \code{factr} and \code{maxit} can be set.
  }
  \item{...}{
  	currently without effect.
  }
}
\details{
  This is a wrapper function for the function \code{\link{aha}}. In the latter the
  Aurenhammer--Hoffmann--Aronov (1998) method is used with the multiscale approach presented in Mérigot (2011). 
  
  This function is not usually called by the user since it is automatically called by \code{\link{transport}} if the
  first argument is of class \code{pgrid} and the second argument is of class \code{wpp}.
}
\value{
  An object of class \code{power_diagram} describing the optimal transport partition for \code{a} and \code{b}. Such an
  object has components \code{sites} and \code{cells}.
  
  The former is a data.frame with columns \code{xi}, \code{eta} and \code{w} containing
  the \eqn{x}- and \eqn{y}-coordinates of the mass locations of \code{b} and the weight vector \eqn{w} that describes the
  associated power tesselation uniquely. The latter is a list with as many components as mass locations of \code{b}, each
  giving a matrix containing in the columns the \eqn{x}- and \eqn{y}-coordinates of the polygonal cell associated with
  that location or \code{NULL} if the location has no cell.
  
  Plotting methods exist both for objects of class \code{\link{power_diagram}} and for \link[=plot.pgrid]{optimal transport maps
  represented by power diagrams}.
}
\references{
F. Aurenhammer, F. Hoffmann and B. Aronov (1998). Minkowski-type theorems and least-squares clustering. Algorithmica 20(1), 61--76.

Q. Mérigot (2011). A multiscale approach to optimal transport. Computer Graphics Forum 30(5), 1583--1592. \href{http://dx.doi.org/10.1111/j.1467-8659.2011.02032.x}{doi:10.1111/j.1467-8659.2011.02032.x}
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de} \cr  
  Björn Bähre \email{bjobae@gmail.com} \cr  
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=plot.pgrid]{plot}}, \code{\link{transport}}, \code{\link{aha}}
}
\examples{
##  See examples for function transport}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{semi-dicrete}
\keyword{transport partition}
\keyword{transport tesselation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
