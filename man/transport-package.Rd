\name{transport-package}
\alias{transport-package}
\docType{package}
\encoding{UTF-8}
\title{
Optimal Transport in Various Forms
}
\description{
Solve optimal transport problems. Compute Wasserstein distances (a.k.a. Kantorovitch, Fortet--Mourier, Mallows, Earth Mover's, or minimal \eqn{L_p} distances), return the corresponding transport plans, and display them graphically. Objects that can be compared include grey-scale images, (weighted) point patterns, and mass vectors.
}
\details{
\tabular{ll}{
Package: \tab transport\cr
Type: \tab Package\cr
Version: \tab 0.9-4\cr
Date: \tab 2017-10-03\cr
License: \tab GPL (>=2) \cr
LazyData: \tab yes\cr
}

The main end-user function is \code{\link{transport}}. It computes optimal transport plans between images (class \code{\link{pgrid}}), point patterns (class \code{\link{pp}}), weighted point patterns (class \code{\link{wpp}}) and mass vectors, based on various algorithms. These transport plans can be \code{\link[=plot.pgrid]{plot}}ed. The function \code{\link{wasserstein}} allows for the numerical computation of \eqn{p}-th order Wasserstein distances.

Most functions in this package are designed for data in two and higher dimensions. A quick tool for computing the \eqn{p}-th order Wasserstein distance between univariate samples is \code{\link{wasserstein1d}}.

}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de} \cr  
  Björn Bähre \email{bjobae@gmail.com} (code for \code{\link{aha}}-method) \cr  
  Carsten Gottschlich \email{gottschlich@math.uni-goettingen.de} \cr __(original java code for \code{shortlist}- and \code{revsimplex}-methods) \cr
  Florian Heinemann \email{florian.heinemann@stud.uni-goettingen.de} (\code{\link{transport_track}}-function) \cr
  Bernhard Schmitzer \email{schmitzer@uni-muenster.de} (\code{\link{shielding}}-method)

Maintainer: Dominic Schuhmacher \email{dominic.schuhmacher@mathematik.uni-goettingen.de}
}
\references{
  See help page for the function \code{\link{transport}}.
}
\keyword{package}
\examples{
##  See examples for function transport
}
