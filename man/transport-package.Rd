\name{transport-package}
\alias{transport-package}
\docType{package}
\encoding{UTF-8}
\title{
Optimal Transport in Various Forms
}
\description{
Solve optimal transport problems in statistics. Compute Wasserstein metrics (a.k.a. Kantorovitch, Fortet--Mourier, Mallows, Earth Mover's, or minimal L_p metrics), return the corresponding transference plans, and display them graphically. Objects that can be compared include grey-scale images, point patterns, and mass vectors.
}
\details{
\tabular{ll}{
Package: \tab transport\cr
Type: \tab Package\cr
Version: \tab 0.6-1\cr
Date: \tab 2014-05-30\cr
License: \tab GPL (>=2) \cr
LazyData: \tab yes\cr
}

The main end-user function is \code{\link{transport}}. It computes optimal transference plans between images (class \code{\link{pgrid}}), point patterns (class \code{\link{pp}}), or mass vectors, based on various algorithms. These transference plans can \code{\link[=plot.pgrid]{plot}}ed. The function \code{\link{wasserstein}} allows for the numerical computation
of \eqn{p}-th order Wasserstein distances.
}
\author{
Dominic Schuhmacher <dominic.schuhmacher@mathematik.uni-goettingen.de> with substantial contributions of code by Björn Bähre and Carsten Gottschlich

Maintainer: Dominic Schuhmacher \email{dominic.schuhmacher@mathematik.uni-goettingen.de}
}
\references{
  See help page for the function \code{\link{transport}}.
}
\keyword{package}
\examples{
##  See examples for function transport
}
