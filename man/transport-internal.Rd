\name{transport-internal}
\alias{dedegenerate}
\alias{findblocks}
\alias{fudge}
\alias{image2}
\alias{image3}
\alias{refinesol}
\alias{triangulate}

\title{
  Internal Functions in Package \sQuote{transport}
}
\description{
  Internal functions in package \sQuote{transport}
}
\usage{
dedegenerate(basis)
% Fix a degenerate starting solution for function transport. I.e. if the (m x n)-basis matrix has fewer
% than the required m+n-1 ones add the missing number of ones in such a way that no cycle (loop) is 
% generated 
findblocks(tbasis)
% finds the connected blocks in a basis matrix and returns them as list of lists of row and col indices
% (loops indefinitely if not every row and every col contains at least one 1)
fudge(temp, N = 1e+09)
% transforms vectors with natural entries to vectors with natural entries of the fixed total mass N 
% by randomly adding removing units of mass at entries that have positive mass
image2(x, y, z, rot = FALSE, ...)
image3(z, x=1:dim(z)[1], y=1:dim(z)[2], rot=TRUE, ...)
% like image, but if rot=TRUE, matrices are plotted in the same orientation as they are written
refinesol(a1, b1, a2, b2, assig1, basis1, mult = 2)
% for the multiscale approach. Computes a starting solution on the finer grid (finer by factor mult)
% based on the solution given by assig1 and basis1 on the coarser grid
triangulate(basis)
% compute triangulation of m x n basis matrix
}
%- maybe also 'usage' for other objects documented here.
\details{
  These internal functions should usually not be called directly by the user.
  They perform more elementary tasks, do little to no checks of their inputs,
  and may change from one package version to the next without notice.
  Also, most of these functions are not exported, so call them by
  \code{transport:::name_of_fun} if you must.
}
\author{
  Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de}
}
\keyword{internal}

