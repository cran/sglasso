\name{Kh}
\alias{Kh}
\title{Extract Sparse Structured Precision Matrix}
\description{
Function \code{Kh} computes the sequence of sparse structured precision matrices estimated by the \code{sglasso} function.
}
\usage{
Kh(object, rho)
}
\arguments{
	\item{object}{fitted \code{sglasso} object;}
	\item{rho}{a vector specifying a subset of the values of the tuning parameter used in \code{sglasso} to compute the solution path. By default, the entire sequence of estimated sparse structured precision matrices is returned.}
}
\value{
\code{Kh} returns a named list containing the sequence of estimated sparse structured precision matrices.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{\code{\link{sglasso}} function.}
\examples{
N <- 100
p <- 5
X <- matrix(rnorm(N * p), N, p)
S <- crossprod(X)/N
mask <- outer(1:p, 1:p, function(i, j) 0.5^abs(i-j))
out.sglasso_path <- sglasso(S, mask, nrho = 5, tol = 1.0e-13)
out.sglasso_path
Kh(out.sglasso_path)
rho <- out.sglasso_path$rho[3]
out.sglasso_single <- sglasso(S, mask, nrho = 1, min_rho = rho, tol = 1.0e-13, algorithm = "ccm")
Kh(out.sglasso_single)
}
\keyword{models}
\keyword{graph}
