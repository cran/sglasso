\name{plot.sglasso}
\alias{plot.sglasso}
\title{Plot Method for the Weighted l1-Penalized RCON(V, E) Model}
\description{
\code{plot.sglasso} produces two plots to study the sequence of models estimates by \code{sglasso}.
}
\usage{
\method{plot}{sglasso}(x, \dots)
}
\arguments{
	\item{x}{fitted \code{sglasso} object;}
	\item{\dots}{other parameters to be passed through the plotting functions.}
}
\details{
This function produces two different plots. The first one shows the path of the estimated parameters as function of the tuning parameter \eqn{\rho}{rho}. In the same way, the second plot shows the path of the weighted scores as function of \eqn{\rho}{rho}.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{sglasso}} function and \code{\link{summary.sglasso}} method.
}
\examples{
N <- 100
p <- 5
X <- matrix(rnorm(N * p), N, p)
S <- crossprod(X)/N
mask <- outer(1:p, 1:p, function(i,j) 0.5^abs(i-j))
mask[1,5] <- mask[1,4] <- mask[2,5] <- NA
mask[5,1] <- mask[4,1] <- mask[5,2] <- NA
out.sglasso_path <- sglasso(S, mask, tol = 1.0e-13)
plot(out.sglasso_path)
}
\keyword{models}
