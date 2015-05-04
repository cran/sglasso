\name{plot.klcv}
\alias{plot.klcv}
\title{Plot Method for Leave-One-Out Cross-Validated Kullback-Leibler Divergence}
\description{
\code{plot.klcv} produces a plot to study the sequence of leave-one-out cross-validated Kullback-Leibler divergences computed by \code{klcv}.
}
\usage{
\method{plot}{klcv}(x, \dots)
}
\arguments{
	\item{x}{fitted \code{klcv} object;}
	\item{\dots}{other parameters to be passed through the plotting function.}
}
\details{
This method function produces a plot showing the sequence of leave-one-out cross-validated Kullback-Leibler as function of the tuning parameter \eqn{rho}. The optimal value of the tuning parameter is identified by a vertical dashed line.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{klcv}} function.
}
\examples{
N <- 100
p <- 5
X <- matrix(rnorm(N * p), N, p)
S <- crossprod(X) / N
mask <- outer(1:p, 1:p, function(i,j) 0.5^abs(i-j))
mask[1,5] <- mask[1,4] <- mask[2,5] <- NA
mask[5,1] <- mask[4,1] <- mask[5,2] <- NA
out.sglasso_path <- sglasso(S, mask, tol = 1.0e-13)
out.klcv <- klcv(out.sglasso_path, X)
plot(out.klcv)
}
\keyword{models}
