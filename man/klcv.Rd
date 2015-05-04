\name{klcv}
\alias{klcv}
\alias{print.klcv}
\title{Cross-Validated Kullback-Leibler Divergence}
\description{
Model selection criterion based on the leave-one-out cross-validated Kullback-Leibler divergence.
}
\usage{
klcv(object, X, scale = 1)
}
\arguments{
  \item{object}{fitted \code{sglasso}/\code{fglasso} object;}
  \item{X}{the matrix used to compute the empirical variance/covariance matrix. Its dimension is \code{N} \eqn{\times}{x} \code{p}, where \code{p} is the number of random variables and \code{N} is the samlpe size;}
  \item{scale}{scalar value used to scale the estimated degrees-of-freedom. See below for more details.}
}
\details{
\code{klcv} function implements the leave-one-out cross-validate Kullback-Leibler divergence criterion proposed in Vujacic et al. (2015). For \eqn{l_1}-penalized Gaussian Graphical Models this measure of goodness-of-fit has the following form \deqn{klcv(\rho) = -\frac{\ell(\hat K(\rho))}{N} + \frac{\code{scale}}{2N} gdf(\hat K(\rho)),} where \eqn{\hat K(\rho)} is the glasso estimate of the concentration matrix, \eqn{\ell(\hat K(\rho))} is the corresponding value of the log-likelihood function, \code{scale} is a scale factor for the complexity part, i.e. \eqn{gdf(\hat K(\rho))}, which is defined as \deqn{gdf(\hat K(\rho)) = \frac{1}{N-1}\sum_{k=1}^N vec\{(\hat K(\rho)^{-1} - S_k)\circ 1_\rho\}'vec[\hat K(\rho)\{(S-S_k)\circ 1_\rho\}\hat K(\rho)].} In the previous expression \eqn{S} is the empirical variance/covariance matrix, \eqn{S_k = X_k X_k'}, \eqn{1_\rho} is a matrix with entries \eqn{I(\hat k_{ij}(\rho)\ne 0)} and \eqn{\circ} is the Hadamard product operator.
}
\value{
	\code{klcv} returns an S3 object with calls \code{klcv}, i.e. a named list with the following components:
		\item{klcv}{the vector with the leave-one-out cross-validated Kullback-Leibler divergence;}
		\item{rho}{the rho-values used to compute the leave-one-out cross-validated Kullback-Leibler divergence;}
		\item{loglik}{a vector with the log-likelihood computed for the sequence of weighted l1-penalized \emph{RCON(V, E);}}
		\item{gdf}{a vector returning the generalized degrees-of-freedom;}
		\item{scale}{the scale value used to define the leave-one-out cross-validated Kullback-Leibler divergence;}
		\item{min.klcv}{minimum value of the leave-one-out cross-validated Kullback-Leibler divergence;}
		\item{rho.opt}{the rho-value corresponding to minimum leave-one-out cross-validated Kullback-Leibler divergence;}
		\item{rhoid}{the index of the rho-value identified by the leave-one-out cross-validated Kullback-Leibler divergence.}
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\references{
Vujacic, I., Abbruzzo, A. and Wit, E. C. (2015) A computationally fast alternative to cross-validation in penalized Gaussian graphical models. \emph{J. Stat. Comput. Simul.}
}
\seealso{
\code{\link{sglasso}}, \code{\link{loglik}} functions and \code{\link{plot.klcv}} method.
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
out.klcv
}
\keyword{models}
