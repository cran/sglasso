\name{loglik}
\alias{loglik}
\title{Extract Log-Likelihood}
\description{
This function extracts the log-likelihood for the sequence of weighted l1-penalized \emph{RCON(V, E)} models estimated by the \code{sglasso} function.
}
\usage{
loglik(object, size = 2)
}
\arguments{
	\item{object}{a fitted \code{sglasso} object;}
	\item{size}{sample size. Default value is 2 to remove the constant term in the log-likelihood function. See below for more details.}
}
\details{
Denoted with \eqn{\psi = (\eta^T,\theta^T)^T} the parameter vector of the structured concentration matrix \eqn{K(\psi)}, the log-likelihood function of the \emph{RCON(V, E)} model is equal, up to a constant, to the following expression \deqn{\ell(\psi) = \frac{N}{2}[\log det K(\psi) - tr\{S K(\psi)\}],}{l(\psi) = N/2[log det K(\psi) - tr{S K(\psi)}],} where \eqn{S = N^{-1}\sum_{i=1}^NX_iX_i^T}{S = N^{-1}\sum_{i=1}^N X_i X_i^T}, \eqn{N} is the sample size and \eqn{X_i} is the \eqn{i}th observed \eqn{p}-dimensional vector. Denoted with \eqn{\hat\psi = (\hat\eta^T,\hat\theta^T)^T}{hat{\psi} = (hat{\eta}^T, hat{\theta}^T)^T} the sglasso estimates, straightforward algebra shows that \deqn{\ell(\hat\psi) = \frac{N}{2}[\log det K(\hat\psi) - p + \rho\sum_{m=1}^S w_m|\hat\theta_m|],}{l(hat{\psi}) = N/2[log det K(hat{\psi}) - p + \rho\sum_{m=1}^S w_m |hat{\theta_m|}],} where \eqn{\rho} is the tuning parameter and \eqn{w_m} are the weights used to define the weighted l1-norm.
}
\value{
\code{loglik} returns a vector containing the log-likelihood computed for the sequence of weighted l1-penalized \emph{RCON(V, E)}.
}
\author{Luigi Augugliaro\cr
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{summary.sglasso}} method and \code{\link{sglasso}} function.
}
\examples{
N <- 100
p <- 5
X <- matrix(rnorm(N * p), N, p)
S <- crossprod(X)/N
mask <- outer(1:p, 1:p, function(i, j) 0.5^abs(i-j))
out.sglasso_path <- sglasso(S, mask, nrho = 5, tol = 1.0e-13)
out.sglasso_path
loglik(out.sglasso_path, size = N)
rho <- out.sglasso_path$rho[3]
out.sglasso_single <- sglasso(S, mask, nrho = 1, min_rho = rho, tol = 1.0e-13, algorithm = "ccm")
loglik(out.sglasso_single, size = N)
}
\keyword{models}
\keyword{graph}
