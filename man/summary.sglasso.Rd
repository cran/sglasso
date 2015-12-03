\name{summary.sglasso}
\alias{summary.sglasso}
\title{Summarizing sglasso Fits}
\description{
\code{summary} method for class \code{"sglasso"}.
}
\usage{
\method{summary}{sglasso}(object, N, k = c("bic","aic"), 
        digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{object}{fitted \code{sglasso} object;}
  \item{N}{sample size;}
  \item{k}{character/numeric argument used to specify the 'weight' of the complexity part in the measure of goodness-of-fit used to select the best model (see below for more details). Default is \code{k = "bic"};}
  \item{digits}{significant digits in printout;}
  \item{\dots}{additional print arguments.}
}
\details{
	\code{summary.sglasso} gives us information about the sequence of models estimated by the sglasso estimator. To select the best model, summary method uses a measure of Goodness-of-Fit (GoF) defined as follows: \deqn{-2\ell(\hat\psi) + k \times df,}{-2 l(hat{\psi}) + k df,} where \eqn{\ell(\hat\psi)}{l(hat{\psi})} is the log-likelihood of the estimated weighted l1-penalized \emph{RCON(V, E)} model, \eqn{df} is the number of nonzero estimated parameters and \eqn{k} is a non-negative value used to weight the complexity part in the measure of goodness-of-fit. By default the summary method computes the BIC criterion to select the best model (\code{k = "bic"}). The AIC criterion can be easily computed setting \code{k = "aic"}. The user can also define other measures of goodness-of-fit specifying \code{k} as any non-negative value.
	
	The output of the summary method is divided in two sections. First section shows the call producing the argument \code{object} followed by a \code{data.frame}. The column named \code{rho} shows the sequence of the \eqn{\rho} values used to compute the solution curve, while the column \code{log-lik} shows the corresponding values of the log-likelihood function. The remaining columns show the number of estimated non-zero parameters, the values of the GoF and the asscoated ranking of the estimated models. Finally, the second section shows the estimated parameters of the best model identified by the used GoF criterion. Informations about the algorithm and the corresponding convergence are also provided.
}
\value{
	A list with components \code{table} and \code{theta_gof} is silently returned. The \code{table} component is the \code{data.frame} previously described while the component \code{theta_gof} is the vector of the estimated parameters corresponding to the best models identified by the GoF criterion.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{sglasso}} and \code{\link{loglik}} functions.
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
summary(out.sglasso_path, N)
rho <- out.sglasso_path$rho[20]
out.sglasso <- sglasso(S, mask, nrho = 1, min_rho = rho, tol = 1.0e-13)
summary(out.sglasso, N)
}
\keyword{models}
