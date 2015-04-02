\name{sglasso}
\alias{sglasso}
\alias{print.sglasso}
\title{Lasso Method for the RCON(V, E) Models}
\description{
Fit the weighted l1-penalized \emph{RCON(V, E)} models using a cyclic coordinate algorithm.
}
\usage{
sglasso(S, mask, w = NULL, min_rho = 1.0e-02, nrho = 50, nstep = 1.0e+05, 
        algorithm = c("ccd","ccm"), tol = 1.0e-05)
}
\arguments{
	\item{S}{the empirical variance/covariance matrix;}
	\item{mask}{a symmetric matrix used to define the structure of the concentration matrix. See the example bellow for more details;}
	\item{w}{a vector specifying the weights used to compute the weighted l1-norm of the parameters of the \emph{RCON(V, E)} model;}
	\item{min_rho}{last value of the sequence of tuning parameters used to compute the sglasso solution path. If the argument \code{nrho = 1}, then \code{min_rho} is the value used to compute the sglasso estimate. Default value is 1.0e-02;}
	\item{nrho}{number of tuning parameters. Default is 50;}
	\item{nstep}{nonnegative integer used to specify the maximun number of iterations of the two cyclic coordinate algorithms. Default is 1.0e+05;}
	\item{algorithm}{character by means of to specify the algorithm used to fit the model, i.e., a cyclic coordinate descente (\code{ccd}) algorithm or a cyclic coordinate minimization (\code{ccm}) algorithm. Default is \code{ccd};}
	\item{tol}{value used for convergence. Default value is 1.0e-05.}
}
\details{
	The \emph{RCON(V, E)} model (Hojsgaard et al., 2008) is a kind of restriction of the Gaussian Graphical Model defined using a coloured graph to specify a set of equality constraints on the entries of the concentration matrix. Roughly speaking, a coloured graph implies a partition of the vertex set into \eqn{R} disjoint subsets, called vertex colour classes, and a partition of the edge set into \eqn{S} disjoint subsets, called edge colour classes. At each vertex/edge colour class is associated a specific colour. If we denote by \eqn{K = (k_{ij})} the concentration matrix, i.e. the inverse of the variance/covariance matrix \eqn{\Sigma}, the coloured graph implies the following equality constraints:\enumerate{
	\item \eqn{k_{ii} = \eta_n} for any index \eqn{i} belonging to the \eqn{n}th vertex colour class;
	\item \eqn{k_{ij} = \theta_m} for any pair \eqn{(i,j)} belonging to the \eqn{m}th edge colour class.
	} Denoted with \eqn{\psi = (\eta',\theta')'} the \eqn{(R+S)}-dimensional parameter vector, the concentration matrix can be defined as \deqn{K(\psi) = \sum_{n=1}^R\eta_nD_n + \sum_{m=1}^S\theta_mT_m,}{K(\psi) = sum_{n=1}^R eta_n D_n + sum_{m=1}^S theta_m T_m,} where \eqn{D_n} is a diagonal matrix with entries \eqn{D^n_{ii} = 1} if the index \eqn{i} belongs to the \eqn{n}th vertex colour class and zero otherwise. In the same way, \eqn{T_m} is a symmetrix matrix with entries \eqn{T^m_{ij} = 1} if the pair \eqn{(i,j)} belongs to the \eqn{m}th edge colour class. Using the previous specification of the concentration matrix, the structured graphical lasso (sglasso) estimator (Abbruzzo et al., 2014) is defined as \deqn{\hat\psi = \arg\max_{\psi} \log det K(\psi) - tr\{Sk(\psi)\} - \rho\sum_{m=1}^Sw_m|\theta_m|,}{hat{\psi} = argmax_{\psi} log det K(\psi) - tr{Sk(\psi)} - \rho \sum_{m=1}^S w_m |\theta_m|,} where \eqn{S} is the empirical variance/covariance matrix, \eqn{\rho} is the tuning parameter used to control the ammount of shrinkage and \eqn{w_m} are weights used to define the weighted \eqn{\ell_1}-norm. By default, the \code{sglasso} function sets the weights equal to the cardinality of the edge colour classes.
}
\value{
\code{sglasso} returns an obejct with S3 class \code{"sglasso"}, i.e., a named list containing the following components:
	\item{call}{the call that produced this object;}
	\item{nv}{number of vertex colour classes;}
	\item{ne}{number of edge colour classes;}
	\item{theta}{the matrix of the sglasso estimates. The first \code{nv} rows correspond to the unpenalized parameters while the remaining rows correspond to the weighted l1-penalized parameters;}
	\item{w}{the vector of weights used to define the weighted l1-norm;}
	\item{df}{\code{nrho}-dimensional vector of the number of estimated nonzero parameters;}
	\item{rho}{\code{nrho}-dimensional vector of the sequence of tuning parameters;}
	\item{grd}{the matrix of the scores;}
	\item{nstep}{nonnegative integer used to specify the maximum number of iterations of the algorithms;}
	\item{nrho}{number of tuning parameters used to compute the sglasso solution path;}
	\item{algorithm}{the algorithm used to fit the model;}
	\item{tol}{a nonnegative value used to define the convergence of the algorithms;}
	\item{S}{the empirical variace/covariance matrix used to compute the sglasso solution path;}
	\item{mask}{the \code{mask} used to define the structure of the concentration matrix;}
	\item{n}{number of interations of the algorithm;}
	\item{conv}{an integer value used to encode the warnings related to the algorihtms. If \code{conv = 0} the convergence has been achieved otherwise the maximum number of iterations has been achieved.}
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\references{
Abbruzzo, A., Augugliaro, L., Mineo, A.M. and Wit, E.C. (2014) 
\emph{Cyclic coordinate for penalized Gaussian Graphical Models with symmetry restrictions}. In Proceeding of COMPSTAT 2014 - 21th International Conference on Computational Statistics, Geneva, August 19-24, 2014.
\cr\cr

Hojsgaard, S. and Lauritzen, S.L. (2008) 
\emph{Graphical gaussian models with edge and vertex symmetries}. J. Roy. Statist. Soc. Ser. B., Vol. 70(5), 1005--1027.
}
\seealso{
\code{\link{summary.sglasso}} and \code{\link{plot.sglasso}} methods.

The function \code{\link{Kh}} extracts the estimated sparse structured concentration matrices.
}
\examples{
########################################################
# sglasso solution path
#
## structural zeros:
## there are two ways to specify structural zeros which are 
## related to the kind of mask. If mask is a numeric matrix
## NA is used to identify the structural zero. If mask is a
## character matrix then the structural zeros are specified
## using NA or ".".
N <- 100
p <- 5
X <- matrix(rnorm(N * p), N, p)
S <- crossprod(X) / N
mask <- outer(1:p, 1:p, function(i,j) 0.5^abs(i-j))
mask[1,5] <- mask[1,4] <- mask[2,5] <- NA
mask[5,1] <- mask[4,1] <- mask[5,2] <- NA
mask

out.sglasso_path <- sglasso(S, mask, tol = 1.0e-13)
out.sglasso_path

rho <- out.sglasso_path$rho[20]
out.sglasso <- sglasso(S, mask, nrho = 1, min_rho = rho, tol = 1.0e-13, algorithm = "ccm")
out.sglasso

out.sglasso_path$theta[, 20]
out.sglasso$theta[, 1]
}
\keyword{models}
\keyword{graph}
\keyword{multivariate}
