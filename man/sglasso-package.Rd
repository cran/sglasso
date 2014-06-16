\name{sglasso-package}
\alias{sglasso-package}
\docType{package}
\title{
Lasso Method for RCON(V, E) Models
}
\description{
	\emph{RCON(V, E)} models (Hojsgaard, et al., 2008) are a kind of restriction of the Gaussian Graphical Models defined by a set of equality constraints on the entries of the concentration matrix. \code{sglasso} package implements the structured graphical lasso (sglasso) estimator proposed in Augugliaro et al. (2014) for the weighted l1-penalized \emph{RCON(V, E)} model. Two cyclic coordinate algorithms are implemented to compute the sglasso estimator, i.e. a cyclic coordinate minimization (CCM) algorithm and a cyclic coordinate descent (CCD) algorithm.
}
\details{
\tabular{ll}{
Package: \tab sglasso\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2014-06-16\cr
License: \tab GPL (>=2)
}
}
\author{
Luigi Augugliaro\cr\cr
Maintainer:\cr
Luigi Augugliaro <luigi.augugliaro@unipa.it>
}
\references{
Augugliaro, L., Abbruzzo, A., Mineo, A.M. and Wit, E.C. (2014) 
\emph{The structured graphical lasso estimator}. Technical Report.
\cr\cr
Hojsgaard, S. and Lauritzen, S.L. (2008) 
\emph{Graphical gaussian models with edge and vertex symmetries}. J. Roy. Statist. Soc. Ser. B., Vol. 70(5), 1005--1027.
}
\keyword{multivariate}
\keyword{graph}
\keyword{package}
