\name{sglasso-package}
\alias{sglasso-package}
\docType{package}
\title{
Lasso Method for RCON(V, E) Models
}
\description{
\emph{RCON(V, E)} models (Hojsgaard, et al., 2008) are a kind of restriction of the Gaussian Graphical Models defined by a set of equality constraints on the entries of the concentration matrix. \code{sglasso} package implements the structured graphical lasso (sglasso) estimator proposed in Abbruzzo et al. (2014) for the weighted l1-penalized \emph{RCON(V, E)} model. Two cyclic coordinate algorithms are implemented to compute the sglasso estimator, i.e. a cyclic coordinate minimization (CCM) and a cyclic coordinate descent (CCD) algorithm.
}
\details{
\tabular{ll}{
Package: \tab sglasso\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-05-04\cr
License: \tab GPL (>=2)
}
}
\author{
Luigi Augugliaro\cr\cr
Maintainer:\cr
Luigi Augugliaro <luigi.augugliaro@unipa.it>
}
\references{
Abbruzzo, A., Augugliaro, L., Mineo, A. M. and Wit, E. C. (2014)
Cyclic coordinate for penalized Gaussian Graphical Models with symmetry restrictions. In \emph{Proceeding of COMPSTAT 2014 - 21th International Conference on Computational Statistics}, Geneva, August 19-24, 2014.
\cr

Hojsgaard, S. and Lauritzen, S. L. (2008)
Graphical gaussian models with edge and vertex symmetries. \emph{J. Roy. Statist. Soc. Ser. B.}, Vol. \bold{70(5)}, 1005--1027.
}
\keyword{multivariate}
\keyword{graph}
\keyword{package}
