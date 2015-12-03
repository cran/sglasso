\name{gplot.sglasso}
\alias{gplot.sglasso}
\title{Plotting Sparse Graphs}
\description{
\code{gplot.sglasso} shows the sequence of graphs estimated by \code{sglasso}.
}
\usage{
\method{gplot}{sglasso}(object, rhoid, layout = layout.circle, \dots)
}
\arguments{
	\item{object}{fitted \code{sglasso} object;}
	\item{rhoid}{vector of integers used to specificy the \eqn{\rho}{rho}-values used to fit the sglasso model. By default \code{gplot.sglasso} shows the sequence of graphs estimated by \code{sglasso}. Only topologically different graphs are plotted;}
	\item{layout}{a function or a matrix used to specify the layout of the graphs that will be plotted. By default the \code{layout.circle} function is used;}
	\item{\dots}{further graphical parameters used to plot the graphs. See package \pkg{igraph} for more details.}
}
\details{
\code{gplot.sglasso} shows the sequence of topologically different graphs estimated by \code{sglasso}. To specify the layout of the graphs, the user can use any layout function available in the R package \pkg{igraph}. The user can also specify the placement of the vertices by a matrix with two columns and the same number of rows as the number of vertices.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}
\seealso{
\code{\link{sglasso}} function.
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
gplot(out.sglasso_path)
gplot(out.sglasso_path, rhoid = 1:5)
}
\keyword{models}
