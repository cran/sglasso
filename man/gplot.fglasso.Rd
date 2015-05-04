\name{gplot.fglasso}
\alias{gplot.fglasso}
\title{Plotting Sparse Factorial Dynamic Gaussian Graphical Model}
\description{
\code{gplot.fglasso} shows the sequence of graphs estimated by \code{fglasso}.
}
\usage{
\method{gplot}{fglasso}(object, rhoid, tp = c(1, 2), sub.tp1, sub.tp2, cex.sub = 1, 
    k = 1.5, layout = layout.circle, \dots)
}
\arguments{
	\item{object}{fitted \code{fglasso} object;}
	\item{rhoid}{an integer used to specificy the \eqn{\rho}{rho}-value used to fit the fglasso model;}
	\item{tp}{a vector of length equal to two used to specify the time points of the two graphs that will be compared. By default the first two time points are used;}
	\item{sub.tp1}{sub title for the graph estimated at time point \code{tp[1]};}
	\item{sub.tp2}{sub title for the graph estimated at time point \code{tp[2]};}
	\item{cex.sub}{a numerical value giving the amount by which plotting sub titles should be magnified relateve to the default;}
	\item{k}{value used to specify the distance between the two graphs;}
	\item{layout}{a function or a matrix used to specify the layout of the graphs that will be plotted. By default the \code{layout.circle} function is used;}
	\item{\dots}{further graphical parameters used to plot the graphs. See package \pkg{igraph} for more details.}
}
\details{
For a given value of the tuning parameter, specified by the argument \code{rhoid}, \code{gplot.fglasso} shows the graphs estimated at the time points \code{tp[1]} and \code{tp[2]}. By convention, the graph associated to the sub matrix \eqn{K_{tp[1],tp[2]}} is represented by a directed graph where a directed edge is drawn by an arrow from a vertex in the first graph pointing forwards a vertex in the second graph.
}
\value{
\code{gplot.fglasso} returns a list with components:
	\item{graph.tp1}{an object with class \code{igraph} representing the undirected graph estimated at the time point \code{tp[1]};}
	\item{graph.tp2}{an object with class \code{igraph} representing the undirected graph estimated at the time point \code{tp[2]};}
	\item{graph.net}{an object with class \code{igraph} representing the directed graph associated to the submatrix \eqn{K_{tp[1],tp[2]}};}
	\item{layout}{the matrix used to specify the placement of the vertices.}
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}
\seealso{
\code{\link{fglasso}} function.
}
\examples{
N <- 50
tp <- 3
p <- 3
X <- matrix(rnorm(N * p * tp), N, tp * p)
S <- crossprod(X) / N
model <- list(lag0 = c(s = "c", n = "ut"), lag1 = c(s = "t", n = "t"))
out.fglasso <- fglasso(S = S, model = model, tp = tp, p = p)
gplot(out.fglasso, rhoid = 50, sub.tp1 = "First graph", 
   sub.tp2 = "Second graph")
}
\keyword{models}
