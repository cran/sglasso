\name{gplot}
\alias{gplot}
\title{Plotting Sparse Graph}
\description{
\code{gplot} is a generic function for plotting sparse graphs.
}
\usage{
gplot(object, \dots)
}
\arguments{
	\item{object}{fitted \code{sglasso}/\code{fglasso} object;}
	\item{\dots}{other parameters passed to \code{gplot.sglasso} or \code{gplot.fglasso}.}
}
\details{
\code{gplot} is a generic function used to plot a graph estimated by \code{sglasso} or \code{fglasso}. See the method function \code{gplot.sglasso} or \code{gplot.fglasso} for more details about the specific arguments.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\seealso{
\code{\link{gplot.sglasso}} and \code{\link{gplot.fglasso}} method functions.
}
\keyword{models}
