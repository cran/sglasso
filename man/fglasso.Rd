\name{fglasso}
\alias{fglasso}
\title{Factorial Graphical Lasso Model}
\description{
Fit the factorial graphical lasso model.
}
\usage{
fglasso(S, model, tp, p, ...)
}
\arguments{
	\item{S}{the empirical variance/covariance matrix;}
	\item{model}{a list or a matrix used to specify the factorial Gaussian Graphical Model (see Details);}
	\item{tp}{number of time points;}
	\item{p}{number of random variables observed for each time point;}
	\item{...}{further arguments passed to \code{sglasso}.}
}
\details{
The factorial graphical model (Abbruzzo et al., to appear) is a special kind of the Graphical Gaussian Model proposed to study dynamic networks. Let \eqn{X_t = (X_{it},\ldots,X_{it})'}{X_t = (X_{it},...,X_{it})'} be a \eqn{p}{p}-dimensional random variable at time \eqn{t}{t}. Assuming that \eqn{X = (X'_1,\ldots,X'_T)}{X = (X'_1,...,X'_T)'} follows a multivariate normal distribution, the concentration matrix \eqn{K}{K} has the following block structure
\deqn{
K = \left(
\begin{array}{cccc}
K_{1,1} & K_{1,2} & \ldots & K_{1,T}\\
K_{2,1} & K_{2,2} & \ldots & K_{2,T}\\
\vdots & \vdots & \ddots & \vdots\\
K_{T,1} & K_{T,2} & \ldots & K_{T,T}
\end{array}\right),
}
where \eqn{K_{t,t}} give information about the conditinal independence structure among the \eqn{p} random variables at time \eqn{t}, and \eqn{K_{t,t + h}} give information about the conditional independence structure between \eqn{X_t} and \eqn{X_{t + h}}. An interpretation of the elements of the submatrices \eqn{K_{t,t + h}} brings to the notion of natural structure, i.e.,
\deqn{
K_{t,t + h} = \left(
\begin{array}{cccc}
k_{1t,1(t+h)} & 0 & \ldots & 0\\
0 & k_{2t,2(t+h)} & \ldots & 0\\
\vdots & \vdots & \ddots & \vdots\\
0 & 0 & \ldots & k_{pt,p(t+h)}
\end{array}\right) + \left(
\begin{array}{cccc}
0 & k_{1t,2(t+h)} & \ldots & k_{1t,p(t+h)}\\
k_{2t,1(t+h)} & 0 & \ldots & 0\\
\vdots & \vdots & \ddots & \vdots\\
k_{pt,1(t+h)} & k_{pt,2(t+h)} & \ldots & 0
\end{array}\right).
}
The entries of the first matrix are called \bold{self-self conditinal dependences} at temporal lag \eqn{h} and represents the (negative) self-similarity of a given random variable across different time points. The entries of the second matrix are the \bold{conditional dependence} among the \eqn{p} random variables. To make the interpretation of the results more relevant and, at the same time, reduce the number of parameters, the authors propose the following equality constraints:
\tabular{clccclcc}{
\tab self-self dependence \tab effect \tab \code{R} code \tab \tab conditinal dependence \tab effect \tab \code{R} code \cr
i.	\tab \eqn{k_{it,i(t+h)} = 0} \tab zero \tab \code{"."} \tab iv. \tab \eqn{k_{it,j(t+h)} = 0} \tab zero \tab \code{"."} \cr
ii.	\tab \eqn{k_{it,i(t+h)} = s^h} \tab costant \tab \code{"c"} \tab ii. \tab \eqn{k_{it,j(t+h)} = n^h} \tab costant \tab \code{"c"} \cr
iii.	\tab \eqn{k_{it,i(t+h)} = s^h_i} \tab unit \tab \code{"u"} \tab iii. \tab \eqn{k_{it,j(t+h)} = n^h_i} \tab unit \tab \code{"u"} \cr
iv.	\tab \eqn{k_{it,i(t+h)} = s^{t,h}} \tab time \tab \code{"t"} \tab iv. \tab \eqn{k_{it,j(t+h)} = n^{t,h}} \tab time \tab \code{"t"} \cr
v.	\tab \eqn{k_{it,i(t+h)} = s^{t,h}_{i,t}} \tab interaction \tab \code{"ut"} \tab v. \tab \eqn{k_{it,j(t+h)} = n^{t,h}_{i,t}} \tab interaction \tab \code{"ut"}
}
Argument \code{model} is used to specify the restrinctions previously describted. This argument can be a named list or a matrix with dimension \eqn{nlag\times 2}, where \eqn{nlag\le\code{tp}}. To gain more insight, suppose that we want to model only the sub-matrices \eqn{K_{t,t}} and \eqn{K_{t,t+1}}, i.e., the sub-matrices corresponding to the temporal lag zero and one. A possible \code{R code} is\cr

\code{> model.mat <- matrix("", nrow = 2, ncol = 2)}\cr
\code{> rownames(model.mat) <- c("lag0", "lag1")}\cr
\code{> colnames(model.mat) <- c("s", "n")}\cr
\code{> model.mat[1, ] <- c("c", "ut")}\cr
\code{> model.mat[2, ] <- c("t", ".")}\cr

In this example we are modelling the diagonal elements of the sub-matrices \eqn{K_{t,t}} with the constant effect while the off-diagonal elements are modelled by the interaction effect. In the same way, the diagonal elements of the sub-matrices \eqn{K_{t,t+1}} are modelled by the time effect while the remaning element are equal to zero. The \code{fglasso} function passes the matrix \code{model.mat} to the internal function \code{fglasso_model2mask}, i.e.,\cr 

\code{> mask <- fglasso_model2mask(model.mat, tp = 3, p = 3)}\cr
\code{> mask}\cr

which returns the mask used in \code{sglasso} to fit the specified factorial Gaussian Graphical Model. The same model can be specified by the following named list\cr

\code{> model.list <- list(lag0 = c(s = "c", n = "ut"), lag1 = c(s = "t", n = "."))}\cr

See the example below for more details.
}
\value{
\code{fglasso} returns an obejct with S3 class \code{"sglasso"}. See the corresponding manual for more details.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}}
\references{
Abbruzzo, A. and Wit, E.C. (2014) 
\emph{Dynamic factorial graphical models for dynamic networks}. To appear in Network Science.
\cr\cr

Abbruzzo, A., Augugliaro, L., Mineo, A.M. and Wit, E.C. (2014) 
\emph{Cyclic coordinate for penalized Gaussian Graphical Models with symmetry restrictions}. In Proceeding of COMPSTAT 2014 - 21th International Conference on Computational Statistics. Geneva, August 19-24, 2014.
}
\seealso{
\code{\link{sglasso}} function.
}
\examples{
#######################
# fglasso solution path
#######################

N <- 50
tp <- 3
p <- 3
X <- matrix(rnorm(N * p * tp), N, tp * p)
S <- crossprod(X) / N
model <- list(lag0 = c(s = "c", n = "ut"), lag1 = c(s = "t", n = "."))
out.fglasso <- fglasso(S = S, model = model, tp = tp, p = p)
out.fglasso
}
\keyword{models}
\keyword{graph}
\keyword{multivariate}
