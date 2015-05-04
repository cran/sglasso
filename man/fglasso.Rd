\name{fglasso}
\alias{fglasso}
\title{L1-penalized Factorial Graphical Lasso Model}
\description{
Fit the weight l1-penlized factorial dynamic Gaussian Graphical Model.
}
\usage{
fglasso(S, model, tp, p, \dots)
}
\arguments{
	\item{S}{the empirical variance/covariance matrix;}
	\item{model}{a list or a matrix used to specify the factorial dynamic Gaussian Graphical Model (see Details);}
	\item{tp}{number of time points;}
	\item{p}{number of random variables observed for each time point;}
	\item{\dots}{further arguments passed to \code{sglasso}.}
}
\details{
The factorial dynamic Gaussian Graphical Model (Abbruzzo et al., 2015) is a special kind of \emph{RCON(V, E)} model (Hojsgaard, et al.,2008) proposed to study dynamic networks. Let \eqn{X_t = (X_{it},\ldots,X_{it})'}{X_t = (X_{it},\dots,X_{it})'} be a \eqn{p}{p}-dimensional random variable at time \eqn{t}{t}. Assuming that \eqn{X = (X'_1,\ldots,X'_T)}{X = (X'_1,...,X'_T)'} follows a multivariate normal distribution, the concentration matrix \eqn{K}{K} has the following block structure
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
k_{1,1}^{t,t+h} & 0 & \ldots & 0\\
0 & k_{2,2}^{t,t+h} & \ldots & 0\\
\vdots & \vdots & \ddots & \vdots\\
0 & 0 & \ldots & k_{p,p}^{t,t+h}
\end{array}\right) + \left(
\begin{array}{cccc}
0 & k_{1,2}^{t,t+h} & \ldots & k_{1,p}^{t,t+h}\\
k_{2,1}^{t,t+h} & 0 & \ldots & 0\\
\vdots & \vdots & \ddots & \vdots\\
k_{p,1}^{t,t+h} & k_{p,2}^{t,t+h} & \ldots & 0
\end{array}\right).
}
The entries of the first matrix are called \bold{self-self conditinal dependences} at temporal lag \eqn{h} and represents the (negative) self-similarity of a given random variable across different time points. The entries of the second matrix are the \bold{conditional dependence} among the \eqn{p} random variables. To make the interpretation of the results more relevant and, at the same time, reduce the number of parameters, the authors propose the following equality constraints:
\tabular{clccclcc}{
\tab \eqn{k_{i,i}^{t,t+h}} \tab effect \tab \code{R} code \tab \tab \eqn{k_{i,j}^{t,t+h}} \tab effect \tab \code{R} code \cr
i.	\tab \eqn{0} \tab zero \tab \code{"."} \tab iv. \tab \eqn{0} \tab zero \tab \code{"."} \cr
ii.	\tab \eqn{s^h} \tab costant \tab \code{"c"} \tab ii. \tab \eqn{n^h} \tab costant \tab \code{"c"} \cr
iii.\tab \eqn{s^h_i} \tab unit \tab \code{"u"} \tab iii. \tab \eqn{n^h_i} \tab unit \tab \code{"u"} \cr
iv.	\tab \eqn{s^{t,h}} \tab time \tab \code{"t"} \tab iv. \tab \eqn{n^{t,h}} \tab time \tab \code{"t"} \cr
v.	\tab \eqn{s^{t,h}_i} \tab interaction \tab \code{"ut"} \tab v. \tab \eqn{n^{t,h}_{i,j}} \tab interaction \tab \code{"ut"}
}
Argument \code{model} is used to specify the restrinctions previously describted. This argument can be a named list or a matrix with dimension \eqn{nlag\times 2}, where \eqn{nlag\le\code{tp}}. To gain more insight, suppose that we want to model only the sub-matrices \eqn{K_{t,t}} and \eqn{K_{t,t+1}}, i.e., the sub-matrices corresponding to the temporal lag zero and one. A possible \code{R code} is\cr

\code{model.mat <- matrix("", nrow = 2, ncol = 2)}\cr
\code{rownames(model.mat) <- c("lag0", "lag1")}\cr
\code{colnames(model.mat) <- c("s", "n")}\cr
\code{model.mat[1, ] <- c("c", "ut")}\cr
\code{model.mat[2, ] <- c("t", ".")}\cr

In this example we are modelling the diagonal elements of the sub-matrices \eqn{K_{t,t}} with the constant effect while the off-diagonal elements are modelled by the interaction effect. In the same way, the diagonal elements of the sub-matrices \eqn{K_{t,t+1}} are modelled by the time effect while the remaning elements are equal to zero. The \code{fglasso} function passes the matrix \code{model.mat} to the internal function \code{fglasso_model2mask}, i.e.,\cr 

\code{mask <- fglasso_model2mask(model.mat, tp = 3, p = 3)}\cr

which returns the mask used in \code{sglasso} to fit the specified factorial dynamic Gaussian Graphical model. The same model can be specified by the following named list\cr

\code{model.list <- list(lag0 = c(s = "c", n = "ut"), lag1 = c(s = "t", n = "."))}\cr

See the example below for more details.
}
\value{
\code{fglasso} returns an obejct with S3 class \code{"sglasso"}. See the corresponding manual for more details.
}
\author{Luigi Augugliaro\cr 
Maintainer: Luigi Augugliaro \email{luigi.augugliaro@unipa.it}
}
\references{
Wit, E. C. and Abbruzzo, A.(2015) Dynamic factorial graphical models for dynamic networks. \emph{Network Science}, Vol. \bold{3(1)}, 37-- 57
\cr

Abbruzzo, A., Augugliaro, L., Mineo, A.M. and Wit, E.C. (2014) Cyclic coordinate for penalized Gaussian Graphical Models with symmetry restrictions. In \emph{Proceeding of COMPSTAT 2014 - 21th International Conference on Computational Statistics}. Geneva, August 19-24, 2014.
\cr

Hojsgaard, S. and Lauritzen, S.L. (2008) Graphical gaussian models with edge and vertex symmetries. \emph{J. Roy. Statist. Soc. Ser. B.}, Vol. \bold{70(5)}, 1005--1027.
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
