\name{ordglm}
\title{Generalized Linear Ordinal Regression Models}
\alias{ordglm}
\usage{
ordglm(formula, data=parent.frame(), link="logit", maxiter=10, weights=1)
}
\description{
\code{ordglm} fits linear regression functions with
logistic or probit link to ordinal response data by proportional odds.
}
\arguments{
\item{formula}{A model formula. The response must be integers numbered
from zero to one less than the number of ordered categories.}
\item{data}{An optional data frame containing the variables in the model.}
\item{link}{Logit or probit link function.}
\item{maxiter}{Maximum number of iterations allowed.}
\item{weights}{A vector containing the frequencies for grouped data.}
}
\value{
A list of class ordglm is returned.
The printed output includes the -log likelihood, the corresponding
AIC, the deviance, the maximum likelihood estimates, standard errors,
and correlations.}
\references{
Jansen, J. (1991) Fitting regression models to ordinal
data. Biometrical Journal 33, 807-815.

Johnson, V.E. and Albert, J.H. (1999) Ordinal Data Modeling. Springer-Verlag.
}
\author{J.K. Lindsey, adapted and heavily modified from Matlab code
(ordinalMLE) by J.H. Albert.}
\seealso{
\code{\link{glm}}, \code{\link[gnlm]{nordr}}
}
\examples{
# McCullagh (1980) JRSS B42, 109-142
# tonsil size: 2x3 contingency table
y <- c(0:2,0:2)
carrier <- gl(2,3,6)
wt <- c(19,29,24,497,560,269)
ordglm(y~carrier, weights=wt)
}
\keyword{models}
