\name{rs2}
\title{Two-factor Box-Tidwell Nonlinear Response Surface Models}
\alias{rs2}
\usage{
rs2(y, x1, x2, power=c(1,1), weight=rep(1,length(y)), family=normal,
	iterlim=20)
}
\description{
\code{rs2} fits a two-covariate power-transformed response surface by
iterating the function, \code{\link{glm}}.}
\value{
A list of class, \code{rs}, is returned containing the model and the power
estimates.
}
\arguments{
\item{y}{Response variable}
\item{x1}{First covariate}
\item{x2}{Second covariate}
\item{power}{Initial estimates of the two power transformations}
\item{weight}{Weight vector}
\item{family}{glm family}
\item{iterlim}{Iteration limit}
}
\author{J.K. Lindsey}
\seealso{
\code{\link{lm}}, \code{\link{glm}}, \code{\link[gnlm]{gnlr}},
\code{\link[gnlm]{gnlr3}}, \code{\link[gnlm]{rs3}}
}
\examples{
x1 <- rep(1:4,5)
x2 <- rep(1:5,rep(4,5))
y <- rpois(20,1+2*sqrt(x1)+3*log(x2)+4*x1+log(x2)^2+2*sqrt(x1)*log(x2))
rs2(y, x1, x2, family=poisson)
}
\keyword{models}
