% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nordr.r
\name{nordr}
\alias{nordr}
\title{Nonlinear Ordinal Regression Models}
\usage{
nordr(
  y = NULL,
  distribution = "proportional",
  mu = NULL,
  linear = NULL,
  pmu = NULL,
  pintercept = NULL,
  weights = NULL,
  envir = parent.frame(),
  print.level = 0,
  ndigit = 10,
  gradtol = 1e-05,
  steptol = 1e-05,
  fscale = 1,
  iterlim = 100,
  typsize = abs(p),
  stepmax = 10 * sqrt(p \%*\% p)
)
}
\arguments{
\item{y}{A vector of ordinal responses, integers numbered from zero to one
less than the number of categories or an object of class, \code{response}
(created by \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
\code{\link[rmutil]{rmna}}) or \code{\link[rmutil]{lvna}}). If the
\code{repeated} data object contains more than one response variable, give
that object in \code{envir} and give the name of the response variable to be
used here.}

\item{distribution}{The ordinal distribution: proportional odds,
continuation ratio, or adjacent categories.}

\item{mu}{User-specified function of \code{pmu}, and possibly \code{linear},
giving the logistic regression equation. This must contain the first
intercept. It may contain a linear part as the second argument to the
function. It may also be a formula beginning with ~, specifying a logistic
regression function for the location parameter, either a linear one using
the Wilkinson and Rogers notation or a general function with named unknown
parameters. If it contains unknown parameters, the keyword \code{linear} may
be used to specify a linear part. If nothing is supplied, the location is
taken to be constant unless the linear argument is given.}

\item{linear}{A formula beginning with ~ in W&R notation, specifying the
linear part of the logistic regression function.}

\item{pmu}{Vector of initial estimates for the regression parameters,
including the first intercept. If \code{mu} is a formula with unknown
parameters, their estimates must be supplied either in their order of
appearance in the expression or in a named list.}

\item{pintercept}{Vector of initial estimates for the contrasts with the
first intercept parameter (difference in intercept for successive
categories): two less than the number of different ordinal values.}

\item{weights}{Weight vector for use with contingency tables.}

\item{envir}{Environment in which model formulae are to be interpreted or a
data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
name of the response variable should be given in \code{y}. If \code{y} has
class \code{repeated}, it is used as the environment.}

\item{print.level}{Arguments controlling \code{\link{nlm}}.}

\item{ndigit}{Arguments controlling \code{\link{nlm}}.}

\item{gradtol}{Arguments controlling \code{\link{nlm}}.}

\item{steptol}{Arguments controlling \code{\link{nlm}}.}

\item{fscale}{Arguments controlling \code{\link{nlm}}.}

\item{iterlim}{Arguments controlling \code{\link{nlm}}.}

\item{typsize}{Arguments controlling \code{\link{nlm}}.}

\item{stepmax}{Arguments controlling \code{\link{nlm}}.}
}
\value{
A list of class nordr is returned that contains all of the relevant
information calculated, including error codes.
}
\description{
\code{nordr} fits arbitrary nonlinear regression functions (with logistic
link) to ordinal response data by proportional odds, continuation ratio, or
adjacent categories.
}
\details{
Nonlinear regression models can be supplied as formulae where parameters are
unknowns in which case factor variables cannot be used and parameters must
be scalars. (See \code{\link[rmutil]{finterp}}.)

The printed output includes the -log likelihood (not the deviance), the
corresponding AIC, the maximum likelihood estimates, standard errors, and
correlations.
}
\examples{

# McCullagh (1980) JRSS B42, 109-142
# tonsil size: 2x3 contingency table
y <- c(0:2,0:2)
carrier <- c(rep(0,3),rep(1,3))
carrierf <- gl(2,3,6)
wt <- c(19,29,24,
	497,560,269)
pmu <- c(-1,0.5)
mu <- function(p) c(rep(p[1],3),rep(p[1]+p[2],3))
# proportional odds
# with mean function
nordr(y, dist="prop", mu=mu, pmu=pmu, weights=wt, pintercept=1.5)
# using Wilkinson and Rogers notation
nordr(y, dist="prop", mu=~carrierf, pmu=pmu, weights=wt, pintercept=1.5)
# using formula with unknowns
nordr(y, dist="prop", mu=~b0+b1*carrier, pmu=pmu, weights=wt, pintercept=1.5)
# continuation ratio
nordr(y, dist="cont", mu=mu, pmu=pmu, weights=wt, pintercept=1.5)
# adjacent categories
nordr(y, dist="adj", mu=~carrierf, pmu=pmu, weights=wt, pintercept=1.5)
#
# Haberman (1974) Biometrics 30, 589-600
# institutionalized schizophrenics: 3x3 contingency table
y <- rep(0:2,3)
fr <- c(43,6,9,
	16,11,18,
	3,10,16)
length <- gl(3,3)
\dontrun{
# fit continuation ratio model with nordr and as a logistic model
nordr(y, mu=~length, weights=fr, pmu=c(0,-1.4,-2.3), pint=0.13,
	dist="cont")
	}
# logistic regression with reconstructed table
frcr <- cbind(c(43,16,3,49,27,13),c(6,11,10,9,18,16))
lengthord <- gl(3,1,6)
block <- gl(2,3)
summary(glm(frcr~lengthord+block,fam=binomial))
# note that AICs and deviances are different

}
\seealso{
\code{\link[rmutil]{finterp}}, \code{\link[gnlm]{fmr}},
\code{\link{glm}}, \code{\link[repeated]{glmm}},
\code{\link[repeated]{gnlmm}}, \code{\link[gnlm]{gnlr}},
\code{\link[gnlm]{gnlr3}}, \code{\link[gnlm]{nlr}},
\code{\link[gnlm]{ordglm}}
}
\author{
J.K. Lindsey
}
\keyword{models}
