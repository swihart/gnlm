% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnlr3.r
\name{gnlr3}
\alias{gnlr3}
\title{Generalized Nonlinear Regression Models for Three Parameter Distributions}
\usage{
gnlr3(
  y = NULL,
  distribution = "normal",
  mu = NULL,
  shape = NULL,
  family = NULL,
  linear = NULL,
  pmu = NULL,
  pshape = NULL,
  pfamily = NULL,
  exact = FALSE,
  wt = 1,
  common = FALSE,
  delta = 1,
  envir = parent.frame(),
  print.level = 0,
  typsize = abs(p),
  ndigit = 10,
  gradtol = 1e-05,
  stepmax = 10 * sqrt(p \%*\% p),
  steptol = 1e-05,
  iterlim = 100,
  fscale = 1
)
}
\arguments{
\item{y}{The response vector for uncensored data, two columns for censored
data, with the second being the censoring indicator (1: uncensored, 0: right
censored, -1: left censored.), or an object of class, \code{response}
(created by \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
\code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}). If the
\code{repeated} data object contains more than one response variable, give
that object in \code{envir} and give the name of the response variable to be
used here.}

\item{distribution}{Either a character string containing the name of the
distribution or a function giving the -log likelihood and calling the
location, shape, and family functions. Distributions are Box-Cox transformed
normal, generalized inverse Gauss, generalized logistic, Hjorth, generalized
gamma, Burr, generalized Weibull, power exponential, Student t, generalized
extreme value, power variance function Poisson, and skew Laplace. (For
definitions of distributions, see the corresponding [dpqr]distribution
help.)}

\item{mu}{A user-specified function of \code{pmu}, and possibly
\code{linear}, giving the regression equation for the location. This may
contain a linear part as the second argument to the function. It may also be
a formula beginning with ~, specifying either a linear regression function
for the location parameter in the Wilkinson and Rogers notation or a general
function with named unknown parameters. If it contains unknown parameters,
the keyword \code{linear} may be used to specify a linear part. If nothing
is supplied, the location is taken to be constant unless the linear argument
is given.}

\item{shape}{A user-specified function of \code{pshape}, and possibly
\code{linear}, giving the regression equation for the dispersion or shape
parameter. This may contain a linear part as the second argument to the
function. It may also be a formula beginning with ~, specifying either a
linear regression function for the shape parameter in the Wilkinson and
Rogers notation or a general function with named unknown parameters. If it
contains unknown parameters, the keyword \code{linear} may be used to
specify a linear part. If nothing is supplied, this parameter is taken to be
constant unless the linear argument is given. This parameter is the
logarithm of the usual one.}

\item{family}{A user-specified function of \code{pfamily}, and possibly
\code{linear}, for the regression equation of the third (family) parameter
of the distribution. This may contain a linear part that is the second
argument to the function. It may also be a formula beginning with ~,
specifying either a linear regression function for the family parameter in
the Wilkinson and Rogers notation or a general function with named unknown
parameters. If neither is supplied, this parameter is taken to be constant
unless the linear argument is given. In most cases, this parameter is the
logarithm of the usual one.}

\item{linear}{A formula beginning with ~ in W&R notation, specifying the
linear part of the regression function for the location parameters or list
of three such expressions for the location, shape, and/or family parameters.}

\item{pmu}{Vector of initial estimates for the location parameters. If
\code{mu} is a formula with unknown parameters, their estimates must be
supplied either in their order of appearance in the expression or in a named
list.}

\item{pshape}{Vector of initial estimates for the shape parameters. If
\code{shape} is a formula with unknown parameters, their estimates must be
supplied either in their order of appearance in the expression or in a named
list.}

\item{pfamily}{Vector of initial estimates for the family parameters. If
\code{family} is a formula with unknown parameters, their estimates must be
supplied either in their order of appearance in the expression or in a named
list.}

\item{exact}{If TRUE, fits the exact likelihood function for continuous data
by integration over intervals of observation given in \code{delta}, i.e.
interval censoring.}

\item{wt}{Weight vector.}

\item{common}{If TRUE, at least two of \code{mu}, \code{shape}, and
\code{family} must both be either functions with, as argument, a vector of
parameters having some or all elements in common between them so that
indexing is in common between them or formulae with unknowns. All parameter
estimates must be supplied in \code{pmu}. If FALSE, parameters are distinct
between the two functions and indexing starts at one in each function.}

\item{delta}{Scalar or vector giving the unit of measurement (always one for
discrete data) for each response value, set to unity by default - for
example, if a response is measured to two decimals, \code{delta=0.01}. If
the response is transformed, this must be multiplied by the Jacobian. The
transformation cannot contain unknown parameters. For example, with a log
transformation, \code{delta=1/y}. (The delta values for the censored
response are ignored.)}

\item{envir}{Environment in which model formulae are to be interpreted or a
data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
name of the response variable should be given in \code{y}. If \code{y} has
class \code{repeated}, it is used as the environment.}

\item{print.level}{Arguments controlling \code{\link{nlm}}.}

\item{typsize}{Arguments controlling \code{\link{nlm}}.}

\item{ndigit}{Arguments controlling \code{\link{nlm}}.}

\item{gradtol}{Arguments controlling \code{\link{nlm}}.}

\item{stepmax}{Arguments controlling \code{\link{nlm}}.}

\item{steptol}{Arguments controlling \code{\link{nlm}}.}

\item{iterlim}{Arguments controlling \code{\link{nlm}}.}

\item{fscale}{Arguments controlling \code{\link{nlm}}.}
}
\value{
A list of class \code{gnlm} is returned that contains all of the
relevant information calculated, including error codes.
}
\description{
\code{gnlr3} fits user specified nonlinear regression equations to one, two,
or all three parameters of three parameter distributions. Continuous data
may be left, right, and/or interval censored.
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

sex <- c(rep(0,10),rep(1,10))
sexf <- gl(2,10)
age <- c(8,10,12,12,8,7,16,7,9,11,8,9,14,12,12,11,7,7,7,12)
y <- cbind(c(9.2, 7.3,13.0, 6.9, 3.9,14.9,17.8, 4.8, 6.4, 3.3,17.2,
	14.4,17.0, 5.0,17.3, 3.8,19.4, 5.0, 2.0,19.0),
	c(0,1,0,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1))
# y <- cbind(rweibull(20,2,2+2*sex+age),rbinom(20,1,0.7))
# log linear regression with the generalized Weibull distribution
mu <- function(p) exp(p[1]+p[2]*sex+p[3]*age)
gnlr3(y, dist="Weibull", mu=mu, pmu=c(3,1,0), pshape=2, pfamily=-2)
# or equivalently
mu1 <- function(p,linear) exp(linear)
gnlr3(y, dist="Weibull", mu=mu1, linear=~sex+age, pmu=c(3,1,0),
	pshape=2, pfamily=-2)
# or
gnlr3(y, dist="Weibull", mu=~exp(b0+b1*sex+b2*age),
	pmu=list(b0=3,b1=1,b2=0), pshape=2, pfamily=-2)
#
# include regression for the shape parameter with same mu function
shape <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr3(y, dist="Weibull", mu=mu, shape=shape,
	pmu=c(3,1,0), pshape=c(2,0,0), pfamily=-2)
# or equivalently
gnlr3(y, dist="Weibull", mu=mu1, linear=list(~sexf+age,~sex+age,NULL),
	pmu=c(3,1,0), pshape=c(2,0,0), pfamily=-2)
# or
gnlr3(y, dist="Weibull", mu=~exp(b0+b1*sex+b2*age),
	shape=~c0+c1*sex+c2*age, pmu=c(3,1,0),
	pshape=list(c0=2,c1=0,c2=0), pfamily=-2)
# include regression for the family parameter with same mu
# and shape functions
family <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr3(y, dist="Weibull", mu=mu1, linear=~sexf+age, shape=shape,
	family=family, pmu=c(2.5,1,0), pshape=c(2,0,0), pfamily=c(-2,0,0))
# or equivalently
gnlr3(y, dist="Weibull", mu=mu1, linear=list(~sex+age,~sex+age,~sex+age),
	pmu=c(2.5,1,0), pshape=c(2,0,0), pfamily=c(-2,0,0))
# or
gnlr3(y, dist="Weibull", mu=~exp(b0+b1*sex+b2*age),
	shape=~c0+c1*sex+c2*age, family=~d0+d1*sex+d2*age,
	pmu=list(b0=2.5,b1=1,b2=0), pshape=list(c0=2,c1=0,c2=0),
	pfamily=list(d0=-2,d1=0,d2=0))
#
# common parameters
mu <- function(p) exp(p[1]+p[2]*sex+p[3]*age)
shape <- function(p) p[4]+p[5]*sex+p[3]*age
family <- function(p) p[6]+p[7]*sex+p[3]*age
gnlr3(y, dist="Weibull", mu=mu, shape=shape, family=family,
	pmu=c(2.5,1,0,1,0,1,0), common=TRUE)
# or
gnlr3(y, dist="Weibull", mu=~exp(a+b*sex+c*age), shape=~d+e*sex+c*age,
	family=~f+g*sex+c*age, pmu=c(2.5,1,0,1,0,1,0), common=TRUE)

}
\seealso{
\code{\link[rmutil]{finterp}}, \code{\link[gnlm]{fmr}},
\code{\link{glm}}, \code{\link[gnlm]{gnlr}}, \code{\link{lm}},
\code{\link[gnlm]{nlr}}, \code{\link[stats]{nls}}.
}
\author{
J.K. Lindsey
}
\keyword{models}
