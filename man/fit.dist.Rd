% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitdist.r
\name{fit.dist}
\alias{fit.dist}
\title{Fit Probability Distributions to Frequency Data}
\usage{
fit.dist(
  y,
  ni,
  distribution = "normal",
  breaks = FALSE,
  delta = 1,
  censor = FALSE,
  exact = TRUE,
  plot = FALSE,
  add = FALSE,
  xlab = deparse(substitute(y)),
  ylab = "Probability",
  xlim = range(y),
  main = paste("Histogram of", deparse(substitute(y))),
  ...
)
}
\arguments{
\item{y}{Vector of observations.}

\item{ni}{Corresponding vector of frequencies.}

\item{distribution}{Character string specifying the distribution.}

\item{breaks}{If TRUE, \code{y} contains breaks between categories instead
of mid-points.}

\item{delta}{Scalar or vector giving the unit of measurement (always one for
discrete data) for each response value, set to unity by default. For
example, if a response is measured to two decimals, delta=0.01.}

\item{censor}{If TRUE, the last category is right censored.}

\item{exact}{If FALSE, uses the approximations for certain distributions in
Lindsey (1995).}

\item{plot}{If TRUE, plots the histogram of observed frequencies and the
frequency polygon of fitted values.}

\item{add}{If TRUE, adds a new frequency polygon of fitted values without
replotting the histogram.}

\item{xlab}{Plotting control options.}

\item{ylab}{Plotting control options.}

\item{xlim}{Plotting control options.}

\item{main}{Plotting control options.}

\item{...}{Plotting control options.}
}
\description{
\code{fit.dist} fits the distributions in Chapter 4 of Lindsey (1995, 2003
2nd edn): binomial, beta-binomial, Poisson, negative binomial, geometric,
zeta, normal, log normal, inverse Gauss, logistic, Laplace, Cauchy, Student
t, exponential, Pareto, gamma, and Weibull to frequency (histogram) data,
possibly plotting the frequency polygon of fitted values with the histogram.
}
\examples{

f <- c(215, 1485, 5331, 10649, 14959, 11929, 6678, 2092, 342)
y <- seq(0,8)
fit.dist(y, f, "binomial", plot=TRUE, xlab="Number of males",
	main="Distribution of males in families of 8 children")
#
f <- c(1,1,6,3,4,3,9,6,5,16,4,11,6,11,3,4,5,6,4,4,5,1,1,4,1,2,
	0,2,0,0,1)
y <- seq(1100,4100,by=100)
fit.dist(y, f, "normal", delta=100, plot=TRUE,
	xlab="Monthly salary (dollars)",
	main="Distribution of women mathematicians' salaries")
fit.dist(y, f, "log normal", delta=100, plot=TRUE, add=TRUE, lty=3)
fit.dist(y, f, "logistic", delta=100, exact=FALSE, plot=TRUE, add=TRUE, lty=2)

}
\references{
Lindsey, J.K. (1995) Introductory Statistics: A Modelling
Approach. Oxford: Oxford University Press.
}
\author{
J.K. Lindsey
}
\keyword{models}
