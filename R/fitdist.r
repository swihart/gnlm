#
#  gnlm : A Library of Special Functions for Nonlinear Regression
#  Copyright (C) 1998, 1999, 2000, 2001 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     fit.dist(y, ni, distribution="normal", breaks=FALSE, delta=1,
#       censor=FALSE, exact=TRUE, plot=FALSE, add=FALSE, main, xlab, ...)
#
#  DESCRIPTION
#
#    A function to fit common distributions to frequency data and to
# plot the resulting curve along with the histogram



#' Fit Probability Distributions to Frequency Data
#' 
#' \code{fit.dist} fits the distributions in Chapter 4 of Lindsey (1995, 2003
#' 2nd edn): binomial, beta-binomial, Poisson, negative binomial, geometric,
#' zeta, normal, log normal, inverse Gauss, logistic, Laplace, Cauchy, Student
#' t, exponential, Pareto, gamma, and Weibull to frequency (histogram) data,
#' possibly plotting the frequency polygon of fitted values with the histogram.
#' 
#' 
#' @param y Vector of observations.
#' @param ni Corresponding vector of frequencies.
#' @param distribution Character string specifying the distribution.
#' @param breaks If TRUE, \code{y} contains breaks between categories instead
#' of mid-points.
#' @param delta Scalar or vector giving the unit of measurement (always one for
#' discrete data) for each response value, set to unity by default. For
#' example, if a response is measured to two decimals, delta=0.01.
#' @param censor If TRUE, the last category is right censored.
#' @param exact If FALSE, uses the approximations for certain distributions in
#' Lindsey (1995).
#' @param plot If TRUE, plots the histogram of observed frequencies and the
#' frequency polygon of fitted values.
#' @param add If TRUE, adds a new frequency polygon of fitted values without
#' replotting the histogram.
#' @param xlab Plotting control options.
#' @param ylab Plotting control options.
#' @param xlim Plotting control options.
#' @param main Plotting control options.
#' @param ... Plotting control options.
#' @author J.K. Lindsey
#' @references Lindsey, J.K. (1995) Introductory Statistics: A Modelling
#' Approach. Oxford: Oxford University Press.
#' @keywords models
#' @examples
#' 
#' f <- c(215, 1485, 5331, 10649, 14959, 11929, 6678, 2092, 342)
#' y <- seq(0,8)
#' fit.dist(y, f, "binomial", plot=TRUE, xlab="Number of males",
#' 	main="Distribution of males in families of 8 children")
#' #
#' f <- c(1,1,6,3,4,3,9,6,5,16,4,11,6,11,3,4,5,6,4,4,5,1,1,4,1,2,
#' 	0,2,0,0,1)
#' y <- seq(1100,4100,by=100)
#' fit.dist(y, f, "normal", delta=100, plot=TRUE,
#' 	xlab="Monthly salary (dollars)",
#' 	main="Distribution of women mathematicians' salaries")
#' fit.dist(y, f, "log normal", delta=100, plot=TRUE, add=TRUE, lty=3)
#' fit.dist(y, f, "logistic", delta=100, exact=FALSE, plot=TRUE, add=TRUE, lty=2)
#' 
#' @export fit.dist
fit.dist <- function(y, ni, distribution="normal", breaks=FALSE, delta=1,
	censor=FALSE, exact=TRUE, plot=FALSE, add=FALSE,
	xlab=deparse(substitute(y)), ylab="Probability", xlim=range(y),
	main=paste("Histogram of",deparse(substitute(y))), ...){
#
# plotting function
#
cum.histo <- function(breaks, freq, prob=FALSE,
	main=main, xlab=xlab, ylab=ylab, xlim=xlim, ...){
	if(prob){
		freq <- freq/(sum(freq)*diff(breaks))
		if(missing(ylab))ylab <- "Relative Frequency"}
	else if(missing(ylab))ylab <- "Frequency"
	if(xlim[1]>min(y)||xlim[2]<max(y)){
		n1 <- which(yi>=xlim[1]&yi<=xlim[2])
		breaks <- breaks[c(n1,n1[length(n1)]+1)]
		freq <- freq[n1]}
	plot(breaks,c(freq,0),type="n",main=main,xlab=xlab,ylab=ylab,
		xlim=c(xlim[1]-0.5,xlim[2]+0.5),...)
	rect(breaks[-length(breaks)],0,breaks[-1],freq,border=par("fg"))}
#
# check that correct data supplied
#
if(!is.vector(y,mode="numeric"))stop("y must be a numeric vector")
if(!is.vector(ni,mode="numeric")||any(ni<0))
	stop("ni must be a numeric vector with non-negative values")
if(any(diff(y)<=0)||all(ni==1))stop("grouped frequency data must be supplied")
n <- length(ni)
if(length(delta)==1)delta <- rep(delta,n)
if(breaks){
	if(length(y)!=n+1)
		stop("Breaks vector must be one longer than frequency vector")
	yi <- (y[1:n]+y[2:(n+1)])/2
	delta <- diff(y)}
else yi <- y
#
# check distribution
#
distribution <- match.arg(distribution,c("binomial","beta binomial",
	"Poisson","negative binomial","geometric","zeta","normal","log normal",
	"inverse Gauss","logistic","Laplace","Cauchy","Student t",
	"exponential","Pareto","gamma","Weibull"))
if(distribution=="zeta"&&any(y<1))stop("y must be positive")
if(distribution=="Cauchy"&&!exact)
	stop("Cauchy distribution can only be fitted exactly")
if(distribution=="Student t"&&!exact)
	stop("Student t distribution can only be fitted exactly")
if(distribution!="normal"&&distribution!="Cauchy"&&distribution!="Laplace"&&
	distribution!="logistic"&&distribution!="Student t"&&any(y<0))
	stop("y must be non-negative")
#
# calculate relative frequencies and empirical mean and variance
#
pi.hat <- ni/sum(ni)
ybar <- weighted.mean(yi,ni)
s2 <- weighted.mean((yi-ybar)^2,ni)
#
# calculate parameter values
#
switch(distribution,
binomial={
	m <- length(yi)-1
	nu <- ybar/m
	lpi.tilde <- dbinom(yi,m,nu,log=TRUE)
	param <- nu
	names(param) <- "nu.hat"
	p <- 1},
"beta binomial"={
	m <- length(yi)-1
	nu <- ybar/m
	gam <- (s2-m*nu*(1-nu))/(m*(m-1)*nu*(1-nu))
	if(exact){
		fcn <- function(p){
			gam1 <- p[1]*(1-p[2])/p[2]
			gam2 <- (1-p[1])*(1-p[2])/p[2]
			-sum(ni*(lbeta(yi+gam1,m-yi+gam2)-lbeta(gam1,gam2)))}
		z <- nlm(fcn,p=c(nu,gam),stepmax=sqrt(nu^2+gam^2)/2,
			 typsize=abs(c(nu,gam)),ndigit=10,gradtol=0.00001,
			 steptol=0.00001,fscale=1)
		nu <- z$estimate[1]
		gam <- z$estimate[2]}
	gam1 <- nu*(1-gam)/gam
	gam2 <- (1-nu)*(1-gam)/gam
	lpi.tilde <- lchoose(m,yi)+lbeta(yi+gam1,m-yi+gam2)-lbeta(gam1,gam2)
	param <- c(nu,gam)
	names(param) <- c("nu.hat","rho.hat")
	p <- 2},
Poisson={
	lpi.tilde <- dpois(yi,ybar,log=TRUE)
	param <- ybar
	names(param) <- "mu.hat"
	p <- 1},
geometric={
	nu <- 1/(1+ybar)
	lpi.tilde <- dnbinom(yi,1,prob=nu,log=TRUE)
	param <- nu
	names(param) <- "nu.hat"
	p <- 1},
"negative binomial"={
	if(ybar>s2)
		stop("underdispersed data not suitable for negative binomial")
	nu <- ybar/s2
	gam <- ybar^2/(s2-ybar)
	if(exact){
		nu <- log(nu/(1-nu))
		fcn <- function(p)
			-sum(ni*dnbinom(yi,p[2],prob=1/(1+exp(-p[1])),log=TRUE))
		z <- nlm(fcn,p=c(nu,gam),stepmax=sqrt(nu^2+gam^2)/2,
			 typsize=abs(c(nu,gam)),ndigit=10,gradtol=0.00001,
			 steptol=0.00001,fscale=1)
		nu <- 1/(1+exp(-z$estimate[1]))
		gam <- z$estimate[2]}
	lpi.tilde <- dnbinom(yi,gam,prob=nu,log=TRUE)
	param <- c(nu,gam)
	names(param) <- c("nu.hat","gamma.hat")
	p <- 2},
zeta={
	pi.tilde <- 1/yi
	nu <- sum(pi.tilde)
	pi.tilde <- pi.tilde/nu
	rho <- round(pi.hat[1]/pi.tilde[1]+0.1)
	if(exact){
		fcn <- function(p) {
			if(censor) const <- sum(yi^(-p[1]))
			else const <- sum(seq(1,30)^(-p[1]))
			sum(ni*(p[1]*log(yi)+log(const)))}
		z <- nlm(fcn, p=rho, stepmax=1)
		rho <- z$estimate[1]}
	lpi.tilde <- -rho*log(yi)
	nu <- sum(exp(lpi.tilde))
	lpi.tilde <- lpi.tilde-log(nu)
	param <- rho
	names(param) <- "rho.hat"
	p <- 1},
normal={
	mu.hat <- ybar
	sigma2.hat <- s2
	lpi.tilde <- dnorm(yi,mu.hat,sqrt(sigma2.hat),log=TRUE)
	param <- c(mu.hat,sigma2.hat)
	names(param) <- c("mu.hat","sigma2.hat")
	p <- 2},
"log normal"={
	mu.hat <- weighted.mean(log(yi),ni)
	sigma2.hat <- weighted.mean((log(yi)-mu.hat)^2,ni)
	lpi.tilde <- dlnorm(yi,mu.hat,sqrt(sigma2.hat),log=TRUE)
	param <- c(mu.hat,sigma2.hat)
	names(param) <- c("mu.hat","sigma2.hat")
	p <- 2},
"inverse Gauss"={
	mu.hat <- ybar
	sigma2.hat <- weighted.mean(1/yi,ni)-(1/ybar)
	lpi.tilde <- -(yi-mu.hat)^2/(2*yi*sigma2.hat*mu.hat^2)-
		0.5*log(2*pi*yi^3*sigma2.hat)
	param <- c(mu.hat,sigma2.hat)
	names(param) <- c("mu.hat","sigma2.hat")
	p <- 2},
logistic={
	mu.hat <- ybar
	sigma <- sqrt(s2)
	if(exact){
		fcn <- function(p)
			-sum(ni*dlogis(yi,p[1],p[2]*sqrt(3)/pi,log=TRUE))
		z <- nlm(fcn, p=c(mu.hat,sigma), stepmax=10)
		mu.hat <- z$estimate[1]
		sigma <- z$estimate[2]}
	lpi.tilde <- dlogis(yi,mu.hat,sigma*sqrt(3)/pi,log=TRUE)
	param <- c(mu.hat,sigma)
	names(param) <- c("mu.hat","sigma.hat")
	p <- 2},
Laplace={
	mu.hat <- yi[which(cumsum(ni)>sum(ni)/2)[1]]
	sigma <- sum(ni*abs(yi-mu.hat))/sum(ni)
	lpi.tilde <- dlaplace(yi,mu.hat,sigma,log=TRUE)
	param <- c(mu.hat,sigma)
	names(param) <- c("mu.hat","sigma.hat")
	p <- 2},
Cauchy={
	mu.hat <- ybar
	sigma <- sqrt(s2)/2
	if(exact){
		fcn <- function(p) -sum(ni*dcauchy(yi,p[1],exp(p[2]),log=TRUE))
		z <- nlm(fcn, p=c(mu.hat,log(sigma)), stepmax=10)
		mu.hat <- z$estimate[1]
		sigma <- exp(z$estimate[2])}
	lpi.tilde <- dcauchy(yi,mu.hat,sigma,log=TRUE)
	param <- c(mu.hat,sigma)
	names(param) <- c("mu.hat","sigma.hat")
	p <- 2},
"Student t"={
	mu.hat <- ybar
	sigma <- sqrt(s2)/2
	df <- 50
	if(exact){
		fcn <- function(p)
			-sum(ni*(dt((yi-p[1])/exp(p[2]),exp(p[3]),log=TRUE)-
			p[2]))
		z <- nlm(fcn, p=c(mu.hat,log(sigma),log(df)), stepmax=10)
		mu.hat <- z$estimate[1]
		sigma <- exp(z$estimate[2])
		df <- exp(z$estimate[3])}
	lpi.tilde <- dt((yi-mu.hat)/sigma,df,log=TRUE)-log(sigma)
	param <- c(mu.hat,sigma,df)
	names(param) <- c("mu.hat","sigma.hat","df")
	p <- 3},
exponential={
	lpi.tilde <- dexp(yi,1/ybar,log=TRUE)
	param <- ybar
	names(param) <- "mu.hat"
	p <- 1},
Pareto={
	delta.hat <- yi[1]-delta[1]/2
	if(delta.hat<=0)stop("smallest observation must be >0")
	alpha.hat <- sum(ni)/sum(ni*log(yi/delta.hat))
	lpi.tilde <- log(alpha.hat)+alpha.hat*log(delta.hat)-
		(alpha.hat+1)*log(yi)
	param <- c(alpha.hat,delta.hat)
	names(param) <- c("alpha.hat","delta.hat")
	p <- 2},
gamma={
	alpha.hat <- ybar^2/s2
	mu.hat <- ybar
	if(exact){
		fcn <- function(p)
			-sum(ni*dgamma(yi,p[2],scale=p[1]/p[2],log=TRUE))
		z <- nlm(fcn, p=c(mu.hat,alpha.hat), stepmax=10)
		mu.hat <- z$estimate[1]
		alpha.hat <- z$estimate[2]}
	lpi.tilde <- dgamma(yi,alpha.hat,scale=mu.hat/alpha.hat,log=TRUE)
	param <- c(alpha.hat,mu.hat)
	names(param) <- c("alpha.hat","mu.hat")
	p <- 2},
Weibull={
	temp <- ybar^2/(s2+ybar^2)
	Alpha.Weibull.fn <- function(y){
		alpha.trans.fn <- function(al)
			gamma(1+1/al)*gamma(1+1/al)/gamma(1+2/al)
		tol <- 0.001
		al.start <- 0.0001
		al.end <- 50
		al.mid <- 0.5*(al.start+al.end)
		y.tmp <- alpha.trans.fn(al.mid)
		while (abs(y.tmp-y)>tol){
			if ((y.tmp-y)>0) al.end <- al.mid
			else al.start <- al.mid
			al.mid <- 0.5*(al.start+al.end)
			y.tmp <- alpha.trans.fn(al.mid)}
		al.mid}
	alpha.hat <- Alpha.Weibull.fn(temp)
	mu.hat <- ybar/gamma(1+1/alpha.hat)
	if(exact){
		fcn <- function(p)
			-sum(ni*dweibull(yi,p[2],p[1],log=TRUE))
		z <- nlm(fcn, p=c(mu.hat,alpha.hat), stepmax=10)
		mu.hat <- z$estimate[1]
		alpha.hat <- z$estimate[2]}
	lpi.tilde <- dweibull(yi,alpha.hat,mu.hat,log=TRUE)
	param <- c(alpha.hat,mu.hat)
	names(param) <- c("alpha.hat","mu.hat")
	p <- 2})
#
# calculate theoretical probabilities
#
lpi.tilde <- lpi.tilde+log(delta)
pi.tilde <- exp(lpi.tilde)
if(censor){
	pi.tilde[length(pi.tilde)] <- 1-sum(pi.tilde[1:(length(pi.tilde)-1)])
	if(pi.tilde[length(pi.tilde)]==0)pi.tilde[length(pi.tilde)] <- 1e-20}
#
# calculate likelihood function, AIC, and residuals
#
like.comp <- rep(0,length(pi.tilde))
like.comp[ni>0] <- -ni[ni>0]*(log(pi.tilde[ni>0])-log(pi.hat[ni>0]))
loglike <- sum(like.comp)
AIC <- loglike+p
resid <- as.vector((ni-sum(ni)*pi.tilde)/sqrt(sum(ni)*pi.tilde))
#
# print results
#
result1.output <- c(ybar,s2,param)
result2.output <- c(loglike,AIC)
names(result1.output) <- c("mean","variance",names(param))
names(result2.output) <- c("-log likelihood","AIC")
cat("\n",distribution," distribution,","  n = ",sum(ni),"\n\n",sep="")
print(result1.output)
cat("\n")
print(result2.output)
cat("\n")
#
# if required, plot results
#
if(plot){
	n1 <- if(xlim[1]>min(y)||xlim[2]<max(y))which(yi>=xlim[1]&yi<=xlim[2])
		 else 1:n
	if(censor&&xlim[2]==max(y))n1 <- n1[-length(n1)]
	if(add)lines(yi[n1],pi.tilde[n1]/delta[n1],...)
	else {
		cum.histo(c(yi-delta/2,yi[n]+delta[n]/2),
			ni,prob=TRUE,main,xlab,ylab=ylab,xlim,...)
		lines(yi[n1], pi.tilde[n1]/delta[n1],...)}}
#
# return dataframe of results
#
if(length(unique(delta))==1)
	data.frame(yi,ni,pi.hat,pi.tilde,like.comp,resid)
else data.frame(yi,ni,delta,pi.hat,pi.tilde,like.comp,resid)}
