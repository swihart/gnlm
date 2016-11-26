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
#     nlr(y=NULL, mu=NULL, pmu=NULL, distribution="normal", wt=1, delta=1,
#	envir=parent.frame(), print.level=0, typsize=abs(pmu),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(pmu%*%pmu),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models for distributions
# in the exponential family.



#' Nonlinear Normal, Gamma, and Inverse Gaussian Regression Models
#' 
#' \code{nlr} fits a user-specified nonlinear regression equation by least
#' squares (normal) or its generalization for the gamma and inverse Gauss
#' distributions.
#' 
#' A nonlinear regression model can be supplied as a formula where parameters
#' are unknowns in which case factor variables cannot be used and parameters
#' must be scalars. (See \code{\link[rmutil]{finterp}}.)
#' 
#' The printed output includes the -log likelihood (not the deviance), the
#' corresponding AIC, the parameter estimates, standard errors, and
#' correlations.
#' 
#' 
#' @param y The response vector or an object of class, \code{response} (created
#' by \code{\link[rmutil]{restovec}}) or \code{repeated} (created by
#' \code{\link[rmutil]{rmna}} or \code{\link[rmutil]{lvna}}).
#' @param mu A function of \code{p} giving the regression equation for the mean
#' or a formula beginning with ~, specifying either a linear regression
#' function in the Wilkinson and Rogers notation or a general nonlinear
#' function with named unknown parameters.
#' @param pmu Vector of initial estimates of the parameters. If \code{mu} is a
#' formula with unknown parameters, their estimates must be supplied either in
#' their order of appearance in the expression or in a named list.
#' @param distribution The distribution to be used: normal, gamma, or inverse
#' Gauss.
#' @param wt Weight vector.
#' @param delta Scalar or vector giving the unit of measurement for each
#' response value, set to unity by default. For example, if a response is
#' measured to two decimals, \code{delta=0.01}. If the response is transformed,
#' this must be multiplied by the Jacobian. For example, with a log
#' transformation, \code{delta=1/y}.
#' @param envir Environment in which model formulae are to be interpreted or a
#' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}. If
#' \code{y} has class \code{repeated}, it is used as the environment.
#' @param print.level Arguments controlling \code{\link{nlm}}.
#' @param typsize Arguments controlling \code{\link{nlm}}.
#' @param ndigit Arguments controlling \code{\link{nlm}}.
#' @param gradtol Arguments controlling \code{\link{nlm}}.
#' @param stepmax Arguments controlling \code{\link{nlm}}.
#' @param steptol Arguments controlling \code{\link{nlm}}.
#' @param iterlim Arguments controlling \code{\link{nlm}}.
#' @param fscale Arguments controlling \code{\link{nlm}}.
#' @return A list of class \code{nlr} is returned that contains all of the
#' relevant information calculated, including error codes.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{finterp}}, \code{\link[gnlm]{fmr}},
#' \code{\link{glm}}, \code{\link[repeated]{glmm}},
#' \code{\link[repeated]{gnlmm}}, \code{\link[gnlm]{gnlr}},
#' \code{\link[gnlm]{gnlr3}}, \code{\link{lm}}, \code{\link[stats]{nls}}.
#' @keywords models
#' @examples
#' 
#' x <- c(3,5,0,0,0,3,2,2,2,7,4,0,0,2,2,2,0,1,3,4)
#' y <- c(5.8,11.6,2.2,2.7,2.3,9.4,11.7,3.3,1.5,14.6,9.6,7.4,10.7,6.9,
#' 	2.6,17.3,2.8,1.2,1.0,3.6)
#' # rgamma(20,2,scale=0.2+2*exp(0.1*x))
#' # linear least squares regression
#' mu1 <- function(p) p[1]+p[2]*x
#' summary(lm(y~x))
#' nlr(y, mu=mu1, pmu=c(3,0))
#' # or
#' nlr(y, mu=~x, pmu=c(3,0))
#' # or
#' nlr(y, mu=~b0+b1*x, pmu=c(3,0))
#' # linear gamma regression
#' nlr(y, dist="gamma", mu=~x, pmu=c(3,0))
#' # nonlinear regression
#' mu2 <- function(p) p[1]+p[2]*exp(p[3]*x)
#' nlr(y, mu=mu2, pmu=c(0.2,3,0.2))
#' # or
#' nlr(y, mu=~b0+c0*exp(c1*x), pmu=list(b0=0.2,c0=3,c1=0.2))
#' # with gamma distribution
#' nlr(y, dist="gamma", mu=~b0+c0*exp(c1*x), pmu=list(b0=0.2,c0=3,c1=0.2))
#' 
#' @export nlr
nlr <- function(y=NULL, mu=NULL, pmu=NULL, distribution="normal", wt=1,
	delta=1, envir=parent.frame(), print.level=0, typsize=abs(pmu),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(pmu%*%pmu),
	steptol=0.00001, iterlim=100, fscale=1){
call <- sys.call()
#
# check distribution, data, and initial parameter estimates
#
distribution <- match.arg(distribution,c("normal","inverse Gauss","gamma"))
if(is.null(pmu))stop("Initial parameter estimates must be supplied")
np <- length(pmu)
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("nlr only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("nlr does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL
#
# if a data object was supplied, modify formula or function to read from it
#
mu2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	else if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
		        tmp <- parse(text=deparse(mu)[-1])
		        mu <- if(respenv)fnenvir(mu,.envir=y,.name=envname)
		        	else fnenvir(mu,.envir=envir,.name=envname)
		        mu2 <- mu
			attr(mu2,"model") <- tmp}
		else mu2 <- mu}}
#
# transform formula to function and check number of parameters
#
if(inherits(mu,"formula")){
	mu <- if(respenv)finterp(mu,.envir=y,.name=envname)
		else finterp(mu,.envir=envir,.name=envname)
	npt1 <- length(attr(mu,"parameters"))
	if(is.character(attr(mu,"model"))){
	# W&R formula
		if(length(attr(mu,"model"))==1){
		# intercept model
			tmp <- attributes(mu)
			mu <- function(p) p[1]*rep(1,n)
			attributes(mu) <- tmp}}
	else {
	# formula with unknowns
		if(np!=npt1){
			cat("\nParameters are ")
			cat(attr(mu,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}}
#
# give appropriate attributes to mu for printing
#
if(is.null(mu)||!is.function(mu))
	stop("A mean function or formula must be supplied")
if(is.null(attr(mu,"parameters"))){
	attributes(mu) <- if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=y))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&any(!is.na(y$response$wt)))
		wt <- as.vector(y$response$wt)
	if(!is.null(y$response$delta))
		delta <- as.vector(y$response$delta)
	type <- y$response$type
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("nlr does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$wt))wt <- as.vector(envir$response$wt)
	if(!is.null(envir$response$delta))
		delta <- as.vector(envir$response$delta[,col])}
else if(inherits(envir,"data.frame"))y <- envir[[deparse(substitute(y))]]
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("nlr only handles univariate responses")
	if(!is.null(y$wt)&&any(!is.na(y$wt)))wt <- as.vector(y$wt)
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	type <- y$type
	y <- response(y)}
if(distribution=="normal"&&type!="unknown"&&type!="continuous"&&
	type!="duration")stop("continuous data required")
if((distribution=="gamma"||distribution=="inverse Gauss")&&
	type!="unknown"&&type!="continuous"&&type!="duration")
	stop("duration data required")
if(!is.vector(y,mode="numeric"))stop("y must be a vector")
if(any(is.na(y)))stop("NAs in y - use rmna")
n <- length(y)
if(length(delta)==1)delta <- rep(delta,length(y))
if(length(wt)==1)wt <- rep(wt,length(y))
#
# check that correct number of estimates was supplied
#
nlp <- if(is.function(mu))length(attr(mu,"parameters"))
       else npt1
if(nlp!=np)stop(paste("pmu should have",nlp,"initial estimates"))
#
# set up likelihood function without dispersion parameter and optimize
#
fn <- switch(distribution,
	normal=function(p) sum(wt*(y-mu(p))^2),
	gamma=function(p) -sum(wt*(log(y/mu(p))-(y-mu(p))/mu(p))),
	"inverse Gauss"=function(p) sum(wt*((y-mu(p))^2)/(y*mu(p)^2)))
if(fscale==1)fscale <- fn(pmu)
if(is.na(fn(pmu)))
	stop("Non-numerical function value: probably invalid initial values")
z0 <- nlm(fn, p=pmu, hessian=TRUE, print.level=print.level, typsize=typsize,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
n <- sum(wt)
disp <- z0$minimum/n
p <- z0$estimate
#
# calculate log likelihood
#
switch(distribution,
	normal=maxlike <- length(y)*(log(2*pi*disp)+1)/2,
	gamma=maxlike <- (sum(wt*(y/mu(p)+log(mu(p))-log(y)))+n*log(disp))/
		disp+n*lgamma(1/disp)+sum(log(y)*wt),
	"inverse Gauss"=maxlike <- (sum(wt)*(log(disp*2*pi)+1)+
		 3*sum(log(y)*wt))/2)
maxlike <- maxlike-sum(log(delta))
#
# calculate fitted values and residuals
#
fitted.values <-  as.vector(mu(z0$estimate))
residuals <-  y-fitted.values
#
# calculate se's
#
if(np==1)cov <- 1/z0$hessian
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)}
cov <- 2*cov*z0$minimum/sum(wt)
se <- sqrt(diag(cov))
z1 <- list(
	call=call,
	distribution=distribution,
	delta=delta,
	mu=mu,
	prior.weights=wt,
	maxlike=maxlike,
	dispersion=disp,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=maxlike+np+1,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	np=np,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "nlr"
return(z1)}

### standard methods
###

weights.nlr <- function(z) z$prior.weights

df.residual.nlr <- function(z) z$df

deviance.nlr <- function(z) 2*z$maxlike

### print method
###
print.nlr <- function(z,digits=max(4,.Options$digits-3),correlation=TRUE){
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat(z$distribution,"distribution\n\n")
cat("Mean function:\n")
if(!is.null(attr(z$mu,"formula")))cat(deparse(attr(z$mu,"formula")),sep="\n")
else if(!is.null(attr(z$mu,"model"))){
	t <- deparse(attr(z$mu,"model"))
	t[1] <- sub("expression\\(","",t[1])
	t[length(t)] <- sub("\\)$","",t[length(t)])
	cat(t,sep="\n")}
cat("\n-Log likelihood   ",z$maxlike,"\n")
cat("Degrees of freedom",z$df,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n\n")
cat("Mean parameters:\n")
coef.table <- cbind(z$coefficients[1:z$np], z$se[1:z$np])
if(inherits(z$mu,"formulafn"))
	cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
		else attr(z$mu,"parameters")
else cname <- seq(1,z$np)
dimnames(coef.table) <- list(cname, c("estimate", "se"))
print.default(coef.table, digits=digits, print.gap=2)
cat("\nDispersion estimate:",z$dispersion,"\n")
if(z$np>1&&correlation){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,z$np),seq(1,z$np))
	print.default(z$corr, digits=digits)}
invisible(z)}
