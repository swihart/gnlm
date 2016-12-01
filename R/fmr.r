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
#     fmr(y=NULL, distribution="normal", mu=NULL, mix=NULL, linear=NULL,
#	pmu=NULL, pmix=NULL, pshape=NULL, censor="right", exact=FALSE,
#	wt=1, delta=1, common=FALSE, envir=parent.frame(),
#	print.level=0, typsize=abs(p), ndigit=10, gradtol=0.00001,
#	stepmax=10*sqrt(p%*%p), steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with a variety of
# distributions and a mixture in the tail(s).



#' Generalized Nonlinear Regression Models with Two or Three Point Mixtures
#' 
#' \code{fmr} fits user specified nonlinear regression equations to the
#' location parameter of the common one and two parameter distributions. (The
#' log of the scale parameter is estimated to ensure positivity.)
#' 
#' For the Poisson and related distributions, the mixture involves the zero
#' category. For the binomial and related distributions, it involves the two
#' extreme categories. For all other distributions, it involves either left or
#' right censored individuals. A user-specified -log likelihood can also be
#' supplied for the distribution.
#' 
#' Nonlinear regression models can be supplied as formulae where parameters are
#' unknowns in which case factor variables cannot be used and parameters must
#' be scalars. (See \code{\link[rmutil]{finterp}}.)
#' 
#' The printed output includes the -log likelihood (not the deviance), the
#' corresponding AIC, the maximum likelihood estimates, standard errors, and
#' correlations.
#' 
#' 
#' @param y A response vector for uncensored data, a two column matrix for
#' binomial data or censored data, with the second column being the censoring
#' indicator (1: uncensored, 0: right censored, -1: left censored), or an
#' object of class, \code{response} (created by \code{\link[rmutil]{restovec}})
#' or \code{repeated} (created by \code{\link[rmutil]{rmna}} or
#' \code{\link[rmutil]{lvna}}). If the \code{repeated} data object contains
#' more than one response variable, give that object in \code{envir} and give
#' the name of the response variable to be used here.
#' @param distribution Either a character string containing the name of the
#' distribution or a function giving the -log likelihood and calling the
#' location and mixture functions. Distributions are binomial, beta binomial,
#' double binomial, multiplicative binomial, Poisson, negative binomial, double
#' Poisson, multiplicative Poisson, gamma count, Consul, geometric, normal,
#' inverse Gauss, logistic, exponential, gamma, Weibull, extreme value, Pareto,
#' Cauchy, Student t, Laplace, and Levy. (For definitions of distributions, see
#' the corresponding [dpqr]distribution help.)
#' @param mu A user-specified function of \code{pmu}, and possibly
#' \code{linear}, giving the regression equation for the location. This may
#' contain a linear part as the second argument to the function. It may also be
#' a formula beginning with ~, specifying either a linear regression function
#' for the location parameter in the Wilkinson and Rogers notation or a general
#' function with named unknown parameters. If it contains unknown parameters,
#' the keyword \code{linear} may be used to specify a linear part. If nothing
#' is supplied, the location is taken to be constant unless the linear argument
#' is given.
#' @param mix A user-specified function of \code{pmix}, and possibly
#' \code{linear}, giving the regression equation for the mixture parameter.
#' This may contain a linear part as the second argument to the function. It
#' may also be a formula beginning with ~, specifying either a linear
#' regression function for the mixture parameter in the Wilkinson and Rogers
#' notation or a general function with named unknown parameters. If it contains
#' unknown parameters, the keyword \code{linear} may be used to specify a
#' linear part. If nothing is supplied, this parameter is taken to be constant.
#' This parameter is the logit of the mixture probability.
#' @param linear A formula beginning with ~ in W&R notation, or list of two
#' such expressions, specifying the linear part of the regression function for
#' the location or location and mixture parameters.
#' @param pmu Vector of initial estimates for the location parameters. If
#' \code{mu} is a formula with unknown parameters, their estimates must be
#' supplied either in their order of appearance in the expression or in a named
#' list.
#' @param pshape An initial estimate for the shape parameter.
#' @param pmix Vector of initial estimates for the mixture parameters. If
#' \code{mix} is a formula with unknown parameters, their estimates must be
#' supplied either in their order of appearance in the expression or in a named
#' list.
#' @param censor \code{right}, \code{left}, or \code{both} indicating where the
#' mixing distribution is placed. \code{both} is only possible for binomial
#' data.
#' @param exact If TRUE, fits the exact likelihood function for continuous data
#' by integration over intervals of observation given in \code{delta}, i.e.
#' interval censoring.
#' @param wt Weight vector.
#' @param delta Scalar or vector giving the unit of measurement (always one for
#' discrete data) for each response value, set to unity by default - for
#' example, if a response is measured to two decimals, \code{delta=0.01}. If
#' the response is transformed, this must be multiplied by the Jacobian. The
#' transformation cannot contain unknown parameters. For example, with a log
#' transformation, \code{delta=1/y}.
#' @param common If TRUE, \code{mu} and \code{mix} must both be either
#' functions with, as argument, a vector of parameters having some or all
#' elements in common between them so that indexing is in common between them
#' or formulae with unknowns. All parameter estimates must be supplied in
#' \code{pmu}. If FALSE, parameters are distinct between the two functions and
#' indexing starts at one in each function.
#' @param envir Environment in which model formulae are to be interpreted or a
#' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
#' name of the response variable should be given in \code{y}. If \code{y} has
#' class \code{repeated}, it is used as the environment.
#' @param print.level Arguments controlling \code{\link{nlm}}.
#' @param typsize Arguments controlling \code{\link{nlm}}.
#' @param ndigit Arguments controlling \code{\link{nlm}}.
#' @param gradtol Arguments controlling \code{\link{nlm}}.
#' @param stepmax Arguments controlling \code{\link{nlm}}.
#' @param steptol Arguments controlling \code{\link{nlm}}.
#' @param iterlim Arguments controlling \code{\link{nlm}}.
#' @param fscale Arguments controlling \code{\link{nlm}}.
#' @return A list of class \code{gnlm} is returned that contains all of the
#' relevant information calculated, including error codes.
#' @author J.K. Lindsey
#' @seealso \code{\link[rmutil]{finterp}}, \code{\link{glm}},
#' \code{\link[gnlm]{gnlr}}, \code{\link[gnlm]{gnlr3}}, \code{\link{lm}}.
#' @keywords models
#' @examples
#' 
#' sex <- c(rep(0,10),rep(1,10))
#' sexf <- gl(2,10)
#' age <- c(8,10,12,12,8,7,16,7,9,11,8,9,14,12,12,11,7,7,7,12)
#' y <- cbind(c(9.2, 7.3,13.0, 6.9, 3.9,14.9,17.8, 4.8, 6.4, 3.3,17.2,
#' 	14.4,17.0, 5.0,17.3, 3.8,19.4, 5.0, 2.0,19.0),
#' 	c(0,1,0,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1))
#' # y <- cbind(rweibull(20,2,2+2*sex+age),rbinom(20,1,0.7))
#' # log linear regression with Weibull distribution with a point mass
#' #   for right censored individuals
#' mu <- function(p) exp(p[1]+p[2]*sex+p[3]*age)
#' fmr(y, dist="Weibull", mu=mu, pmu=c(4,0,0), pmix=0.5, pshape=1)
#' # or equivalently
#' fmr(y, dist="Weibull", mu=function(p,linear) exp(linear),
#' 	linear=~sexf+age, pmu=c(4,0,0), pmix=0.5, pshape=1)
#' # or
#' fmr(y, dist="Weibull", mu=~exp(b0+b1*sex+b2*age), pmu=list(b0=4,b1=0,b2=0),
#' 	pmix=0.5, pshape=1)
#' #
#' # include logistic regression for the mixture parameter
#' mix <- function(p) p[1]+p[2]*sex
#' fmr(y, dist="Weibull", mu=~exp(a+b*age), mix=mix, pmu=c(4,0),
#' 	pmix=c(10,0), pshape=0.5)
#' # or equivalently
#' fmr(y, dist="Weibull", mu=function(p,linear) exp(linear),
#' 	linear=list(~age,~sexf), pmu=c(4,0), pmix=c(10,0), pshape=0.5)
#' # or
#' fmr(y, dist="Weibull", mu=~exp(b0+b1*age), mix=~c0+c1*sex,
#' 	pmu=list(b0=4,b1=0), pmix=list(c0=10,c1=0), pshape=0.5)
#' #
#' # generate zero-inflated negative binomial data
#' x1 <- rpois(50,4)
#' x2 <- rpois(50,4)
#' ind <- rbinom(50,1,1/(1+exp(-1-0.1*x1)))
#' y <- ifelse(ind,rnbinom(50,3,mu=exp(1+0.2*x2)),0)
#' # standard Poisson models
#' gnlr(y, dist="Poisson", mu=~exp(a), pmu=1)
#' gnlr(y, dist="Poisson", mu=~exp(linear), linear=~x2, pmu=c(1,0.2))
#' # zero-inflated Poisson ZIP
#' fmr(y, dist="Poisson", mu=~exp(a), pmu=1, pmix=0)
#' fmr(y, dist="Poisson", mu=~exp(linear), linear=~x2, pmu=c(1,0.2), pmix=0)
#' fmr(y, dist="Poisson", mu=~exp(a), mix=~x1, pmu=1, pmix=c(1,0))
#' fmr(y, dist="Poisson", mu=~exp(linear), linear=~x2, mix=~x1, pmu=c(1,0.2),
#' 	pmix=c(1,0))
#' # zero-inflated negative binomial
#' fmr(y, dist="negative binomial", mu=~exp(a), pmu=1, pshape=0, pmix=0)
#' fmr(y, dist="negative binomial", mu=~exp(linear), linear=~x2, pmu=c(1,0.2),
#' 	pshape=0, pmix=0)
#' fmr(y, dist="negative binomial", mu=~exp(a), mix=~x1, pmu=1, pshape=0,
#'        pmix=c(1,0))
#' fmr(y, dist="negative binomial", mu=~exp(linear), linear=~x2, mix=~x1,
#' 	pmu=c(1,0.2), pshape=0, pmix=c(1,0))
#' 
#' @export fmr
fmr <- function(y=NULL, distribution="normal", mu=NULL, mix=NULL, linear=NULL,
	pmu=NULL, pmix=NULL, pshape=NULL, censor="right", exact=FALSE,
	wt=1, delta=1, common=FALSE, envir=parent.frame(), print.level=0,
	typsize=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
#
# inverse Gaussian cdf
#
pinvgauss <- function(y,m,s){
	t <- y/m
	v <- sqrt(y*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
#
# Laplace cdf
#
plaplace <- function(y){
	t <- exp(-abs(y))/2
	ifelse(y<0,t,1-t)}
#
# Levy cdf
#
plevy <- function(y, m, s)
	.C("plevy",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(1),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		#DUP=FALSE,
		PACKAGE="gnlm")$res

call <- sys.call()
#
# check distribution
#
if(is.function(distribution)){
	fcn <- distribution
	distribution <- "own"}
else distribution <- match.arg(distribution,c("binomial","beta binomial",
	"double binomial","mult binomial","Poisson","negative binomial",
	"double Poisson","mult Poisson","gamma count","Consul","geometric",
	"normal","inverse Gauss","logistic","exponential","gamma","Weibull",
	"extreme value","Pareto","Cauchy","Student t","Laplace","Levy"))
#
# check for parameters common to location and mixture functions
#
if(common){
	if(!is.function(mu)&&!inherits(mu,"formula"))
		stop("with common parameters, mu must be a function or formula")
	if(!is.function(mix)&&!inherits(mix,"formula"))
		stop("with common parameters, mix must be a function or formula")
	if(!is.null(linear))stop("linear cannot be used with common parameters")}
#
# count number of parameters
#
npl <- length(pmu)
npm <- length(pmix)
sht <- distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"
if(!sht)pshape <- NULL
if(sht&&is.null(pshape))
	stop("An estimate of the shape parameter must be given")
np <- npl+npm+sht
#
# find number of observations now for creating null functions
#
n <- if(inherits(envir,"repeated")||inherits(envir,"response"))sum(nobs(envir))
	else if(inherits(envir,"data.frame"))dim(envir)[1]
	else if(is.vector(y,mode="numeric"))length(y)
	else if(is.matrix(y))dim(y)[1]
	else sum(nobs(y))
if(n==0)stop(paste(deparse(substitute(y)),"not found or of incorrect type"))
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("fmr only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("fmr does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(inherits(envir,"repeated")||inherits(envir,"response"))
		deparse(substitute(envir))
	else NULL
#
# find linear part of each regression and save model for printing
#
lin1 <- lin2 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]}
else lin1 <- linear
if(inherits(lin1,"formula")&&is.null(mu)){
	mu <- lin1
	lin1 <- NULL}
if(inherits(lin2,"formula")&&is.null(mix)){
	mix <- lin2
	lin2 <- NULL}
if(inherits(lin1,"formula")){
	lin1model <- if(respenv){
		if(!is.null(attr(finterp(lin1,.envir=y,.name=envname),"parameters")))
			attr(finterp(lin1,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(lin1,.envir=envir,.name=envname),"parameters")))
			attr(finterp(lin1,.envir=envir,.name=envname),"model")}}
else lin1model <- NULL
if(inherits(lin2,"formula")){
	lin2model <- if(respenv){
		if(!is.null(attr(finterp(lin2,.envir=y,.name=envname),"parameters")))
			attr(finterp(lin2,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(lin2,.envir=envir,.name=envname),"parameters")))
			attr(finterp(lin2,.envir=envir,.name=envname),"model")}}
else lin2model <- NULL
#
# check if linear contains W&R formula
#
if(inherits(lin1,"formula")){
	tmp <- attributes(if(respenv)finterp(lin1,.envir=y,.name=envname)
		else finterp(lin1,.envir=envir,.name=envname))
	lf1 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	if(length(tmp$model)==1){
		if(is.null(mu))mu <- ~1
		else stop("linear must contain covariates")}
	rm(tmp)}
else lf1 <- 0
if(inherits(lin2,"formula")){
	tmp <- attributes(if(respenv)finterp(lin2,.envir=y,.name=envname)
		else finterp(lin2,.envir=envir,.name=envname))
	lf2 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	if(length(tmp$model)==1){
		if(is.null(mix))mix <- ~1
		else stop("linear must contain covariates")}
	rm(tmp)}
else lf2 <- 0
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu2 <- mixt2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	# modify formulae
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	if(inherits(mix,"formula")){
		mix2 <- if(respenv)finterp(mix,.envir=y,.name=envname)
			else finterp(mix,.envir=envir,.name=envname)}
	# modify functions
	if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
		        tmp <- parse(text=deparse(mu)[-1])
		        mu <- if(respenv)fnenvir(mu,.envir=y,.name=envname)
		        	else fnenvir(mu,.envir=envir,.name=envname)
		        mu2 <- mu
		        attr(mu2,"model") <- tmp}
		else mu2 <- mu}
	if(is.function(mix)){
		if(is.null(attr(mix,"model"))){
		        tmp <- parse(text=deparse(mix)[-1])
		        mix <- if(respenv)fnenvir(mix,.envir=y,.name=envname)
		        	else fnenvir(mix,.envir=envir,.name=envname)
		        mixt2 <- mix
		        attr(mixt2,"model") <- tmp}
		else mixt2 <- mix}}
else {
     if(is.function(mu)&&is.null(attr(mu,"model")))mu <- fnenvir(mu)
     if(is.function(mix)&&is.null(attr(mix,"model")))
		mix <- fnenvir(mix)}
#
# transform location formula to function and check number of parameters
#
if(inherits(mu,"formula")){
	if(npl==0)stop("formula for mu cannot be used if no parameters are estimated")
	linarg <- if(lf1>0) "linear" else NULL
	mu3 <- if(respenv)finterp(mu,.envir=y,.name=envname,.args=linarg)
		else finterp(mu,.envir=envir,.name=envname,.args=linarg)
	npt1 <- length(attr(mu3,"parameters"))
	if(is.character(attr(mu3,"model"))){
	# W&R formula
		if(length(attr(mu3,"model"))==1){
		# intercept model
			tmp <- attributes(mu3)
			mu3 <- function(p) p[1]*rep(1,n)
			attributes(mu3) <- tmp}}
	else {
	# formula with unknowns
		if(npl!=npt1&&!common&&lf1==0){
			cat("\nParameters are ")
			cat(attr(mu3,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu3,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}
	if(npl<npt1)stop("Not enough initial estimates for mu")}
else if(!is.function(mu)){
	mu3 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	mu3 <- mu
	npt1 <- length(attr(mu3,"parameters"))-(lf1>0)}
#
# if linear part, modify location function appropriately
#
if(lf1>0){
	if(is.character(attr(mu3,"model")))
		stop("mu cannot be a W&R formula if linear is supplied")
	dm1 <- if(respenv)wr(lin1,data=y)$design
		else wr(lin1,data=envir)$design
	if(is.null(mu2))mu2 <- mu3
	mu1 <- function(p)mu3(p,dm1%*%p[(npt1+1):(npt1+lf1)])}
else {
	if(lf1==0&&length(mu3(pmu))==1){
		mu1 <- function(p) mu3(p)*rep(1,n)
		attributes(mu1) <- attributes(mu3)}
	else {
		mu1 <- mu3
		rm(mu3)}}
#
# give appropriate attributes to mu1 for printing
#
if(is.null(attr(mu1,"parameters"))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,.envir=y))
			else attributes(fnenvir(mu,.envir=envir))}
		else attributes(mu)}
		else attributes(fnenvir(mu1))}
#
# check that correct number of estimates was supplied
#
nlp <- npt1+lf1
if(!common&&nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
npl1 <- if(common&&!inherits(mix,"formula")) 1 else npl+1
#
# transform mixture formula to function and check number of parameters
#
if(inherits(mix,"formula")){
	if(npm==0&&!common)
		stop("formula for mix cannot be used if no parameters are estimated")
	old <- if(common)mu1 else NULL
	linarg <- if(lf2>0) "linear" else NULL
	mixt4 <- if(respenv)finterp(mix,.envir=y,.start=npl1,.name=envname,.old=old,.args=linarg)
		else finterp(mix,.envir=envir,.start=npl1,.name=envname,.old=old,.args=linarg)
	npt2 <- length(attr(mixt4,"parameters"))
	if(is.character(attr(mixt4,"model"))){
	# W&R formula
		if(length(attr(mixt4,"model"))==1){
		# intercept model
			mixt3 <- function(p)
				1/(1+exp(-p[npl1]*rep(1,n)))
			mixt2 <- fnenvir(function(p)
				1/(1+exp(-p[1]*rep(1,n))))
			rm(mixt4)}
		else {
		# design matrix
			tmp <- attributes(mixt4)
			dm2 <- if(respenv)wr(mix,data=y)$design
				else wr(mix,data=envir)$design
			mixt3 <- function(p)
				1/(1+exp(-dm2%*%p[npl1:(npl1+npt2-1)]))
			attributes(mixt3) <- tmp
			rm(mixt4)}}
	else {
	# formula with unknowns
		if(npm!=npt2&&!common&&lf2==0){
			cat("\nParameters are ")
			cat(attr(mixt4,"parameters"),"\n")
			stop(paste("pmix should have",npt2,"estimates"))}
		mixt3 <- function(p)1/(1+exp(-mixt4(p)))
		attributes(mixt3) <- attributes(mixt4)
		if(is.list(pmix)){
			if(!is.null(names(pmix))){
				o <- match(attr(mixt3,"parameters"),names(pmix))
				pmix <- unlist(pmix)[o]
				if(sum(!is.na(o))!=length(pmix))stop("invalid estimates for mix - probably wrong names")}
			else pmix <- unlist(pmix)}}}
else if(!is.function(mix)){
	mixt3 <- function(p) exp(p[npl1])/(1+exp(p[npl1]))*rep(1,n)
	mixt2 <- fnenvir(function(p) exp(p[1])/(1+exp(p[1]))*rep(1,n))
	npt2 <- 1}
else {
	mixt3 <- function(p) 1/(1+exp(-mix(p[npl1:np])))
	attributes(mixt3) <- attributes(mix)
	npt2 <- length(attr(mixt3,"parameters"))-(lf2>0)}
#
# if linear part, modify mix function appropriately
#
if(lf2>0){
	if(is.character(attr(mixt3,"model")))
		stop("mix cannot be a W&R formula if linear is supplied")
	dm2 <- if(respenv)wr(lin2,data=y)$design
		else wr(lin2,data=envir)$design
	if(is.null(mixt2))mixt2 <- mixt3
	mixt <-mixt3(p,dm2%*%p[(npl1+lf2-1):np])}
else {
	mixt <- mixt3
	rm(mixt3)}
#
# give appropriate attributes to mixt for printing
#
if(is.null(attr(mixt,"parameters"))){
	attributes(mixt) <- if(is.function(mix)){
		if(!inherits(mix,"formulafn")){
			if(respenv)attributes(fnenvir(mix,.envir=y))
			else attributes(fnenvir(mix,.envir=envir))}
		else attributes(mix)}
		else attributes(fnenvir(mixt))}
#
# check that correct number of estimates was supplied
#
nlp <- npt2+lf2
if(!common&&nlp!=npm)stop(paste("pmix should have",nlp,"initial estimates"))
#
# when there are parameters common to location and shape functions,
# check that correct number of estimates was supplied
#
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(mixt,"parameters"))))
	if(nlp!=npl)stop(paste("with a common parameter model, pmu should contain",nlp,"estimates"))}
p <- c(pmu,pmix,pshape)
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
	respname <- colnames(y$response$y)
	y <- response(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("fmr does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	respname <- colnames(envir$response$y)[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$n)&&!all(is.na(envir$response$n[,col])))
		y <- cbind(y,envir$response$n[,col]-y)
	else if(!is.null(envir$response$censor)&&!all(is.na(envir$response$censor[,col])))
		y <- cbind(y,envir$response$censor[,col])
	if(!is.null(envir$response$wt))wt <- as.vector(envir$response$wt)
	if(!is.null(envir$response$delta))
		delta <- as.vector(envir$response$delta[,col])}
else if(inherits(envir,"data.frame")){
	respname <- deparse(substitute(y))
	y <- envir[[deparse(substitute(y))]]}
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("fmr only handles univariate responses")
	if(!is.null(y$wt)&&any(!is.na(y$wt)))wt <- as.vector(y$wt)
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	type <- y$type
	respname <- colnames(y$y)
	y <- response(y)}
else respname <- deparse(substitute(y))
if(any(is.na(y)))stop("NAs in y - use rmna")
#
# set up censoring indicators for likelihood
#
if(distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="double Poisson"||distribution=="mult Poisson"||
	distribution=="gamma count"||distribution=="Consul"){
	# count data
	if(type!="unknown"&&type!="discrete")stop("discrete data required")
	if(!is.vector(y,mode="numeric"))stop("y must be a vector")
	censor <- NULL
	cens <- ifelse(y==0,1,0)}
else {
	if(distribution=="binomial"||distribution=="double binomial"||
		distribution=="beta binomial"||distribution=="mult binomial"){
		# binomial data
		if(type!="unknown"&&type!="nominal")
			stop("nominal data required")
		if(distribution=="binomial"&&(is.vector(y)||(length(dim(y))==2
			&&dim(y)[2]==1))&&all(y==0|y==1))y <- cbind(y,1-y)
		if(any(y<0))stop("All response values must be positive")}
	if(length(dim(y))!=2||dim(y)[2]!=2)
		stop(paste("Two column matrix required for response:",
		if(distribution=="binomial"||distribution=="beta binomial"||
			distribution=="double binomial"||
			distribution=="mult binomial")"successes and failures"
		else "times and censor indicator"))
	else {
		if(distribution=="binomial"||distribution=="beta binomial"||
			distribution=="double binomial"||
			distribution=="mult binomial"){
		# binomial data
			if(is.null(censor))
				stop("Censoring must be left, right, or both")
			if(censor!="left"&&censor!="right"&&censor!="both")
				stop("Censoring must be left, right, or both")
			lcens <- if((censor=="left"|censor=="both")&y[,1]==0)1
				else 0
			rcens <- if((censor=="right"|censor=="both")&y[,2]==0)1
				else 0
			if(censor=="both"){
				lcens <- lcens/2
				rcens <- rcens/2}
			nn <- y[,1]+y[,2]}
		else {
		# continuous data
			if(any(delta<=0&y[,2]==1))
				stop("All deltas for uncensored data must be positive")
			else {
				delta <- ifelse(delta<=0,0.000001,delta)
				delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001
				,delta)}
			y[,2] <- as.integer(y[,2])
			if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
				stop("Censor indicator must be -1, 0, or 1")
			if(censor!="left"&&censor!="right")
				stop("Censoring must be left or right")
			if(censor=="left"&&!any(y[,2]==-1))
				stop("No left censored observations")
			if(censor=="right"&&!any(y[,2]==0))
				stop("No right censored observations")
			cens <- as.integer(y[,2]==1)
			b <- as.integer((censor=="right"&y[,2]==0)|
				(censor=="left"&y[,2]==-1))
			r <- as.integer(censor=="left"&y[,2]==0)
			l <- as.integer(censor=="right"&y[,2]==-1)
			lc <- if(censor=="left")1 else 0
			rc <- if(censor=="right")-1 else 1}}
	if(distribution=="double Poisson"||distribution=="mult Poisson")
				my <- min(3*max(y),100)}
#
# check that data are appropriate for distribution
#
if(distribution!="normal"&&distribution!="logistic"&&distribution!="Cauchy"&&
	distribution!="Laplace"&&distribution!="Student t"&&
	distribution!="Poisson"&&distribution!="negative binomial"&&
	distribution!="Consul"&&distribution!="double Poisson"&&
	distribution!="mult Poisson"&&distribution!="gamma count"&&
	distribution!="binomial"&& distribution!="beta binomial"&&
	distribution!="double binomial"&&distribution!="mult binomial"){
	if(type!="unknown"&&type!="duration"&&type!="continuous")
		stop("duration data required")
	if(any(y[,1]<=0))stop("All response values must be > 0")}
else if(distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="gamma count"||distribution=="double Poisson"||
	distribution=="mult Poisson"||distribution=="Consul"||
	distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial"){
	if(type!="unknown"&&type!="discrete")stop("discrete data required")
	if(any(y<0))stop("All response values must be >= 0")}
else if(type!="unknown"&&type!="continuous"&&type!="duration")
	stop("continuous data required")
#
# prepare weights and unit of measurement
#
if(min(wt)<0)stop("All weights must be non-negative")
if(length(wt)==1)wt <- rep(wt,n)
if(length(delta)==1)delta <- rep(delta,n)
#
# check that location function returns appropriate values
#
if(any(is.na(mu1(pmu))))stop("The location regression returns NAs: probably invalid initial values")
if(distribution=="Levy"&&any(y[,1]<=mu1(p)))
	stop("location parameter must be strictly less than corresponding observation")
#
# check that mixture function returns appropriate values
#
if(any(is.na((mixt(p)))))
	stop("The mix function returns NAs: probably invalid initial values")
#
# create the appropriate likelihood function
s <- NULL
ret <- switch(distribution,
	binomial={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			-sum(wt*log((1-s)*(lcens+rcens)+s*dbinom(y[,1],
				y[,1]+y[,2],m)))}
		const <- 0},
	"beta binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			v <- exp(p[np])
			t <- v*m
			u <- v*(1-m)
			-sum(wt*log((1-s)*(lcens+rcens)+s*
				exp(lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u)+
				lchoose(nn,y[,1]))))}
		const <- 0},
	"double binomial"={
		fcn <- function(p) {
			-sum(wt*log((1-s)*(lcens+rcens)+s*exp(.C("ddb",
				as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(p[np])),
				as.integer(n),as.double(wt),
				res=double(n),#DUP=FALSE,
				PACKAGE="gnlm")$res)))}
		const <- 0},
	"mult binomial"={
		fcn <- function(p) {
			-sum(wt*log((1-s)*(lcens+rcens)+s*exp(.C("dmb",
				as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(p[np])),
				as.integer(n),as.double(wt),
				res=double(n),#DUP=FALSE,
				PACKAGE="gnlm")$res)))}
		const <- 0},
	Poisson={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			-sum(wt*log((1-s)*cens+s*dpois(y,m)))}
		const <- 0},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			t <- exp(p[np])
			-sum(wt*log((1-s)*cens+s*dnbinom(y,t,mu=m)))}
		const <- 0},
	"double Poisson"={
		fcn <- function(p) {
			-sum(wt*log((1-s)*cens+s*exp(.C("ddp",as.integer(y),
				as.integer(my),as.double(mu1(p)),
				as.double(exp(p[np])),as.integer(length(y)),
				as.double(wt),res=double(length(y)),
				#DUP=FALSE,
				PACKAGE="gnlm")$res)))}
		const <- 0},
	"mult Poisson"={
		fcn <- function(p) {
			-sum(wt*log((1-s)*cens+s*exp(.C("dmp",as.integer(y),
				as.integer(my),as.double(mu1(p)),
				as.double(exp(p[np])),as.integer(length(y)),
				as.double(wt),res=double(length(y)),
				#DUP=FALSE,
				PACKAGE="gnlm")$res)))}
		const <- 0},
	"gamma count"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			t <- exp(p[np])
			-sum(wt*log((1-s)*cens+s*ifelse(y==0,1-pgamma(m*t,
				(y+1)*t),pgamma(m*t,y*t+(y==0))-
				pgamma(m*t,(y+1)*t))))}
		const <- 0},
	Consul={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			u <- exp(p[np])
			-sum(wt*log((1-s)*cens+s*exp(log(m)-(m+y*(u-1))/u
				-y*p[np]+(y-1)*log(m+y*(u-1))-lgamma(y+1))))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pn <- pnorm(y[,1],m,t)
				-sum(wt*log(s*cens*(pnorm(y[,1]+delta/2,m,t)-
					pnorm(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pn-lc))*b
					+s*(r+pn*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pn <- pnorm(y[,1],m,t)
				-sum(wt*log(s*cens*dnorm(y[,1],m,t)
					+(1-cens)*((1+s*(rc*pn-lc))*b+
					s*(r+pn*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pit <- pinvgauss(y[,1],m,t)
				-sum(wt*log(s*cens*(pinvgauss(y[,1]+delta/2,
					m,t)-pinvgauss(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pit-lc))*b
					+s*(r+pit*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pit <- pinvgauss(y[,1],m,t)
				-sum(wt*log(s*cens*exp(-(p[np]+(y[,1]-m)^2/
					(y[,1]*t*m^2))/2)
					+(1-cens)*((1+s*(rc*pit-lc))*b
					+s*(r+pit*(l-r)))))}
			const <- wt*cens*(log(2*pi*y[,1]^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])*sqrt(3)/pi
				pl <- plogis(y[,1],m,t)
				-sum(wt*log(s*cens*(plogis(y[,1]+delta/2,m,t)-
					plogis(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])*sqrt(3)/pi
				y1 <- (y[,1]-m)/t
				pl <- plogis(y[,1],m,t)
				-sum(wt*log(s*cens*dlogis(y[,1],m,t)
					+(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ps <- pt(y[,1]-m,t)
				-sum(wt*log(s*cens*(pt(y[,1]+delta/2-m,t)-
					pt(y[,1]-delta/2-m,t))
					+(1-cens)*((1+s*(rc*ps-lc))*b
					+s*(r+ps*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ps <- pt(y[,1]-m,t)
				-sum(wt*log(s*cens*dt(y[,1]-m,t)
					+(1-cens)*((1+s*(rc*ps-lc))*b
					+s*(r+ps*(l-r)))))}
			const <- -wt*cens*log(delta)}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pc <- pcauchy(y[,1],m,t)
				-sum(wt*log(s*cens*(pcauchy(y[,1]+delta/2,m,t)-
					pcauchy(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pc-lc))*b
					+s*(r+pc*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pc <- pcauchy(y[,1],m,t)
				-sum(wt*log(s*cens*dcauchy(y[,1],m,t)
					+(1-cens)*((1+s*(rc*pc-lc))*b
					+s*(r+pc*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plaplace((y[,1]-m)/t)
				-sum(wt*log(s*cens*(plaplace((y[,1]+delta/2-
					m)/t)-plaplace((y[,1]-delta/2-m)/t))+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plaplace((y[,1]-m)/t)
				-sum(wt*log(s*cens*exp(-abs(y[,1]-m)/t-p[np])+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r)))))}
			const <- -wt*cens*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plevy(y[,1],m,t)
				-sum(wt*log(s*cens*(plevy(y[,1]+delta/2,m,t)
					-plevy(y[,1]-delta/2,m,t))+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plevy(y[,1],m,t)
				-sum(wt*log(s*cens*sqrt(t/(2*pi))*log(y[,1]-m)^
					-1.5*exp(-t/(2*(y[,1]-m)))+(1-cens)*
					((1+s*(rc*pl-lc))*b+s*(r+pl*(l-r)))))}
			const <- -wt*cens*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- mixt(p)
				u <- exp(p[np])
				t <- 1/(mu1(p)*u)
				pp <- 1-(1+y[,1]*t)^-u
				-sum(wt*log(s*cens*((1+(y[,1]-delta/2)*t)^-u-
					(1+(y[,1]+delta/2)*t)^-u)
					+(1-cens)*((1+s*(rc*pp-lc))*b
					+s*(r+pp*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- mixt(p)
				u <- exp(p[np])
				t <- 1/(mu1(p)*u)
				pp <- 1-(1+y[,1]*t)^-u
				-sum(wt*log(s*cens*u*t*(1+y[,1]*t)^(-(u+1))+
					(1-cens)*
					((1+s*(rc*pp-lc))*b+s*(r+pp*(l-r)))))}
			const <- -wt*cens*log(delta)}},
	exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				u <- exp(-y[,1]/m)
				-sum(wt*log(s*cens*(-exp(-(y[,1]+delta/2)/m)+
					exp(-(y[,1]-delta/2)/m))
					+(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				u <- exp(-y[,1]/m)
				-sum(wt*log(s*cens*exp(-y[,1]/m)/m
					+(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				u <- m/t
				pg <- pgamma(y[,1],t,scale=u)
				-sum(wt*log(s*cens*(pgamma(y[,1]+delta/2,t,
					scale=u)-pgamma(y[,1]-delta/2,t,
					scale=u))+(1-cens)*((1+s*(rc*pg-lc))*b
					+s*(r+pg*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				u <- m/t
				pg <- pgamma(y[,1],t,scale=u)
				-sum(wt*log(s*cens*dgamma(y[,1],t,scale=u)
					+(1-cens)*((1+s*(rc*pg-lc))*b
					+s*(r+pg*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pw <- pweibull(y[,1],t,m)
				-sum(wt*log(s*cens*(pweibull(y[,1]+delta/2,t,m)
					-pweibull(y[,1]-delta/2,t,m))
					+(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pw <- pweibull(y[,1],t,m)
				-sum(wt*log(s*cens*dweibull(y[,1],t,m)+
					(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r)))))}
			const <- -wt*cens*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ey <- exp(y[,1])
				pw <- pweibull(ey,t,m)
				-sum(wt*log(s*cens*(pweibull(ey+ey*delta/2,
					t,m)-pweibull(ey-ey*delta/2,t,m))+
					(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r)))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ey <- exp(y[,1])
				pw <- pweibull(ey,t,m)
				-sum(wt*log(s*cens*dweibull(ey,t,m)+
					(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r)))))}
			const <- -wt*cens*log(delta)}},
	own={const <- 0})
#
# check that the likelihood returns an appropriate value and optimize
#
if(fscale==1)fscale <- fcn(p)
if(is.na(fcn(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fcn, p=p, hessian=TRUE, print.level=print.level, typsize=typsize,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
#
# calculate fitted values and raw residuals
#
fitted.values <- if(distribution=="binomial"||distribution=="beta binomial"||
		distribution=="double binomial"||distribution=="mult binomial")
		as.vector((y[,1]+y[,2])*mu1(z0$estimate))
	else as.vector(mu1(z0$estimate))
residuals <- if(distribution!="Poisson"&&distribution!="negative binomial"&&
	distribution!="Consul"&&distribution!="double Poisson"&&
	distribution!="mult Poisson"&&distribution!="gamma count")
		y[,1]-fitted.values
	else y-fitted.values
#
# calculate se's
#
if(np==1)cov <- 1/z0$hessian
else {
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf))0
		else qr(z0$hessian)$rank
	if(a==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)}
se <- sqrt(diag(cov))
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
if(!is.null(mixt2))mixt <- mixt2
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	respname=respname,
	mu=mu1,
	mix=mixt,
	linear=list(lin1,lin2),
	linmodel=list(lin1model,lin2model),
	common=common,
	prior.weights=wt,
	censor=censor,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=z0$minimum+np,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	npl=npl,
	npm=npm,
	nps=as.numeric(sht),
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlm"
return(z1)}
