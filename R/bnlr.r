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
#     bnlr(y=NULL, link="logit", mu=NULL, linear=NULL, pmu=NULL, pshape=NULL,
#	wt=1, envir=parent.frame(), print.level=0, typsize=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
#	iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit binomial nonlinear regression models with a variety of
# link functions.



#' Binomial Nonlinear Regression Models
#' 
#' \code{bnlr} fits user-specified nonlinear regression equations to binomial
#' data with various link functions (\code{logit}, \code{probit}, \code{comp
#' log log}, \code{log log}, \code{Cauchy}, \code{Student t}, \code{stable}, or
#' \code{mixture}). The mixture link is a logistic link with extra probability
#' mass for \code{y=0} and \code{y=n}.
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
#' @param y A two column matrix of binomial data or censored data or an object
#' of class, \code{response} (created by \code{\link[rmutil]{restovec}}) or
#' \code{repeated} (created by \code{\link[rmutil]{rmna}} or
#' \code{\link[rmutil]{lvna}}). If the \code{repeated} data object contains
#' more than one response variable, give that object in \code{envir} and give
#' the name of the response variable to be used here.
#' @param link A character string containing the name of the link function. The
#' \code{Student t}, \code{stable}, and \code{mixture} links contain an unknown
#' parameter to be estimated, respectively the logarithm of the degrees of
#' freedom, the tail parameter transformed by log(tail/(2-tail)), and logit of
#' the mixture probability, so that they lie on the whole real line.
#' @param mu A user-specified function of \code{pmu}, and possibly
#' \code{linear}, giving the regression equation for the location. This may
#' contain a linear part as the second argument to the function. It may also be
#' a formula beginning with ~, specifying either a linear regression function
#' for the location parameter in the Wilkinson and Rogers notation or a general
#' function with named unknown parameters. If it contains unknown parameters,
#' the keyword \code{linear} may be used to specify a linear part. If nothing
#' is supplied, the location is taken to be constant unless the linear argument
#' is given.
#' @param linear A formula beginning with ~ in W&R notation, specifying the
#' linear part of the regression function for the location parameter or list of
#' two such expressions for the location and/or shape parameters.
#' @param pmu Vector of initial estimates for the location parameters. If
#' \code{mu} is a formula with unknown parameters, their estimates must be
#' supplied either in their order of appearance in the expression or in a named
#' list.
#' @param pshape If the \code{link} is \code{Student t}, an initial estimate of
#' the degrees of freedom; if it is \code{stable}, an estimate of the tail
#' parameter; if it is \code{mixture}, an estimate of the mixture probability.
#' @param wt Weight vector.
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
#' \code{\link[gnlm]{gnlr}}, \code{\link[gnlm]{gnlr3}}
#' @keywords models
#' @examples
#' 
#' # assay to estimate LD50
#' y <- c(9,9,10,4,1,0,0)
#' y <- cbind(y,10-y)
#' dose <- log10(100/c(2.686,2.020,1.520,1.143,0.860,0.647,0.486))
#' 
#' summary(glm(y~dose, family=binomial))
#' bnlr(y, mu=~dose, pmu=c(1,1))
#' summary(glm(y~dose, family=binomial(link=probit)))
#' bnlr(y, link="probit", mu=~dose, pmu=c(1,1))
#' \dontrun{
#' bnlr(y, link="log log", mu=~dose, pmu=c(1,1))
#' bnlr(y, link="comp log log", mu=~dose, pmu=c(1,1))
#' bnlr(y, link="Cauchy", mu=~dose, pmu=c(60,-30))
#' bnlr(y, link="Student", mu=~dose, pmu=c(60,-30), pshape=0.1)
#' bnlr(y, link="stable", mu=~dose, pmu=c(20,-15), pshape=0, stepmax=1)
#' bnlr(y, link="mixture", mu=~dose, pmu=c(60,-30), pshape=-2.5)
#' #
#' mu <- function(p) -p[1]*(log10(p[2])-dose)
#' bnlr(y, mu=mu, pmu=c(1,100))
#' bnlr(y, link="probit", mu=mu, pmu=c(1,100))
#' }
#' @export bnlr
bnlr <- function(y=NULL, link="logit", mu=NULL, linear=NULL, pmu=NULL,
	pshape=NULL, wt=1, envir=parent.frame(), print.level=0,
	typsize=abs(p),ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){
#
# stable cdf
#
pstable <- function(y,tail){
	z <- .C("pstable",
		as.integer(length(y)),
		y=y,
		skew=rep(0,n),
		tail=rep(tail,n), 
		eps=1.0e-6,
		err=integer(1),
		ffy=double(length(y)),
		PACKAGE="gnlm")
	z$ffy}
#
call <- sys.call()
#
# check link
#
link <- match.arg(link,c("logit","probit","comp log log","log log",
	"Cauchy","Student t","stable","mixture"))
#
# count number of parameters
#
if(link=="Student t"||link=="stable"||link=="mixture"){
	if(is.null(pshape)||length(pshape)!=1)stop("pshape must be a scalar")}
else pshape <- NULL
npl <- length(pmu)
#
# find number of observations now for creating null functions
#
n <- if(inherits(envir,"repeated")||inherits(envir,"response"))sum(nobs(envir))
	else if(inherits(envir,"data.frame"))dim(envir)[1]
	else if(is.vector(y,mode="numeric"))
		stop("y must be a two-column matrix")
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
		stop("bnlr only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("bnlr does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(inherits(envir,"repeated")||inherits(envir,"response"))
		deparse(substitute(envir))
	else NULL
#
# find linear part of regression and save model for printing
#
if(inherits(linear,"formula")&&is.null(mu)){
	mu <- linear
	linear <- NULL}
if(inherits(linear,"formula")){
	lin1model <- if(respenv){
		if(!is.null(attr(finterp(linear,.envir=y,.name=envname),"parameters")))
			attr(finterp(linear,.envir=y,.name=envname),"model")}
	else {if(!is.null(attr(finterp(linear,.envir=envir,.name=envname),"parameters")))
			attr(finterp(linear,.envir=envir,.name=envname),"model")}}
else lin1model <- NULL
#
# if a data object was supplied, modify formulae or functions to read from it
#
mu2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	# modify formulae
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
	# modify functions
	if(is.function(mu)){
		if(is.null(attr(mu,"model"))){
			tmp <- parse(text=deparse(mu)[-1])
			mu <- if(respenv)fnenvir(mu,.envir=y,.name=envname)
				else fnenvir(mu,.envir=envir,.name=envname)
			mu2 <- mu
			attr(mu2,"model") <- tmp}
		else mu2 <- mu}}
else if(is.function(mu)&&is.null(attr(mu,"model")))mu <- fnenvir(mu)
#
# check if linear contains W&R formula
#
if(inherits(linear,"formula")){
	tmp <- attributes(if(respenv)finterp(linear,.envir=y,.name=envname)
		else finterp(linear,.envir=envir,.name=envname))
	lf1 <- length(tmp$parameters)
	if(!is.character(tmp$model))stop("linear must be a W&R formula")
	else if(length(tmp$model)==1)stop("linear must contain covariates")
	rm(tmp)}
else lf1 <- 0
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
			tmp <- attributes(mu1)
			mu1 <- function(p) p[1]*rep(1,n)
			attributes(mu1) <- tmp}}
	else {
	# formula with unknowns
		if(npl!=npt1&&lf1==0){
			cat("\nParameters are ")
			cat(attr(mu3,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu3,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}}
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
	dm1 <- if(respenv)wr(linear,data=y)$design
		else wr(linear,data=envir)$design
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
if(nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
if(link=="Student t"||link=="stable"||link=="mixture")
	sh1 <- fnenvir(function(p) p[1]*rep(1,n))
else sh1 <- NULL
p <- c(pmu,pshape)
np <- length(p)
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(is.null(y$response$n)||all(is.na(y$response$n)))
		stop("these are not binomial data")
	type <- y$response$type
	respname <- colnames(y$response$y)
	y <- cbind(response(y),y$response$n)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("bnlr does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	respname <- colnames(envir$response$y)[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$n)&&!all(is.na(envir$response$n[,col])))
		y <- cbind(y,envir$response$n[,col]-y)}
else if(inherits(envir,"data.frame")){
	respname <- deparse(substitute(y))
	y <- envir[[deparse(substitute(y))]]}
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)stop("bnlr only handles univariate responses")
	if(is.null(y$n)||all(is.na(y$n)))stop("these are not binomial data")
	type <- y$type
	respname <- colnames(y$y)
	y <- cbind(response(y),y$n)}
else respname <- deparse(substitute(y))
if(any(is.na(y)))stop("NAs in y - use rmna")
#
# check that data are appropriate for distribution
#
if(length(dim(y))!=2||dim(y)[2]!=2)
	stop(paste("Two column matrix required for response: successes and failures"))
if(type!="unknown"&&type!="nominal")stop("nominal data required")
if(any(y<0))stop("All response values must be positive")
#
# prepare weights and unit of measurement
#
if(length(wt)==1)wt <- rep(wt,n)
else if(length(wt)!=n)stop("wt must be the same length as the other variables")
if(min(wt)<0)stop("All weights must be non-negative")
nn <- y[,1]+y[,2]
#
# check that location function returns appropriate values
#
if(any(is.na(mu1(pmu))))stop("The location regression returns NAs: probably invalid initial values")
#
# create the appropriate likelihood function
#
	ret <- switch(link,
	        logit=fcn <- function(p) {
	        	m <- plogis(mu1(p))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        probit=fcn <- function(p) {
	        	m <- pnorm(mu1(p))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        "comp log log"=fcn <- function(p) {
	        	m <- 1-exp(-exp(mu1(p)))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        "log log"=fcn <- function(p) {
	        	m <- exp(-exp(mu1(p)))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        Cauchy=fcn <- function(p) {
	        	m <- pcauchy(mu1(p))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        "Student t"=fcn <- function(p) {
	        	m <- pt(mu1(p),exp(p[np]))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        "stable"=fcn <- function(p) {
	        	m <- pstable(mu1(p),2/(1+exp(-p[np])))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))},
	        "mixture"=fcn <- function(p){
	        	pi <- exp(p[np])/(1+exp(p[np]))
	        	m <- pi*(y[,1]==0||y[,2]==0)+(1-pi)*plogis(mu1(p))
	        	-sum(wt*(y[,1]*log(m)+y[,2]*log(1-m)))})
#
# check that the likelihood returns an appropriate value and optimize
#
if(fscale==1)fscale <- fcn(p)
if(is.na(fcn(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
if(np>0){
	z0 <- nlm(fcn,p=p,hessian=TRUE,print.level=print.level,typsize=typsize,
		ndigit=ndigit,gradtol=gradtol,stepmax=stepmax,steptol=steptol,
		iterlim=iterlim,fscale=fscale)
	z0$minimum <- z0$minimum-sum(wt*lchoose(nn,y[,1]))}
else z0 <- list(minimum=fscale-sum(wt*lchoose(nn,y[,1])),estimate=p,code=0,iterations=0)
#
# calculate fitted values and raw residuals
#
fitted.values <- as.vector((y[,1]+y[,2])*mu1(z0$estimate))
residuals <- y[,1]-fitted.values
#
# calculate se's
#
if(np==0)cov <- NULL
else if(np==1)cov <- 1/z0$hessian
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
z1 <- list(
	call=call,
	distribution="binomial",
	link=link,
	likefn=fcn,
	respname=respname,
	mu=mu1,
	shape=sh1,
	linear=list(linear,NULL),
	linmodel=list(lin1model,NULL),
	censor=FALSE,
	common=FALSE,
	prior.weights=wt,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=z0$minimum+np,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	npl=npl,
	npm=0,
	nps=np-npl,
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlm"
return(z1)}
