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
#     nordr(y=NULL, distribution="proportional", mu=NULL, linear=NULL,
#	pmu=NULL, pintercept=NULL, weights=NULL, envir=parent.frame(),
#	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1,
#	iterlim=100, typsize=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models for ordinal responses.

nordr <- function(y=NULL, distribution="proportional", mu=NULL, linear=NULL,
	pmu=NULL, pintercept=NULL, weights=NULL, envir=parent.frame(),
	print.level=0, ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1,
	iterlim=100, typsize=abs(p), stepmax=10*sqrt(p%*%p)){
#
# likelihood function for proportional odds and continuation ratio
#
lf <- function(p){
	g <- exp(mu1(p[1:npl])+block%*%p[npl1:np])
	g <- g/(1+g)
	if(mdl==1){
		g <- c(g,ext)
		g <- g[1:nlen]/g[nrows1:nlenr]
		g <- ifelse(g>=1,0.99,g)}
	else g <- 1-g
	-sum(pwt*(resp*log(g)+(1-resp)*log(1-g)))}
#
# likelihood function for adjacent categories
#
lf3 <- function(p){
	mu <- -mu1(p[1:npl])
	g <- exp(mu*y-resp%*%p[npl1:np])/
	exp(mu%o%(0:my)-matrix(rep(cumsum(c(0,0,p[npl1:np])),nrows),ncol=my+1,byrow=TRUE))%*%ext
	-sum(pwt*log(g))}
#
call <- sys.call()
#
# find number of observations now for creating null functions
#
nrows <- if(inherits(envir,"repeated")||inherits(envir,"response"))
		sum(nobs(envir))
	else if(is.vector(y,mode="numeric"))length(y)
	else if(is.matrix(y))stop("y must be a vector")
	else sum(nobs(y))
if(nrows==0)
	stop(paste(deparse(substitute(y)),"not found or of incorrect type"))
#
# check if a data object is being supplied
#
respenv <- exists(deparse(substitute(y)),env=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("nordr only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("nordr does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(inherits(envir,"repeated")||inherits(envir,"response"))
		deparse(substitute(envir))
	else NULL
#
# check model
#
tmp <- c("proportional odds","continuation ratio","adjacent categories")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
npl <- length(pmu)
npl1 <- npl+1
#
# find linear part of each regression and save model for printing
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
mu2 <- name <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")||inherits(envir,"data.frame")){
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,.envir=y,.name=envname)
			else finterp(mu,.envir=envir,.name=envname)}
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
	linarg <- if(lf1>0) "linear" else NULL
	mu3 <- if(respenv)finterp(mu,.envir=y,.name=envname,.args=linarg)
		else finterp(mu,.envir=envir,.name=envname,.args=linarg)
	npt1 <- length(attr(mu3,"parameters"))
	if(is.character(attr(mu3,"model"))){
	# W&R formula
		if(length(attr(mu3,"model"))==1){
		# intercept model
			tmp <- attributes(mu3)
			mu3 <- function(p) p[1]*rep(1,nrows)
			attributes(mu3) <- tmp}}
	else {
	# formula with unknowns
		if(npl!=npt1){
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
	mu3 <- function(p) p[1]*rep(1,nrows)
	npt1 <- 1}
else {
	mu3 <- mu
	if(is.matrix(mu3(pmu)))stop("mu must return a vector")
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
		mu1 <- function(p) mu3(p)*rep(1,nrows)
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
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&any(!is.na(y$response$wt)))
		weights <- as.vector(y$response$wt)
	type <- y$response$type
	y <- response(y)
	if(is.matrix(y))stop("response is not ordinal")}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("nordr does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$wt))weights <- as.vector(envir$response$wt)}
else if(inherits(envir,"data.frame"))y <- envir[[deparse(substitute(y))]]
else if(inherits(y,"response")){
	if(dim(y$y)[2]>1)
		stop("nordr only handles univariate responses")
	if(!is.null(y$wt)&&any(!is.na(y$wt)))weights <- as.vector(y$wt)
	type <- y$type
	y <- response(y)
	if(is.matrix(y))stop("response is not ordinal")}
if(any(is.na(y)))stop("NAs in y - use rmna")
if(!is.vector(y,mode="numeric")||any(y<0))
	stop("y must be a numeric vector with integral values starting at 0")
K <- length(unique(y))
if(min(y)!=0||max(y)!=K-1)
	stop(paste("ordinal values must be numbered from 0 to ",K-1))
else if(any(y!=trunc(y)))stop("ordinal values must be integers")
else my <- max(y)
if(type!="unknown"&&type!="ordinal")stop("ordinal data required")
#
# set up constants
#
nrows1 <- nrows+1
nlen <- my*nrows
nlenr <- nlen+nrows
if(any(is.na(mu1(pmu))))
	stop("The location regression returns NAs: probably invalid initial values")
if(missing(pintercept)||length(pintercept)!=my-1)
	stop(paste(my-1,"initial values of intercept contrasts must be supplied"))
tmp <- c(pintercept[1],diff(pintercept))
if(!(all(tmp>0)||all(tmp<0)))print("intercepts are not monotone")
#
# set up values for likelihood function
#
if(mdl==1)ext <- rep(1,nrows)
else if(mdl==3)ext <- rep(1,my+1)
if(mdl==3)resp <- NULL
else resp <- matrix(as.integer(y==0),ncol=1)
block <- NULL
pwt <- matrix(as.integer(y<2),ncol=1,nrow=nrows)
for(i in 2:my){
	resp <- cbind(resp,as.integer(y<i))
	block <- cbind(block,as.integer(c(rep(0,nrows*(i-1)),
		rep(1,nrows),rep(0,nrows*(my-i)))))
	pwt <- cbind(pwt,as.integer(y<i+1))}
if(mdl!=1)resp <- 1-resp
if(mdl!=3){
	resp <- as.vector(resp)
	pwt <- as.vector(pwt)}
else pwt <- rep(1,length(y))
#
# set up weights
#
if(!is.null(weights)){
	if(!is.vector(weights,mode="numeric"))stop("weights must be a vector")
	else if(length(weights)!=nrows)
		stop(paste("weights must have length",nrows))
	if(mdl==3)pwt <- weights
	else pwt <- rep(weights,my)*pwt}
#
# check that the likelihood returns an appropriate value and optimize
#
p <- c(pmu,pintercept)
np <- length(p)
if(mdl==3){
	if(fscale==1)fscale <- lf3(p)
	if(is.na(lf3(p)))
		stop("Likelihood returns NAs: probably invalid initial values")
	z <- nlm(lf3, p, hessian=TRUE, print.level=print.level,
		typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
		steptol=steptol, iterlim=iterlim, fscale=fscale)}
else {
	if(fscale==1)fscale <- lf(p)
	if(is.na(lf(p)))
		stop("Likelihood returns NAs: probably invalid initial values")
	z <- nlm(lf, p, hessian=TRUE, print.level=print.level,
		typsize=typsize, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
		steptol=steptol, iterlim=iterlim, fscale=fscale)}
maxlike <- z$minimum
#
# calculate se's
#
a <- if(any(is.na(z$hessian))||any(abs(z$hessian)==Inf))0
	else qr(z$hessian)$rank
if(a==np)cov <- solve(z$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
#
# return appropriate attributes on functions
#
if(!is.null(mu2))mu1 <- mu2
z1 <- list(
   call=call,
   distribution=distribution,
   weights=weights,
   maxlike=maxlike,
   aic=maxlike+np,
   mu=mu1,
   linear=linear,
   linmodel=lin1model,
   coefficients=z$estimate[1:npl],
   np=np,
   npl=npl1-1,
   nrows=nrows,
   intercept=z$estimate[npl1:np],
   cov=cov,
   corr=corr,
   se=se,
   iterations=z$iter,
   code=z$code)
class(z1) <- "nordr"
z1}

### print method
###
print.nordr <- function(z,digits=max(3,.Options$digits-3),correlation=TRUE){
m <- z$states
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
cat(z$distribution,"model\n\n")
if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
cat("-Log likelihood   ",z$maxlike,"\n")
cat("AIC               ",z$aic,"\n")
cat("Iterations        ",z$iterations,"\n")
cat("\nLocation coefficients\n")
if(inherits(z$mu,"formulafn")){
	cat("Location function:\n")
	if(!is.null(attr(z$mu,"formula")))
		cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	if(!is.null(z$linear)){
		cat("Linear part:\n")
		print(z$linear)}}
cname <- if(is.character(attr(z$mu,"model")))attr(z$mu,"model")
	else if(length(grep("linear",attr(z$mu,"parameters")))>0)
	attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
	else attr(z$mu,"parameters")
if(!is.null(z$linmodel))cname <- c(linmodel,cname)
coef.table <- cbind(z$coef,z$se[1:z$npl])
dimnames(coef.table) <- list(cname,c("estimate","s.e."))
print.default(coef.table, digits=digits, print.gap=2)
cat("\nIntercept contrasts\n")
coef.table <- cbind(z$intercept,z$se[(z$npl+1):z$np])
dimnames(coef.table) <- list(paste("b[",2:(z$np-z$npl+1),"]",sep=""),
	c("estimate","s.e."))
print.default(coef.table, digits=digits,print.gap=2)
if(correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}}
