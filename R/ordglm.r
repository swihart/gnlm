#
#  gnlm : A Library of Special Functions for Nonlinear Regression
#  Copyright (C) 1999, 2000, 2001 J.K. Lindsey
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
#     ordglm(formula, data=parent.frame(), link="logit", maxiter=10, weights=1)
#
#  DESCRIPTION
#
#    A function to fit linear regression models for ordinal responses.



#' Generalized Linear Ordinal Regression Models
#' 
#' \code{ordglm} fits linear regression functions with logistic or probit link
#' to ordinal response data by proportional odds.
#' 
#' 
#' @param formula A model formula. The response must be integers numbered from
#' zero to one less than the number of ordered categories.
#' @param data An optional data frame containing the variables in the model.
#' @param link Logit or probit link function.
#' @param maxiter Maximum number of iterations allowed.
#' @param weights A vector containing the frequencies for grouped data.
#' @return A list of class ordglm is returned. The printed output includes the
#' -log likelihood, the corresponding AIC, the deviance, the maximum likelihood
#' estimates, standard errors, and correlations.
#' @author J.K. Lindsey, adapted and heavily modified from Matlab code
#' (ordinalMLE) by J.H. Albert.
#' @seealso \code{\link{glm}}, \code{\link[gnlm]{nordr}}
#' @references Jansen, J. (1991) Fitting regression models to ordinal data.
#' Biometrical Journal 33, 807-815.
#' 
#' Johnson, V.E. and Albert, J.H. (1999) Ordinal Data Modeling.
#' Springer-Verlag.
#' @keywords models
#' @examples
#' 
#' # McCullagh (1980) JRSS B42, 109-142
#' # tonsil size: 2x3 contingency table
#' y <- c(0:2,0:2)
#' carrier <- gl(2,3,6)
#' wt <- c(19,29,24,497,560,269)
#' ordglm(y~carrier, weights=wt)
#' 
#' @export ordglm
ordglm <- function(formula, data=parent.frame(),  link="logit", maxiter=10,
	weights=1){
call <- sys.call()
link <- match.arg(link,c("logit","probit"))
#
# set up response and model matrix
#
mt <- terms(formula)
mf <- model.frame(mt, data, na.action = na.fail)
y <- model.response(mf, "numeric")
if(!is.vector(y,mode="numeric"))
	stop("the response must be a numeric vector with integral values starting at 0")
X <- -model.matrix(mt,mf)
p <- dim(X)[2]
if(p<2)stop("a covariate must be used")
cname <- colnames(X)
n <- length(y)
ncat <- length(unique(y))
np <- ncat-2+p
#
# set up weights
#
if(length(weights)!=1&&length(weights)!=length(y))
	stop("weights must have same length as the response vector")
if(any(weights<0))stop("weights must be non-negative")
if(length(weights)==1)weights <- rep(1,n)
N <- matrix(0,nrow=n,ncol=ncat)
N[cbind(1:n,y+1)] <- weights
#
# check coding of ordinal variable
#
if(ncat<3)stop("ordinal variables must have at least three categories")
if(min(y)!=0||max(y)!=ncat-1)
	stop(paste("ordinal values must be numbered from 0 to ",ncat-1))
else if(any(y!=trunc(y)))stop("ordinal values must be integers")
#
# calculate constants
#
C <- diag(ncat-1)
C[cbind(2:(ncat-1),1:(ncat-2))] <- -1
C <- rbind(C,c(rep(0,ncat-2),-1))
Zt <- diag(ncat-1)
Zmult <- Zt
Zmult[1,1] <- 0
Zt <- Zt[,2:(ncat-1)]
#
# initialize parameters
#
mu <- (weights%o%rep(1,ncat))*(N+0.5)/(weights%o%rep(1,ncat)+ncat/2)+0.1
gam <- matrix(1,nrow=n,ncol=ncat-1)
pi0 <- matrix(0,nrow=n,ncol=ncat)
sumMu <- t(apply(mu+0.1,1,cumsum))
gam[,1:(ncat-1)] <- if(link=="logit")
	log(sumMu[,1:(ncat-1)]/(rep(weights,ncat-1)+ncat/2-sumMu[,1:(ncat-1)]))
	else if(link=="probit")qnorm(sumMu[,1:(ncat-1)]/(rep(weights,ncat-1)+ncat/2))
alpha <- rep(0,np)
alpha[1:(ncat-2)] <- t(mean(gam[,2:(ncat-1)]))
#
# perform iterative weighted least squares
#
logLike0 <- 0
change <- 1
iter <- 0
while(iter<4||(change>0.00001&&iter<=maxiter)){
	A <- matrix(0,nrow=np,ncol=np)
	bb <- rep(0,np)
	for(i in 1:n)if(weights[i]>0){
		tmp1 <- t(cbind(Zt,-rep(1,ncat-1)%o%X[i,]))%*%
			diag(if(link=="logit")exp(gam[i,])/(1+exp(gam[i,]))^2
			else if(link=="probit")dnorm(gam[i,],0,1))%*%t(weights[i]*C)
		tmp2 <- tmp1%*%diag(rep(1,ncat)/mu[i,])
		A <- A + tmp2%*%t(tmp1)
		bb <- bb + tmp2%*%(N[i,]-mu[i,])}
	cov <- solve(A)
	alpha <- alpha + cov%*%bb
	for(i in 1:n)
		gam[i,] <- cbind(Zmult,-rep(1,ncat-1)%o%X[i,])%*%c(0,alpha)
	if(link=="logit"){
		pi0[,1] <- exp(gam[,1])/(1+exp(gam[,1]))
		pi0[,2:(ncat-1)] <- exp(gam[,2:(ncat-1)])/(1+exp(gam[,2:(ncat-1)]))-
			exp(gam[,1:(ncat-2)])/(1+exp(gam[,1:(ncat-2)]))
		pi0[,ncat] <- 1/(1+exp(gam[,ncat-1]))}
	else if(link=="probit"){
		pi0[,1] <- pnorm(gam[,1])
		pi0[,2:(ncat-1)] <- pnorm(gam[,2:(ncat-1)])-pnorm(gam[,1:(ncat-2)])
		pi0[,ncat] <- 1-pnorm(gam[,ncat-1])}
	if(any(pi0<=0)){
#		print("negative probabilities")
		pi0[pi0<=0] <- 1e-2
		pi0 <- pi0/apply(pi0,1,sum)}
	mu <- matrix(rep(weights+0.1,ncat),ncol=ncat)*pi0
	logLike1 <- sum(N*log(pi0+1e-100))
	if(logLike0==0)change <- -logLike1
	else change <- logLike1-logLike0
	logLike0 <- logLike1
	iter <- iter+1}
#
# calculate log likelihood, fitted values, and residuals
#
maxlike <- -logLike0
if(any(weights>1)){
	# check if rectangular contingency table
	rect <- n/ncat==as.integer(n/ncat)
	if(rect){
		ind <- if(y[1]!=y[2])rep(1:(n/ncat),rep(ncat,n/ncat))
			else rep(1:(n/ncat),ncat)
		denom <- capply(weights,ind,sum)
		fitted <- (denom[ind]*pi0)[cbind(1:n,y+1)]
		PearsRes <- (weights-fitted)/sqrt(fitted)
		devRes <- 2*(weights*log(ifelse(weights==0,1,weights)/fitted)-(weights-fitted))
		deviance <- sum(devRes)
		devRes <- sign(weights-fitted)*sqrt(abs(devRes))}
	else fitted <- PearsRes <- devRes <- deviance <- NULL}
else {
	fitted <- pi0
	PearsRes <- devRes <- deviance <- NULL}
#
# calculate se's
#
cov <- cov[c((ncat-1):np,1:(ncat-2)),c((ncat-1):np,1:(ncat-2))]
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
z <- list(
	call=call,
	maxlike=maxlike,
	aic=maxlike+length(alpha),
	coef=alpha[(ncat-1):length(alpha)],
	intercept=alpha[1:(ncat-2)],
	se=se,
	cov=cov,
	corr=corr,
	deviance=deviance,
	residuals=PearsRes,
	dev.res=devRes,
	df=n-length(alpha),
	fitted=fitted,
	iterations=iter,
	states=ncat,
	link=link,
	cname=cname)
class(z) <- "ordglm"
z}

### standard methods
###
fitted.ordglm <- function(z,...) z$fitted

residuals.ordglm <- function(z,type=c("deviance","pearson")){
type <- match.arg(type)
if(type=="deviance")z$dev.res
else z$residuals}

print.ordglm <- function(x,digits=max(3,.Options$digits-3),correlation=TRUE,...){
  z <- x
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
cat("-Log likelihood   ",z$maxlike,"\n")
cat("AIC               ",z$aic,"\n")
if(!is.null(z$deviance)){
	cat("Deviance          ",z$deviance,"\n")
	cat("df                ",z$df,"\n")}
cat("Iterations        ",z$iterations,"\n")
cat("\nLocation coefficients\n")
#coef.table <- cbind(z$coef,z$se[(z$states-1):length(z$se)])
coef.table <- cbind(z$coef,z$se[1:length(z$coef)])
dimnames(coef.table) <- list(z$cname,c("estimate","s.e."))
print.default(coef.table, digits=digits, print.gap=2)
cat("\nIntercept contrasts\n")
#coef.table <- cbind(z$intercept,z$se[1:(z$states-2)])
coef.table <- cbind(z$intercept,z$se[(length(z$coef)+1):length(z$se)])
dimnames(coef.table) <- list(paste("b[",2:(z$states-1),"]",sep=""),
		     c("estimate","s.e."))
print.default(coef.table, digits=digits, print.gap=2)
if(correlation){
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)}}
