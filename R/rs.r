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
#     rs2(y, x1, x2, power=c(1,1), weight=rep(1,length(x1)),
#	family=gaussian, iterlim=20)
#     rs3(y, x1, x2, x3, power=c(1,1,1), weight=rep(1,length(x1)),
#	family=gaussian, iterlim=20)
#
#  DESCRIPTION
#
#    Functions to fit two- and three-covariate power-transformed
# response surface models for glms (Box-Tidwell transformation)

### two covariate model
###


#' Two-factor Box-Tidwell Nonlinear Response Surface Models
#' 
#' \code{rs2} fits a two-covariate power-transformed response surface by
#' iterating the function, \code{\link{glm}}.
#' 
#' 
#' @param y Response variable
#' @param x1 First covariate
#' @param x2 Second covariate
#' @param power Initial estimates of the two power transformations
#' @param weight Weight vector
#' @param family glm family
#' @param iterlim Iteration limit
#' @return A list of class, \code{rs}, is returned containing the model and the
#' power estimates.
#' @author J.K. Lindsey
#' @seealso \code{\link{lm}}, \code{\link{glm}}, \code{\link[gnlm]{gnlr}},
#' \code{\link[gnlm]{gnlr3}}, \code{\link[gnlm]{rs3}}
#' @keywords models
#' @examples
#' 
#' x1 <- rep(1:4,5)
#' x2 <- rep(1:5,rep(4,5))
#' y <- rpois(20,1+2*sqrt(x1)+3*log(x2)+4*x1+log(x2)^2+2*sqrt(x1)*log(x2))
#' rs2(y, x1, x2, family=poisson)
#' 
#' @export rs2
rs2 <- function(y, x1, x2, power=c(1,1), weight=rep(1,length(x1)),
	family=gaussian, iterlim=20){
#
# check estimates and data
#
if(any(c(x1,x2)<0))stop("All covariates must be non-negative")
if(length(power)!=2)
	stop("Two estimates of power parameters must be supplied\n")
a <- power[1]
b <- power[2]
#
# iterate, calling glm, to obtain power parameter estimates
#
test <- TRUE
i <- 0
while(test){
	xx1 <- x1^a
	xx2 <- x2^b
	u <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2,family=family,
		weights=weight)
	z1 <- (u$coef[2]*xx1+2*u$coef[4]*xx1^2+u$coef[6]*xx1*xx2)*
		log(ifelse(x1==0,1,x1))
	z2 <- (u$coef[3]*xx2+2*u$coef[5]*xx2^2+u$coef[6]*xx1*xx2)*
		log(ifelse(x2==0,1,x2))
	if(any(is.na(c(z1,z2))))
		stop(paste("NAs in calculating estimates:",a,b))
	u <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2+z1+z2,
		family=family,weights=weight)
	a <- a+u$coef[6]
	b <- b+u$coef[7]
	if(any(is.na(c(a,b))))stop(paste("NAs in calculating estimates:",a,b))
	i <- i+1
	test <- ((u$coef[6]^2>0.00001)||(u$coef[7]^2>0.00001))&&(i<iterlim)}
#
# set up final results
#
z <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2,family=family,weights=weight)
z$df.residual <- z$df.residual-2
z$aic <- z$aic+4
z$powers <- c(a,b)
z$iterations <- i
class(z) <- c("rs",class(z))
return(z)}

### three covariate model
###


#' Three-factor Box-Tidwell Nonlinear Response Surface Models
#' 
#' \code{rs3} fits a three-covariate power-transformed response surface by
#' iterating the function, \code{\link{glm}}.
#' 
#' 
#' @param y Response variable
#' @param x1 First covariate
#' @param x2 Second covariate
#' @param x3 Third covariate
#' @param power Initial estimates of the three power transformations
#' @param weight Weight vector
#' @param family glm family
#' @param iterlim Iteration limit
#' @return A list of class, \code{rs}, is returned containing the model and the
#' power estimates.
#' @author J.K. Lindsey
#' @seealso \code{\link{lm}}, \code{\link{glm}}, \code{\link[gnlm]{gnlr}},
#' \code{\link[gnlm]{gnlr3}}, \code{\link[gnlm]{rs2}}
#' @keywords models
#' @examples
#' 
#' x1 <- rep(1:4,5)
#' x2 <- rep(1:5,rep(4,5))
#' x3 <- c(rep(1:3,6),1,2)
#' #y <- rpois(20,1+2*sqrt(x1)+3*log(x2)+1/x3+4*x1+log(x2)^2+1/x3^2+
#' #	2*sqrt(x1)*log(x2)+sqrt(x1)/x3+log(x2)/x3)
#' y <- c(9,11,14,33,11,19,20,27,22,32,24,24,20,28,25,41,26,31,37,34)
#' rs3(y, x1, x2, x3, family=poisson)
#' 
#' @export rs3
rs3 <- function(y, x1, x2, x3, power=c(1,1,1), weight=rep(1,length(x1)),
	family=gaussian, iterlim=20){
#
# check estimates and data
#
if(any(c(x1,x2,x3)<0))stop("All covariates must be non-negative")
if(length(power)!=3)
	stop("Three estimates of power parameters must be supplied\n")
a <- power[1]
b <- power[2]
d <- power[3]
#
# iterate, calling glm, to obtain power parameter estimates
#
test <- TRUE
i <- 0
while(test){
	xx1 <- x1^a
	xx2 <- x2^b
	xx3 <- x3^d
	xx12 <- xx1*xx2
	xx13 <- xx1*xx3
	xx23 <- xx2*xx3
	u <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+
		xx12+xx13+xx23,family=family,weights=weight)
	z1 <- (u$coef[2]*xx1+2*u$coef[5]*xx1^2+u$coef[8]*xx12+
		+u$coef[9]*xx13)*log(ifelse(x1==0,1,x1))
	z2 <- (u$coef[3]*xx2+2*u$coef[6]*xx2^2+u$coef[8]*xx12+
		u$coef[10]*xx23)*log(ifelse(x2==0,1,x2))
	z3 <- (u$coef[4]*xx3+2*u$coef[7]*xx3^2+u$coef[9]*xx13+
		u$coef[10]*xx23)*log(ifelse(x3==0,1,x3))
	if(any(is.na(c(z1,z2,z3))))
		stop(paste("NAs in calculating estimates:",a,b,d))
	u <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+
		xx12+xx13+xx23+z1+z2+z3,family=family,weights=weight)
	a <- a+u$coef[11]
	b <- b+u$coef[12]
	d <- d+u$coef[13]
	if(any(is.na(c(a,b,d))))stop(paste("NAs in calculating estimates:",a,b,d))
	i <- i+1
	test <- ((u$coef[11]^2>0.00001)||(u$coef[12]^2>0.00001)||
		(u$coef[13]^2>0.00001))&&(i<iterlim)}
#
# set up final results
#
z <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+xx12+xx13+xx23,
	family=family,weights=weight)
z$df.residual <- z$df.residual-3
z$aic <- z$aic+6
z$powers <- c(a,b,d)
z$iterations <- i
class(z) <- c("rs",class(z))
return(z)}

### print and summary methods
###
print.rs <- function(x,...){
cat("\nPowered transformed response surface\n\n")
cat("Powers:",x$powers,"\n")
cat("Iterations:",x$iterations,"\n")
#print.glm(x,...)
class(x) <- "glm"
print(x)}

print.summary.rs <- function(x,...){
cat("\nPowered transformed response surface\n\n")
cat("Powers:",x$powers,"\n")
cat("Iterations:",x$iterations,"\n")
#print.summary.glm(x,...)
class(x) <- "summary.glm"
print(x)}

summary.rs <- function(z,...){
zz <- summary.glm(z,...)
class(zz) <- c("summary.rs",class(zz))
if(!is.null(z$powers))zz$powers <- z$powers
if(!is.null(z$iterations))zz$iterations <- z$iterations
zz}


