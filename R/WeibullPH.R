##' Weibull distribution in proportional hazards parameterisation
##' 
##' Density, distribution function, hazards, quantile function and random
##' generation for the Weibull distribution in its proportional hazards
##' parameterisation.
##' 
##' The Weibull distribution in proportional hazards parameterisation with
##' `shape' parameter a and `scale' parameter m has density given by
##' 
##' \deqn{f(x) = a m x^{a-1} exp(- m x^a) }
##' 
##' cumulative distribution function \eqn{F(x) = 1 - exp( -m x^a )}, survivor
##' function \eqn{S(x) = exp( -m x^a )}, cumulative hazard \eqn{m x^a} and
##' hazard \eqn{a m x^{a-1}}.
##' 
##' \code{\link{dweibull}} in base R has the alternative 'accelerated failure
##' time' (AFT) parameterisation with shape a and scale b.  The shape parameter
##' \eqn{a} is the same in both versions.  The scale parameters are related as
##' \eqn{b = m^{-1/a}}, equivalently m = b^-a.
##' 
##' In survival modelling, covariates are typically included through a linear
##' model on the log scale parameter.  Thus, in the proportional hazards model,
##' the coefficients in such a model on \eqn{m} are interpreted as log hazard
##' ratios.
##' 
##' In the AFT model, covariates on \eqn{b} are interpreted as time
##' acceleration factors.  For example, doubling the value of a covariate with
##' coefficient \eqn{beta=log(2)} would give half the expected survival time.
##' These coefficients are related to the log hazard ratios \eqn{\gamma} as
##' \eqn{\beta = -\gamma / a}.
##' 
##' @aliases WeibullPH dweibullPH pweibullPH qweibullPH rweibullPH HweibullPH
##' hweibullPH
##' @param x,q Vector of quantiles.
##' @param p Vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param shape Vector of shape parameters.
##' @param scale Vector of scale parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dweibullPH} gives the density, \code{pweibullPH} gives the
##' distribution function, \code{qweibullPH} gives the quantile function,
##' \code{rweibullPH} generates random deviates, \code{HweibullPH} retuns the
##' cumulative hazard and \code{hweibullPH} the hazard.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{dweibull}}
##' @keywords distribution
##' @name WeibullPH
NULL

##' @export
dweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    dweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

##' @export
pweibullPH <- function(q, shape, scale = 1,
                       lower.tail=TRUE, log.p=FALSE) {
    pweibull(q, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}

##' @export
qweibullPH <- function(p, shape, scale = 1, lower.tail=TRUE, log.p=FALSE) {
    qweibull(p, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}

##' @export
hweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

##' @export
HweibullPH <- function(x, shape, scale=1, log=FALSE) {
    Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

##' @export
rweibullPH <- function(n, shape, scale=1) {
    rweibull(n, shape=shape, scale=scale^{-1/shape})
}

##' @export
rmst_weibullPH = function(t, shape, scale=1, start=0){
  rmst_generic(pweibullPH, t, start=start, shape=shape, scale=scale)
}

##' @export
mean_weibullPH = function(shape, scale=1){
  mean_weibull(shape=shape, scale=scale^{-1/shape})
}
