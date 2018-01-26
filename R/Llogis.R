##' The log-logistic distribution
##' 
##' Density, distribution function, hazards, quantile function and
##' random generation for the log-logistic distribution.
##' 
##' The log-logistic distribution with \code{shape} parameter
##' \eqn{a>0} and \code{scale} parameter \eqn{b>0} has probability
##' density function
##' 
##' \deqn{f(x | a, b) = (a/b) (x/b)^{a-1} / (1 + (x/b)^a)^2}
##' 
##' and hazard
##' 
##' \deqn{h(x | a, b) = (a/b) (x/b)^{a-1} / (1 + (x/b)^a)}
##' 
##' for \eqn{x>0}. The hazard is decreasing for shape \eqn{a\leq 1}{a
##' >= 1}, and unimodal for \eqn{a > 1}.
##' 
##' The probability distribution function is \deqn{F(x | a, b) = 1 - 1
##' / (1 + (x/b)^a)}
##' 
##' If \eqn{a > 1}, the mean is \eqn{b c / sin(c)}, and if \eqn{a > 2}
##' the variance is \eqn{b^2 * (2*c/sin(2*c) - c^2/sin(c)^2)}, where
##' \eqn{c = \pi/a}, otherwise these are undefined.
##' 
##' @aliases Llogis dllogis pllogis qllogis hllogis Hllogis rllogis
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the
##' length is taken to be the number required.
##' @param shape,scale vector of shape and scale parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as
##' log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are
##' \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dllogis} gives the density, \code{pllogis} gives the
##' distribution function, \code{qllogis} gives the quantile function,
##' \code{hllogis} gives the hazard function, \code{Hllogis} gives the
##' cumulative hazard function, and \code{rllogis} generates random
##' deviates.
##' @note Various different parameterisations of this distribution are
##' used.  In the one used here, the interpretation of the parameters
##' is the same as in the standard Weibull distribution
##' (\code{\link{dweibull}}).  Like the Weibull, the survivor function
##' is a transformation of \eqn{(x/b)^a} from the non-negative real line to [0,1],
##' but with a different link function.  Covariates on \eqn{b}
##' represent time acceleration factors, or ratios of expected
##' survival.
##' 
##' The same parameterisation is also uqsed in
##' \code{\link[eha]{dllogis}} in the \pkg{eha} package.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{dweibull}}
##' @references Stata Press (2007) Stata release 10 manual: Survival analysis
##' and epidemiological tables.
##' @keywords distribution
##' @name Llogis
NULL

##' @export
dllogis <- function(x, shape=1, scale=1, log = FALSE)
{
    dllogis_work(x, shape, scale, log)
}

##' @export
pllogis <- function(q, shape=1, scale=1, lower.tail = TRUE, log.p = FALSE) {
    pllogis_work(q, shape, scale, lower.tail, log.p)
}

##' @export
qllogis <- function(p, shape=1, scale=1, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("llogis", lower.tail=lower.tail, log=log.p, p=p, shape=shape, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    # lower.tail is handled in dbase() by transforming p to 1-p
    ret[ind] <- exp(qlogis(p, log(scale), 1/shape, lower.tail=TRUE, log.p))
    ret
}

##' @export
rllogis <- function(n, shape=1, scale=1){
    r <- rbase("llogis", n=n, shape=shape, scale=scale)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind] <- qllogis(p=runif(sum(ind)), shape=shape, scale=scale)
    ret
} 

##' @export
hllogis <- function(x, shape=1, scale=1, log = FALSE) 
{
    h <- dbase("llogis", log=log, x=x, shape=shape, scale=scale)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    if (log) ret[ind] <- log(shape) - log(scale) + (shape-1)*(log(x) - log(scale)) - log(1 + (x/scale)^shape)
    else ret[ind] <- (shape/scale) * (x/scale)^{shape-1} / (1 + (x/scale)^shape)
    ret
}

##' @export
Hllogis <- function(x, shape=1, scale=1, log = FALSE) 
{
    ret <- - pllogis(x, shape, scale, lower.tail=FALSE, log.p=TRUE)
    if (log) ret <- log(ret)
    ret
}

DLdllogis <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    ilt <- tss / (1 + tss)
    res[,1] <- 1 + (1 - 2*ilt)*shape*log(t/scale)
    res[,2] <- - shape + 2*ilt*shape
    res
}

DLSllogis <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    ilt <- tss / (1 + tss)    
    res[,1] <- ifelse(t==0, 0, -ilt*log(t/scale)*shape)
    res[,2] <- shape*ilt
    res
}

##' @export
mean_llogis <- function(shape=1, scale=1){
    ifelse(shape > 1,
       {b <- pi/shape
        scale*b / sin(b)},
           NaN)          
}

##' @export
rmst_llogis = function(t, shape=1, scale=1, start=0){
  rmst_generic(pllogis, t, start=start, shape=shape, scale=scale)
}


var.llogis <- function(shape=1, scale=1){
    ifelse(shape > 2,
       {b <- pi/shape
        scale^2 * (2*b/sin(2*b) - b^2/sin(b)^2)},
           NaN)    
}
