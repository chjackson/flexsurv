##' The Gompertz distribution
##' 
##' Density, distribution function, hazards, quantile function and random
##' generation for the Gompertz distribution with unrestricted shape.
##' 
##' The Gompertz distribution with \code{shape} parameter \eqn{a} and
##' \code{rate} parameter \eqn{b}{b} has probability density function
##' 
##' \deqn{f(x | a, b) = be^{ax}\exp(-b/a (e^{ax} - 1))}{f(x | a, b) = b exp(ax)
##' exp(-b/a (exp(ax) - 1))}
##' 
##' and hazard
##' 
##' \deqn{h(x | a, b) = b e^{ax}}{h(x | a, b) = b exp(ax)}
##' 
##' The hazard is increasing for shape \eqn{a>0} and decreasing for \eqn{a<0}.
##' For \eqn{a=0} the Gompertz is equivalent to the exponential distribution
##' with constant hazard and rate \eqn{b}.
##' 
##' The probability distribution function is \deqn{F(x | a, b) = 1 - \exp(-b/a
##' (e^{ax} - 1))}{F(x | a, b) = 1 - exp(-b/a (exp(ax) - 1))}
##' 
##' Thus if \eqn{a} is negative, letting \eqn{x} tend to infinity shows that
##' there is a non-zero probability \eqn{\exp(b/a)}{exp(b/a)} of living
##' forever.  On these occasions \code{qgompertz} and \code{rgompertz} will
##' return \code{Inf}.
##' 
##' @aliases Gompertz dgompertz pgompertz qgompertz hgompertz Hgompertz
##' rgompertz
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param shape,rate vector of shape and rate parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dgompertz} gives the density, \code{pgompertz} gives the
##' distribution function, \code{qgompertz} gives the quantile function,
##' \code{hgompertz} gives the hazard function, \code{Hgompertz} gives the
##' cumulative hazard function, and \code{rgompertz} generates random deviates.
##' @note Some implementations of the Gompertz restrict \eqn{a} to be strictly
##' positive, which ensures that the probability of survival decreases to zero
##' as \eqn{x} increases to infinity.  The more flexible implementation given
##' here is consistent with \code{streg} in Stata.
##' 
##' The functions \code{\link[eha]{dgompertz}} and similar available in the
##' package \pkg{eha} label the parameters the other way round, so that what is
##' called the \code{shape} there is called the \code{rate} here, and what is
##' called \code{1 / scale} there is called the \code{shape} here. The
##' terminology here is consistent with the exponential \code{\link{dexp}} and
##' Weibull \code{\link{dweibull}} distributions in R.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{dexp}}
##' @references Stata Press (2007) Stata release 10 manual: Survival analysis
##' and epidemiological tables.
##' @keywords distribution
##' @name Gompertz
NULL


## consistent with paper and Stata 
## shape pos: inc haz, shape neg: dec haz just like weibull, shape zero=exponential   
## 
## in eha: shape=lam, gamma=1/scale
## log(shape) + x/scale - shape * scale * (exp(x/scale) - 1))
## shape/scale labelled wrong way round.

##' @export
##' @rdname Gompertz
dgompertz <- function(x, shape, rate=1, log=FALSE) {
    dgompertz_work(x, shape, rate, log)
}

##' @export
##' @rdname Gompertz
pgompertz <- function(q, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
    pgompertz_work(q, shape, rate, lower.tail, log.p)
}

##' @export
##' @rdname Gompertz
qgompertz <- function(p, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("gompertz", lower.tail=lower.tail, log=log.p, p=p, shape=shape, rate=rate)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    ret[ind][shape==0] <- qexp(p[shape==0], rate=rate[shape==0])
    sn0 <- shape!=0
    if (any(sn0)) { 
        p <- p[sn0]; shape <- shape[sn0]; rate <- rate[sn0]
        asymp <- 1 - exp(rate/shape)
        immortal <- shape < 0 & p > asymp
        ret[ind][sn0][immortal] <- Inf        
        ret[ind][sn0][!immortal] <- 1 / shape[!immortal] *
            log1p(-log1p(-p[!immortal]) * shape[!immortal] / rate[!immortal])
    }
    ret
}

##' @export
##' @rdname Gompertz
rgompertz <- function(n, shape = 1, rate = 1){
    r <- rbase("gompertz", n=n, shape=shape, rate=rate)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind] <- qgompertz(p=runif(sum(ind)), shape=shape, rate=rate)
    ret
}

##' @export
##' @rdname Gompertz
hgompertz <- function(x, shape, rate = 1, log = FALSE) 
{
    h <- dbase("gompertz", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    if (log) 
        ret[ind] <- log(rate) + (shape * x)
    else
        ret[ind] <- rate * exp(shape * x)
    ret
}

##' @export
##' @rdname Gompertz
Hgompertz <- function(x, shape, rate = 1, log = FALSE) 
{
    h <- dbase("gompertz", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- ifelse(shape==0, rate*x, rate/shape * expm1(shape*x))
    if (log) ret[ind] <- log(ret[ind])
    ret
}

##' @export
##' @rdname means
rmst_gompertz = function(t, shape, rate=1, start=0){
  rmst_generic(pgompertz, t, start=start, shape=shape, rate=rate)
}

##' @export
##' @rdname means
mean_gompertz = function(shape, rate = 1){
  rmst_generic(pgompertz, Inf, start=0, shape=shape, rate=rate)
}
