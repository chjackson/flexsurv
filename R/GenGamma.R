##' Generalized gamma distribution (original parameterisation)
##' 
##' Density, distribution function, hazards, quantile function and
##' random generation for the generalized gamma distribution, using
##' the original parameterisation from Stacy (1962).
##' 
##' If \eqn{w \sim Gamma(k,1)}{w ~ Gamma(k, 1)}, then \eqn{x =
##' \exp(w/shape + \log(scale))}{x = exp(w/shape + log(scale))}
##' follows the original generalised gamma distribution with the
##' parameterisation given here (Stacy 1962).  Defining
##' \code{shape}\eqn{=b>0}, \code{scale}\eqn{=a>0}, \eqn{x} has
##' probability density
##' 
##' \deqn{f(x | a, b, k) = \frac{b}{\Gamma(k)} \frac{x^{bk - 1}}{a^{bk}} }{ f(x
##' | a, b, k) = (b / \Gamma(k)) (x^{bk -1} / a^{bk}) exp(-(x/a)^b)}\deqn{
##' \exp(-(x/a)^b)}{ f(x | a, b, k) = (b / \Gamma(k)) (x^{bk -1} / a^{bk})
##' exp(-(x/a)^b)}
##' 
##' The original generalized gamma distribution simplifies to the
##' gamma, exponential and Weibull distributions with the following
##' parameterisations:
##' 
##' \tabular{lcl}{ \code{dgengamma.orig(x, shape, scale, k=1)} \tab \code{=}
##' \tab \code{\link{dweibull}(x, shape, scale)} \cr \code{dgengamma.orig(x,
##' shape=1, scale, k)} \tab \code{=} \tab \code{\link{dgamma}(x, shape=k,
##' scale)} \cr \code{dgengamma.orig(x, shape=1, scale, k=1)} \tab \code{=}
##' \tab \code{\link{dexp}(x, rate=1/scale)} \cr }
##' 
##' Also as k tends to infinity, it tends to the log normal (as in
##' \code{\link{dlnorm}}) with the following parameters (Lawless,
##' 1980):
##' 
##' \code{dlnorm(x, meanlog=log(scale) + log(k)/shape,
##' sdlog=1/(shape*sqrt(k)))}
##' 
##' For more stable behaviour as the distribution tends to the log-normal, an
##' alternative parameterisation was developed by Prentice (1974).  This is
##' given in \code{\link{dgengamma}}, and is now preferred for statistical
##' modelling.  It is also more flexible, including a further new class of
##' distributions with negative shape \code{k}.
##' 
##' The generalized F distribution \code{\link{GenF.orig}}, and its similar
##' alternative parameterisation \code{\link{GenF}}, extend the generalized
##' gamma to four parameters.
##' 
##' @aliases GenGamma.orig dgengamma.orig pgengamma.orig qgengamma.orig
##' rgengamma.orig Hgengamma.orig hgengamma.orig
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param shape vector of ``Weibull'' shape parameters.
##' @param scale vector of scale parameters.
##' @param k vector of ``Gamma'' shape parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dgengamma.orig} gives the density, \code{pgengamma.orig}
##' gives the distribution function, \code{qgengamma.orig} gives the quantile
##' function, \code{rgengamma.orig} generates random deviates,
##' \code{Hgengamma.orig} retuns the cumulative hazard and
##' \code{hgengamma.orig} the hazard.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{GenGamma}}, \code{\link{GenF.orig}},
##' \code{\link{GenF}}, \code{\link{Lognormal}}, \code{\link{GammaDist}},
##' \code{\link{Weibull}}.
##' @references Stacy, E. W. (1962). A generalization of the gamma
##' distribution.  Annals of Mathematical Statistics 33:1187-92.
##' 
##' Prentice, R. L. (1974). A log gamma model and its maximum likelihood
##' estimation. Biometrika 61(3):539-544.
##' 
##' Lawless, J. F. (1980). Inference in the generalized gamma and log gamma
##' distributions.  Technometrics 22(3):409-419.
##' @keywords distribution
##' @name GenGamma.orig
NULL



##' Generalized gamma distribution
##' 
##' Density, distribution function, hazards, quantile function and
##' random generation for the generalized gamma distribution, using
##' the parameterisation originating from Prentice (1974). Also known
##' as the (generalized) log-gamma distribution.
##' 
##' If \eqn{\gamma \sim Gamma(Q^{-2}, 1)}{g ~ Gamma(Q^{-2}, 1)} , and \eqn{w =
##' log(Q^2 \gamma) / Q}{w = log(Q^2*g) / Q}, then \eqn{x = \exp(\mu + \sigma
##' w)}{x = exp(mu + sigma w)} follows the generalized gamma distribution with
##' probability density function
##' 
##' \deqn{f(x | \mu, \sigma, Q) = \frac{|Q|(Q^{-2})^{Q^{-2}}}{\sigma x
##' \Gamma(Q^{-2})} \exp(Q^{-2}(Qw - \exp(Qw)))}{ f(x | mu, sigma, Q) = |Q|
##' (Q^{-2})^{Q^{-2}} / / (sigma * x * Gamma(Q^{-2})) exp(Q^{-2}*(Q*w -
##' exp(Q*w)))}
##' 
##' This parameterisation is preferred to the original
##' parameterisation of the generalized gamma by Stacy (1962) since it
##' is more numerically stable near to \eqn{Q=0} (the log-normal
##' distribution), and allows \eqn{Q<=0}.  The original is available
##' in this package as \code{\link{dgengamma.orig}}, for the sake of
##' completion and compatibility with other software - this is
##' implicitly restricted to \code{Q}>0 (or \code{k}>0 in the original
##' notation). The parameters of \code{\link{dgengamma}} and
##' \code{\link{dgengamma.orig}} are related as follows.
##' 
##' \code{dgengamma.orig(x, shape=shape, scale=scale, k=k) = }
##' 
##' \code{dgengamma(x, mu=log(scale) + log(k)/shape, sigma=1/(shape*sqrt(k)),
##' Q=1/sqrt(k))}
##' 
##' The generalized gamma distribution simplifies to the gamma,
##' log-normal and Weibull distributions with the following
##' parameterisations:
##' 
##' \tabular{lcl}{ \code{dgengamma(x, mu, sigma, Q=0)} \tab \code{=} \tab
##' \code{dlnorm(x, mu, sigma)} \cr \code{dgengamma(x, mu, sigma, Q=1)} \tab
##' \code{=} \tab \code{dweibull(x, shape=1/sigma, scale=exp(mu))} \cr
##' \code{dgengamma(x, mu, sigma, Q=sigma)} \tab \code{=} \tab \code{dgamma(x,
##' shape=1/sigma^2, rate=exp(-mu) / sigma^2)} \cr } The properties of the
##' generalized gamma and its applications to survival analysis are discussed
##' in detail by Cox (2007).
##' 
##' The generalized F distribution \code{\link{GenF}} extends the generalized
##' gamma to four parameters.
##' 
##' @aliases GenGamma dgengamma pgengamma qgengamma rgengamma Hgengamma
##' hgengamma
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param mu Vector of ``location'' parameters.
##' @param sigma Vector of ``scale'' parameters.  Note the inconsistent
##' meanings of the term ``scale'' - this parameter is analogous to the
##' (log-scale) standard deviation of the log-normal distribution, ``sdlog'' in
##' \code{\link{dlnorm}}, rather than the ``scale'' parameter of the gamma
##' distribution \code{\link{dgamma}}. Constrained to be positive.
##' @param Q Vector of shape parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dgengamma} gives the density, \code{pgengamma} gives the
##' distribution function, \code{qgengamma} gives the quantile function,
##' \code{rgengamma} generates random deviates, \code{Hgengamma} retuns the
##' cumulative hazard and \code{hgengamma} the hazard.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{GenGamma.orig}}, \code{\link{GenF}},
##' \code{\link{Lognormal}}, \code{\link{GammaDist}}, \code{\link{Weibull}}.
##' @references Prentice, R. L. (1974). A log gamma model and its maximum
##' likelihood estimation. Biometrika 61(3):539-544.
##' 
##' Farewell, V. T. and Prentice, R. L. (1977). A study of
##' distributional shape in life testing. Technometrics 19(1):69-75.
##' 
##' Lawless, J. F. (1980). Inference in the generalized gamma and log
##' gamma distributions.  Technometrics 22(3):409-419.
##' 
##' Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007).
##' Parametric survival analysis and taxonomy of hazard functions for
##' the generalized gamma distribution.  Statistics in Medicine
##' 26:4252-4374
##' 
##' Stacy, E. W. (1962). A generalization of the gamma distribution.
##' Annals of Mathematical Statistics 33:1187-92
##' @keywords distribution
##' @name GenGamma
NULL

## Log-gamma or generalized gamma distribution (parameterisation as in Farewell and Prentice, Technometrics 1977)

### FIXME value for x = 0 

##' @export
##' @rdname GenGamma
dgengamma <- function(x, mu=0, sigma=1, Q, log=FALSE) {
    dgengamma_work(x, mu, sigma, Q, log)
}

##' @export
##' @rdname GenGamma
pgengamma <- function(q, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE) {
    pgengamma_work(q, mu, sigma, Q, lower.tail, log.p)
}

##' @export
##' @rdname GenGamma
Hgengamma <- function(x, mu=0, sigma=1, Q)
{
    -pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail=FALSE, log.p=TRUE)
}

##' @export
##' @rdname GenGamma
hgengamma <- function(x, mu=0, sigma=1, Q)
{
    dgengamma(x=x, mu=mu, sigma=sigma, Q=Q) /
        pgengamma(q=x, mu=mu, sigma=sigma, Q=Q, lower.tail=FALSE)
}

##' @export
##' @rdname GenGamma
qgengamma <- function(p, mu=0, sigma=1, Q, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("gengamma", lower.tail=lower.tail, log=log.p, p=p, mu=mu, sigma=sigma, Q=Q)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    p[Q<0] <- 1 - p[Q<0]
    ret[ind] <- numeric(sum(ind))
    ret[ind][Q==0] <- qlnorm(p[Q==0], mu[Q==0], sigma[Q==0])
    qn0 <- Q!=0
    p <- p[qn0]; mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
    ret[ind][qn0] <- exp(mu + sigma*(log(Q^2*qgamma(p, 1/Q^2, 1)) / Q))
    ret
}

##' @export
##' @rdname GenGamma
rgengamma <- function(n, mu=0, sigma=1, Q) {
    r <- rbase("gengamma", n=n, mu=mu, sigma=sigma, Q=Q)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind][Q==0] <- rlnorm(n, mu, sigma)
    qn0 <- Q!=0
    if (any(qn0)) {
        mu <- mu[qn0]; sigma <- sigma[qn0]; Q <- Q[qn0]
        w <- log(Q^2*rgamma(n, 1/Q^2, 1)) / Q
        ret[ind][qn0] <- exp(mu + sigma*w)
    }
    ret
}

##' @export
##' @rdname means
rmst_gengamma = function(t, mu=0, sigma=1, Q, start=0){
  rmst_generic(pgengamma, t, start=start, mu=mu, sigma=sigma, Q=Q)
}

##' @export
##' @rdname means
mean_gengamma = function(mu=0, sigma=1, Q){
  rmst_generic(pgengamma, Inf, start=0, mu=mu, sigma=sigma, Q=Q)
}

### FIXME limiting value for x=0:  0 if bk >1, 1 if b=k=1, ... ? 

##' @export
##' @rdname GenGamma.orig
dgengamma.orig <- function(x, shape, scale=1, k, log=FALSE){
    d <- dbase("gengamma.orig", log=log, x=x, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    logdens <- log(shape) - lgamma(k) + (shape*k - 1)*log(x) - shape*k*log(scale) - (x/scale)^shape
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

##' @export
##' @rdname GenGamma.orig
pgengamma.orig <- function(q, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("gengamma.orig", lower.tail=lower.tail, log=log.p, q=q, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    y <- log(q)
    w <- (y - log(scale))*shape
    prob <- pgamma(exp(w), shape=k)
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
    ret
}

##' @export
##' @rdname GenGamma.orig
Hgengamma.orig <- function(x, shape, scale=1, k)
{
    -log(pgengamma.orig(q=x, shape=shape, scale=scale, k=k, lower.tail=FALSE))
}

##' @export
##' @rdname GenGamma.orig
hgengamma.orig <- function(x, shape, scale=1, k)
{
    dgengamma.orig(x=x, shape=shape, scale=scale, k=k) /
        pgengamma.orig(q=x, shape=shape, scale=scale, k=k, lower.tail=FALSE)
}

##' @export
##' @rdname GenGamma.orig
qgengamma.orig <- function(p, shape, scale=1, k, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("gengamma.orig", lower.tail=lower.tail, log=log.p, p=p, shape=shape, scale=scale, k=k)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- log(qgamma(p, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}

##' @export
##' @rdname GenGamma.orig
rgengamma.orig <- function(n, shape, scale=1, k) {
    r <- rbase("gengamma.orig", n=n, shape=shape, scale=scale, k=k)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    w <- log(rgamma(n, shape=k))
    y <- w / shape  + log(scale)
    ret[ind] <- exp(y)
    ret
}

##' @export
##' @rdname means
rmst_gengamma.orig = function(t, shape, scale=1, k, start=0){
  rmst_generic(pgengamma.orig, t, start=start, shape=shape, scale=scale, k=k)
}

##' @export
##' @rdname means
mean_gengamma.orig = function(shape, scale=1, k){
  rmst_generic(pgengamma.orig, Inf, start=0, shape=shape, scale=scale, k=k)
}

check.gengamma.orig <- function(shape, scale, k){
    ret <- rep(TRUE, length(shape))
    if (missing(shape)) stop("shape parameter \"shape\" not given")
    if (missing(k)) stop("shape parameter \"k\" not given")
    if (any(shape < 0)) {warning("Negative shape parameter \"shape\""); ret[shape<0] <- FALSE}
    if (any(scale < 0)) {warning("Negative scale parameter"); ret[scale<0] <- FALSE}
    if (any(k < 0)) {warning("Negative shape parameter \"k\""); ret[k<0] <- FALSE}
    ret
}
