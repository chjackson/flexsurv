##' Generalized F distribution (original parameterisation)
##' 
##' Density, distribution function, quantile function and random
##' generation for the generalized F distribution, using the less
##' flexible original parameterisation described by Prentice (1975).
##' 
##' If \eqn{y \sim F(2s_1, 2s_2)}{y ~ F(2*s1, 2*s2)}, and \eqn{w = }{w =
##' log(y)}\eqn{ \log(y)}{w = log(y)} then \eqn{x = \exp(w\sigma + \mu)}{x =
##' exp(w*sigma + mu)} has the original generalized F distribution with
##' location parameter \eqn{\mu}{mu}, scale parameter \eqn{\sigma>0}{sigma>0}
##' and shape parameters \eqn{s_1>0,s_2>0}{s1>0,s2>0}.  The probability density
##' function of \eqn{x} is
##' 
##' \deqn{f(x | \mu, \sigma, s_1, s_2) = \frac{(s_1/s_2)^{s_1} e^{s_1 w}}{\sigma x (1 + s_1 e^w/s_2) ^ {(s_1 + s_2)} B(s_1, s_2)}}{f(x | mu, sigma, s_1, s_2) = ((s1/s2)^{s1} e^{s1 w}) / (sigma x (1 + s1 e^w/s2) ^ (s1 + s2) B(s1, s2))}
##'
##' where \eqn{w = (\log(x) - \mu)/\sigma}{w = (log(x) - mu)/sigma}, and
##' \eqn{B(s_1,s_2) = \Gamma(s_1)\Gamma(s_2)/\Gamma(s_1+s_2)}{B(s1,s2) =
##' \Gamma(s1)\Gamma(s2)/\Gamma(s1+s2)} is the beta function.
##' 
##' As \eqn{s_2 \rightarrow \infty}{s2 -> infinity}, the distribution
##' of \eqn{x} tends towards an original generalized gamma
##' distribution with the following parameters:
##' 
##' \code{\link{dgengamma.orig}(x, shape=1/sigma, scale=exp(mu) /
##' s1^sigma, k=s1)}
##' 
##' See \code{\link{GenGamma.orig}} for how this includes several
##' other common distributions as special cases.
##' 
##' The alternative parameterisation of the generalized F
##' distribution, originating from Prentice (1975) and given in this
##' package as \code{\link{GenF}}, is preferred for statistical
##' modelling, since it is more stable as \eqn{s_1}{s1} tends to
##' infinity, and includes a further new class of distributions with
##' negative first shape parameter.  The original is provided here for
##' the sake of completion and compatibility.
##' 
##' @aliases GenF.orig dgenf.orig pgenf.orig qgenf.orig rgenf.orig Hgenf.orig
##' hgenf.orig
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param mu Vector of location parameters.
##' @param sigma Vector of scale parameters.
##' @param s1 Vector of first F shape parameters.
##' @param s2 vector of second F shape parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dgenf.orig} gives the density, \code{pgenf.orig} gives the
##' distribution function, \code{qgenf.orig} gives the quantile function,
##' \code{rgenf.orig} generates random deviates, \code{Hgenf.orig} retuns the
##' cumulative hazard and \code{hgenf.orig} the hazard.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{GenF}}, \code{\link{GenGamma.orig}},
##' \code{\link{GenGamma}}
##' @references R. L. Prentice (1975). Discrimination among some parametric
##' models. Biometrika 62(3):607-614.
##' @keywords distribution
##' @name GenF.orig
NULL


##' Generalized F distribution
##' 
##' Density, distribution function, hazards, quantile function and
##' random generation for the generalized F distribution, using the
##' reparameterisation by Prentice (1975).
##' 
##' If \eqn{y \sim F(2s_1, 2s_2)}{y ~ F(2*s1, 2*s2)}, and \eqn{w = }{w =
##' log(y)}\eqn{ \log(y)}{w = log(y)} then \eqn{x = \exp(w\sigma + \mu)}{x =
##' exp(w*sigma + mu)} has the original generalized F distribution with
##' location parameter \eqn{\mu}{mu}, scale parameter \eqn{\sigma>0}{sigma>0}
##' and shape parameters \eqn{s_1,s_2}{s1,s2}.
##' 
##' In this more stable version described by Prentice (1975),
##' \eqn{s_1,s_2}{s1,s2} are replaced by shape parameters \eqn{Q,P}, with
##' \eqn{P>0}, and
##' 
##' \deqn{s_1 = 2(Q^2 + 2P + Q\delta)^{-1}, \quad s_2 = 2(Q^2 + 2P -
##' Q\delta)^{-1}}{s1 = 2 / (Q^2 + 2P + Q*delta), s2 = 2 / (Q^2 + 2P -
##' Q*delta)} equivalently \deqn{Q = \left(\frac{1}{s_1} -
##' \frac{1}{s_2}\right)\left(\frac{1}{s_1} + \frac{1}{s_2}\right)^{-1/2},
##' \quad P = \frac{2}{s_1 + s_2} }{Q = (1/s1 - 1/s2) / (1/s1 + 1/s2)^{-1/2}, P
##' = 2 / (s1 + s2) }
##' 
##' Define \eqn{\delta = (Q^2 + 2P)^{1/2}}{delta = (Q^2 + 2P)^{1/2}}, and
##' \eqn{w = (\log(x) - \mu)\delta /\sigma}{w = (log(x) - mu)delta / sigma},
##' then the probability density function of \eqn{x} is \deqn{ f(x) =
##' \frac{\delta (s_1/s_2)^{s_1} e^{s_1 w}}{\sigma x (1 + s_1 e^w/s_2) ^ {(s_1
##' + s_2)} B(s_1, s_2)} }{ f(x) = (delta (s1/s2)^{s1} e^{s1 w}) / (sigma t (1
##' + s1 e^w/s2) ^ {(s1 + s2)} B(s1, s2)) }\deqn{ }{ f(x) = (delta (s1/s2)^{s1}
##' e^{s1 w}) / (sigma t (1 + s1 e^w/s2) ^ {(s1 + s2)} B(s1, s2)) } The
##' original parameterisation is available in this package as
##' \code{\link{dgenf.orig}}, for the sake of completion / compatibility.  With
##' the above definitions,
##' 
##' \code{dgenf(x, mu=mu, sigma=sigma, Q=Q, P=P) = dgenf.orig(x, mu=mu,
##' sigma=sigma/delta, s1=s1, s2=s2)}
##' 
##' The generalized F distribution with \code{P=0} is equivalent to the
##' generalized gamma distribution \code{\link{dgengamma}}, so that
##' \code{dgenf(x, mu, sigma, Q, P=0)} equals \code{dgengamma(x, mu, sigma,
##' Q)}.  The generalized gamma reduces further to several common
##' distributions, as described in the \code{\link{GenGamma}} help page.
##' 
##' The generalized F distribution includes the log-logistic distribution (see
##' \code{\link[eha]{Llogis}}) as a further special case:
##' 
##' \code{dgenf(x, mu=mu, sigma=sigma, Q=0, P=1) = \link[eha]{dllogis}(x,
##' shape=sqrt(2)/sigma, scale=exp(mu))}
##' 
##' The range of hazard trajectories available under this distribution are
##' discussed in detail by Cox (2008).  Jackson et al. (2010) give an
##' application to modelling oral cancer survival for use in a health economic
##' evaluation of screening.
##' 
##' @aliases GenF dgenf pgenf qgenf rgenf Hgenf hgenf
##' @param x,q Vector of quantiles.
##' @param p Vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param mu Vector of location parameters.
##' @param sigma Vector of scale parameters.
##' @param Q Vector of first shape parameters.
##' @param P Vector of second shape parameters.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dgenf} gives the density, \code{pgenf} gives the distribution
##' function, \code{qgenf} gives the quantile function, \code{rgenf} generates
##' random deviates, \code{Hgenf} retuns the cumulative hazard and \code{hgenf}
##' the hazard.
##' @note The parameters \code{Q} and \code{P} are usually called \eqn{q} and
##' \eqn{p} in the literature - they were made upper-case in these R functions
##' to avoid clashing with the conventional arguments \code{q} in the
##' probability function and \code{p} in the quantile function.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{GenF.orig}}, \code{\link{GenGamma}}
##' @references R. L. Prentice (1975). Discrimination among some parametric
##' models. Biometrika 62(3):607-614.
##' 
##' Cox, C. (2008). The generalized \eqn{F} distribution: An umbrella for
##' parametric survival analysis.  Statistics in Medicine 27:4301-4312.
##' 
##' Jackson, C. H. and Sharples, L. D. and Thompson, S. G. (2010).  Survival
##' models in health economic evaluations: balancing fit and parsimony to
##' improve prediction.  International Journal of Biostatistics 6(1):Article
##' 34.
##' @keywords distribution
##' @name GenF
NULL


## Generalized F distribution (Prentice 1975 parameterisation)
## For P=0 this is equivalent to the generalized (log-)gamma (Prentice 1974)
## P=Q=0, lognormal
## P=0, Q=1, Weibull
## Equation 2 in C.Cox (2008) is wrong, delta*beta*m1 not beta*m1 in first exponent in numerator

##' @export
##' @rdname GenF
dgenf <- function(x, mu=0, sigma=1, Q, P, log=FALSE) {
    dgenf_work(x, mu, sigma, Q, P, log)
}
##' @export
##' @rdname GenF
pgenf <- function(q, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
{
    pgenf_work(q, mu, sigma, Q, P, lower.tail, log.p)
}

##' @export
##' @rdname GenF
Hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    -log(pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE))
}

##' @export
##' @rdname GenF
hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P) /
        pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE)
}

##' @export
##' @rdname GenF
qgenf <- function(p, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("genf", lower.tail=lower.tail, log=log.p, p=p, mu=mu, sigma=sigma, Q=Q, P=P)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    ret[ind][P==0] <- qgengamma(p[P==0], mu[P==0], sigma[P==0], Q[P==0])
    pn0 <- P!=0
    if (any(pn0)) {
        mu <- mu[pn0]; sigma <- sigma[pn0]; Q <- Q[pn0]; P <- P[pn0]
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        w <- log(qf(p, 2*s1, 2*s2))
        ret[ind][pn0] <- exp(w*sigma/delta + mu)
    }
    ret
}

##' @export
##' @rdname GenF
rgenf <- function(n, mu=0, sigma=1, Q, P)
{
    r <- rbase("genf", n=n, mu=mu, sigma=sigma, Q=Q, P=P)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind][P==0] <- rgengamma(sum(P==0), mu[P==0], sigma[P==0], Q[P==0])
    pn0 <- P!=0
    if (any(pn0)) {
        mu <- mu[pn0]; sigma <- sigma[pn0]; Q <- Q[pn0]; P <- P[pn0]
        tmp <- Q^2 + 2*P
        delta <- sqrt(tmp)
        s1 <- 2 / (tmp + Q*delta)
        s2 <- 2 / (tmp - Q*delta)
        w <- log(rf(sum(pn0), 2*s1, 2*s2))
        ret[ind][pn0] <- exp(w*sigma/delta + mu)
    }
    ret
}

##' @export
##' @rdname means
rmst_genf= function(t, mu, sigma, Q, P, start=0){
  rmst_generic(pgenf, t, start=start, mu=mu, sigma=sigma, Q=Q, P=P)
}

##' @export
##' @rdname means
mean_genf = function(mu, sigma, Q, P){
  rmst_generic(pgenf, Inf, start=0, mu=mu, sigma=sigma, Q=Q, P=P)
}

##' @export
##' @rdname GenF.orig
dgenf.orig <- function(x, mu=0, sigma=1, s1, s2, log=FALSE) {
    d <- dbase("genf.orig", log=log, x=x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- (log(x) - mu)/sigma
    expw <- x^(1/sigma)*exp(-mu/sigma)
    logdens <- -log(sigma*x) + s1*(log(s1) + w - log(s2)) - (s1+s2)*log(1 + s1*expw/s2) - lbeta(s1, s2)
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

##' @export
##' @rdname GenF.orig
pgenf.orig <- function(q, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("genf.orig", lower.tail=lower.tail, log=log.p, q=q, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- (log(q) - mu)/sigma
    prob <- 1 - pbeta(s2/(s2 + s1*exp(w)), s2, s1)
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
    ret
}

##' @export
##' @rdname GenF.orig
Hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    -log(pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE))
}

##' @export
##' @rdname GenF.orig
hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    dgenf.orig(x=x, mu=mu, sigma=sigma, s1=s1, s2=s2) /
        pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE)
}

##' @export
##' @rdname GenF.orig
qgenf.orig <- function(p, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("genf.orig", lower.tail=lower.tail, log=log.p, p=p, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- log(qf(p, 2*s1, 2*s2))
    ret[ind] <- exp(w*sigma + mu)
    ret
}

##' @export
##' @rdname GenF.orig
rgenf.orig <- function(n, mu=0, sigma=1, s1, s2)
{
    r <- rbase("genf.orig", n=n, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    w <- log(rf(n, 2*s1, 2*s2))
    ret[ind] <- exp(w*sigma + mu)
    ret
}

##' @export
##' @rdname means
rmst_genf.orig= function(t, mu, sigma, s1, s2, start=0){
  rmst_generic(pgenf.orig, t, start=start, mu=mu, sigma=sigma, s1=s1, s2=s2)
}

##' @export
##' @rdname means
mean_genf.orig = function(mu, sigma, s1, s2){
  rmst_generic(pgenf.orig, Inf, start=0, mu=mu, sigma=sigma, s1=s1, s2=s2)
}

check.genf.orig <- function(mu, sigma, s1, s2){
    ret <- rep(TRUE, length(mu))
    if (missing(s1)) stop("shape parameter \"s1\" not given")
    if (missing(s2)) stop("shape parameter \"s2\" not given")
    if (any(sigma < 0)) {warning("Negative scale parameter \"sigma\""); ret[sigma<0] <- FALSE}
    if (any(s1 < 0)) {warning("Negative shape parameter \"s1\""); ret[s1<0] <- FALSE}
    if (any(s2 < 0)) {warning("Negative shape parameter \"s2\""); ret[s2<0] <- FALSE}
    ret
}

## Thanks to Skif Pankov
## currently undocumented and unused!
## Only defined for s2 > sigma
## mean for gengamma.orig will follow
## TODO replace integrate version with this, check equal

#mean_genf.orig <- function(mu, sigma, s1, s2){
#    if (s2 <= sigma) NaN else exp(mu) * (s2/s1)^sigma * gamma(s1 + sigma)*gamma(s2 - sigma) / (gamma(s1)*gamma(s2))
#}

