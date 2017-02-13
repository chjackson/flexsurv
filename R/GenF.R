## Generalized F distribution (Prentice 1975 parameterisation)
## For P=0 this is equivalent to the generalized (log-)gamma (Prentice 1974)
## P=Q=0, lognormal
## P=0, Q=1, Weibull
## Equation 2 in C.Cox (2008) is wrong, delta*beta*m1 not beta*m1 in first exponent in numerator

dgenf <- function(x, mu=0, sigma=1, Q, P, log=FALSE) {
    dgenf_work(x, mu, sigma, Q, P, log)
}

pgenf <- function(q, mu=0, sigma=1, Q, P, lower.tail = TRUE, log.p = FALSE)
{
    pgenf_work(q, mu, sigma, Q, P, lower.tail, log.p)
}

Hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    -log(pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE))
}

hgenf <- function(x, mu=0, sigma=1, Q, P)
{
    dgenf(x=x, mu=mu, sigma=sigma, Q=Q, P=P) /
        pgenf(q=x, mu=mu, sigma=sigma, Q=Q, P=P, lower.tail=FALSE)
}

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


dgenf.orig <- function(x, mu=0, sigma=1, s1, s2, log=FALSE) {
    d <- dbase("genf.orig", log=log, x=x, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- (log(x) - mu)/sigma
    expw <- x^(1/sigma)*exp(-mu/sigma)
    logdens <- -log(sigma*x) + s1*(log(s1) + w - log(s2)) - (s1+s2)*log(1 + s1*expw/s2) - lbeta(s1, s2)
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

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

Hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    -log(pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE))
}

hgenf.orig <- function(x, mu=0, sigma=1, s1, s2)
{
    dgenf.orig(x=x, mu=mu, sigma=sigma, s1=s1, s2=s2) /
        pgenf.orig(q=x, mu=mu, sigma=sigma, s1=s1, s2=s2, lower.tail=FALSE)
}

qgenf.orig <- function(p, mu=0, sigma=1, s1, s2, lower.tail = TRUE, log.p = FALSE)
{
    d <- dbase("genf.orig", lower.tail=lower.tail, log=log.p, p=p, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    w <- log(qf(p, 2*s1, 2*s2))
    ret[ind] <- exp(w*sigma + mu)
    ret
}

rgenf.orig <- function(n, mu=0, sigma=1, s1, s2)
{
    r <- rbase("genf.orig", n=n, mu=mu, sigma=sigma, s1=s1, s2=s2)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    w <- log(rf(n, 2*s1, 2*s2))
    ret[ind] <- exp(w*sigma + mu)
    ret
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

mean.genf.orig <- function(mu, sigma, s1, s2){
    if (s2 <= sigma) NaN else exp(mu) * (s2/s1)^sigma * gamma(s1 + sigma)*gamma(s2 - sigma) / (gamma(s1)*gamma(s2))
}

