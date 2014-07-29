
## consistent with paper and Stata 
## shape pos: inc haz, shape neg: dec haz just like weibull, shape zero=exponential   
## 
## in eha: shape=lam, gamma=1/scale
## log(shape) + x/scale - shape * scale * (exp(x/scale) - 1))
## shape/scale labelled wrong way round.

dgompertz <- function(x, shape, rate=1, log=FALSE) {
    d <- dbase("gompertz", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    logdens <- numeric(length(x))
    logdens[shape==0] <- dexp(x[shape==0], rate=rate[shape==0], log=TRUE)
    sn0 <- shape!=0
    if (any(sn0)) { 
        x <- x[sn0]; shape <- shape[sn0]; rate <- rate[sn0]
        logdens[sn0] <- log(rate) + shape*x - rate/shape*(exp(shape*x) - 1)
    }    
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

pgompertz <- function(q, shape, rate=1, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("gompertz", lower.tail=lower.tail, log=log.p, q=q, shape=shape, rate=rate)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    prob <- numeric(length(q))
    prob[shape==0] <- pexp(q[shape==0], rate=rate[shape==0])
    sn0 <- shape!=0
    if (any(sn0)) { 
        q <- q[sn0]; shape <- shape[sn0]; rate <- rate[sn0]
        prob[sn0] <- 1 - exp(-rate/shape*(exp(shape*q) - 1))
    }
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
    ret
}

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

rgompertz <- function(n, shape = 1, rate = 1){
    r <- rbase("gompertz", n=n, shape=shape, rate=rate)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind] <- qgompertz(p=runif(sum(ind)), shape=shape, rate=rate)
    ret
}

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

Hgompertz <- function(x, shape, rate = 1, log = FALSE) 
{
    h <- dbase("gompertz", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- ifelse(shape==0, rate*x, rate/shape * expm1(shape*x))
    if (log) ret[ind] <- log(ret[ind])
    ret
}

check.gompertz <- function(shape, rate=1){
    ret <- rep(TRUE, length(shape))
    if (missing(shape)) stop("shape parameter not given")
    if (any(rate<0)) {warning("Negative rate parameter"); ret[rate<0] <- FALSE}
    ret
}
