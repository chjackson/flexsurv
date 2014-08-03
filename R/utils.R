### Standardised procedure for defining density, cumulative
### distribution, hazard and cumulative hazard functions for
### time-to-event distributions

dbase <- function(dname, lower.tail=TRUE, log=FALSE, ...){
    args <- list(...)
    ## Vectorise all arguments, replicating to length of longest argument
    n <- max(sapply(args, length))
    for (i in seq_along(args)) {
        args[[i]] <- rep(args[[i]], length=n)
    }
    ret <- numeric(n)
    ## Check for parameters out of range, give warning and return NaN
    ## for those
    check.fn <- paste("check.",dname,sep="")
    check.ret <- do.call(check.fn, args[-1])
    ret[!check.ret] <- NaN
    for (i in seq_along(args))
        ret[is.nan(args[[i]])] <- NaN
    ## name of first arg is x for PDF, haz, or cum haz, q for CDF and p for quantile function
    stopifnot( !(names(args)[1]=="x" && lower.tail==FALSE))
    if (names(args)[1] %in% c("x","q")){
        x <- args[[1]]
        ## PDF, CDF, hazard and cumulative hazard is 0 for any negative time
        ret[!is.nan(ret) & (x<0)] <- if (lower.tail) { if (log) -Inf else 0 } else { if (log) 0 else 1 }
    }
    if (names(args)[1] == "p") {
        p <- args[[1]]
        if (log) p <- exp(p)
        if (!lower.tail) p <- 1 - p
        args[[1]] <- p
        ret[p < 0 | p > 1] <- NaN
        ## should be 0,Inf for p=0,1, but hopefully always handled anyway
        ## Result is NA if x or a parameter is NA
    }
    ## Result is NA if x or a parameter is NA
    nas <- rep(FALSE, n)
    for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
    ret[nas] <- NA
    ind <- !is.nan(ret) & !nas
    if (names(args)[1] %in% c("x", "q")) ind <- ind & (x>=0)
    ## Any remaining elements of vector are filled in by standard
    ## formula for hazard
    li <- list(ret=ret, ind=ind)
    for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
    c(li, args)
}

### Standardised procedure for defining random sampling functions

rbase <- function(dname, n, ...){
    ## Vectorise all arguments, replicating to sample length
    if (length(n) > 1) n <- length(n)
    args <- list(...)
    for (i in seq_along(args)) {
        args[[i]] <- rep(args[[i]], length=n)
    }
    ret <- numeric(n)
    ## Check for parameters out of range, give warning and return NaN
    ## for those
    check.fn <- paste("check.",dname,sep="")
    check.ret <- do.call(check.fn, args)
    ret[!check.ret] <- NaN
    for (i in seq_along(args))
        ret[is.nan(args[[i]])] <- NaN
    nas <- rep(FALSE, n)
    for (i in seq_along(args)) nas <- nas | (is.na(args[[i]]) & !is.nan(args[[i]]))
    ret[nas] <- NA
    ind <- !is.nan(ret) & !nas
    li <- list(ret=ret, ind=ind)
    for(i in seq_along(args)) args[[i]] <- args[[i]][ind]
    c(li, args)
}

### Quantile function "qdist" for a generic distribution "dist".
### Works using an interval search for the solution of pdist(q) = p.
### Requires a probability function "pdist" for the same distribution
### in the working environment.

qgeneric <- function(pdist, p, ...)
{
    args <- list(...)
    if (is.null(args$log.p)) args$log.p <- FALSE
    if (is.null(args$lower.tail)) args$lower.tail <- TRUE
    if (is.null(args$lbound)) args$lbound <- -Inf
    if (is.null(args$ubound)) args$ubound <- Inf
    if (args$log.p) p <- exp(p)
    if (!args$lower.tail) p <- 1 - p
    ret <- numeric(length(p))
    ret[p == 0] <- args$lbound
    ret[p == 1] <- args$ubound
    args[c("lower.tail","log.p","lbound","ubound")] <- NULL
    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq(along=p)[ind]
        h <- function(y) {
            args$q <- y
            (do.call(pdist, args) - p)[hind[i]]
        }
        ptmp <- numeric(length(p[ind]))
        for (i in 1:length(p[ind])) {
            interval <- c(-1, 1)
            while (h(interval[1])*h(interval[2]) >= 0) {
                interval <- interval + c(-1,1)*0.5*(interval[2]-interval[1])
            }
            ptmp[i] <- uniroot(h, interval, tol=.Machine$double.eps)$root
        }
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}

### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

hexp <- function(x, rate=1, log=FALSE){
    h <- dbase("exp", log=log, x=x, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- if (log) log(rate) else rate
    ret
}

Hexp <- function(x, rate=1, log=FALSE){
    h <- dbase("exp", log=log, x=x, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- if (log) {log(rate) + log(x)} else rate*x
    ret
}

check.exp <- function(rate=1){
    ret <- rep(TRUE, length(rate))
    if (any(rate<0)) {warning("Negative rate parameter"); ret[rate<0] <- FALSE}
    ret
}

mean.exp <- function(rate=1) { 1/rate }

var.exp <- function(rate=1) { 1/rate^2 }

### Weibull

hweibull <- function (x, shape, scale = 1, log = FALSE)
{
    h <- dbase("weibull", log=log, x=x, shape=shape, scale=scale)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    if (log)
        ret[ind] <- log(shape) + (shape-1)*log(x/scale) - log(scale)
    else
        ret[ind] <- shape * (x/scale)^(shape - 1)/scale
    ret
}

Hweibull <- function (x, shape, scale = 1, log = FALSE)
{
    h <- dbase("weibull", log=log, x=x, shape=shape, scale=scale)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    if (log)
        ret[ind] <- shape * log(x/scale)
    else
        ret[ind] <- (x/scale)^shape
    ret
}

check.weibull <- function(shape, scale=1){
    ret <- rep(TRUE, length(shape))
    if (any(shape<0)) {warning("Negative shape parameter"); ret[shape<0] <- FALSE}
    if (any(scale<0)) {warning("Negative scale parameter"); ret[scale<0] <- FALSE}
    ret
}

mean.weibull <- function(shape, scale=1) { scale * gamma(1 + 1/shape) }

var.weibull <- function(shape, scale=1) { scale^2 * (gamma(1 + 2/shape) - (gamma(1 + 1/shape))^2) }

## Gamma

hgamma <- function(x, shape, rate=1, log=FALSE){
    h <- dbase("gamma", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- dgamma(x, shape, rate) / pgamma(x, shape, rate, lower.tail=FALSE)
    ret
}

Hgamma <- function(x, shape, rate=1, log=FALSE){
    h <- dbase("gamma", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- - pgamma(x, shape, rate, lower.tail=FALSE, log.p=TRUE)
    ret
}

check.gamma <- function(shape, rate=1){
    ret <- rep(TRUE, length(shape))
    if (any(shape<0)) {warning("Negative shape parameter"); ret[shape<0] <- FALSE}
    if (any(rate<0)) {warning("Negative scale parameter"); ret[rate<0] <- FALSE}
    ret
}

mean.gamma <- function(shape, rate=1) {shape / rate}

var.gamma <- function(shape, rate=1) {shape / rate^2}

## lognormal

hlnorm <- function(x, meanlog=0, sdlog=1, log=FALSE){
    h <- dbase("lnorm", log=log, x=x, meanlog=meanlog, sdlog=sdlog)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- dlnorm(x, meanlog, sdlog) / plnorm(x, meanlog, sdlog, lower.tail=FALSE)
    ret
}

Hlnorm <- function(x, meanlog=0, sdlog=1, log=FALSE){
    h <- dbase("lnorm", log=log, x=x, meanlog=meanlog, sdlog=sdlog)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- - plnorm(x, meanlog, sdlog, lower.tail=FALSE, log.p=TRUE)
    ret
}

check.lnorm <- function(meanlog=0, sdlog=1){
    ret <- rep(TRUE, length(meanlog))
    if (any(sdlog<0)) {warning("Negative SD parameter"); ret[sdlog<0] <- FALSE}
    ret
}

mean.lnorm <- function(meanlog=0, sdlog=1){
    exp(meanlog + 0.5*sdlog^2)
}

var.lnorm <- function(meanlog=0, sdlog=1){
    exp(2*meanlog + sdlog^2)*(exp(sdlog^2) - 1)
}

## Don't warn about NaNs when NaNs are produced.  This happens for
## extreme parameter values during optimisation.

dweibull.quiet <- function(x, shape, scale = 1, log = FALSE) {
    ret <- suppressWarnings(dweibull(x=x, shape=shape, scale=scale, log=log))
    ret
}

pweibull.quiet <- function(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    ret <- suppressWarnings(pweibull(q=q, shape=shape, scale=scale, lower.tail=lower.tail, log.p=log.p))
    ret
}

## suppresses NOTE from checker about variables created with "assign"
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ind"))

