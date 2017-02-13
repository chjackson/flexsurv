### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

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

hweibull.quiet <- function(x, shape, scale = 1, log = FALSE) {
  ret <- suppressWarnings(hweibull(x=x, shape=shape, scale=scale, log=log))
  ret
}

rweibull.quiet <- rweibull

hweibull.quiet <- function(x, shape, scale = 1, log = FALSE) {
  ret <- suppressWarnings(hweibull(x=x, shape=shape, scale=scale, log=log))
  ret
}

Hweibull.quiet <- function(x, shape, scale = 1, log = FALSE) {
  ret <- suppressWarnings(Hweibull(x=x, shape=shape, scale=scale, log=log))
  ret
}
