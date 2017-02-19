### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

##' @export
##' @rdname hazard
hgamma <- function(x, shape, rate=1, log=FALSE){
    h <- dbase("gamma", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- dgamma(x, shape, rate) / pgamma(x, shape, rate, lower.tail=FALSE)
    ret
}

##' @export
##' @rdname hazard
Hgamma <- function(x, shape, rate=1, log=FALSE){
    h <- dbase("gamma", log=log, x=x, shape=shape, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- - pgamma(x, shape, rate, lower.tail=FALSE, log.p=TRUE)
    ret
}

check.gamma <- function(shape, rate=1){
    ret <- rep(TRUE, length(shape))
    if (any(shape<0)) {
        warning("Negative shape parameter")
        ret[shape<0] <- FALSE
    }
    if (any(rate<0)) {
        warning("Negative scale parameter")
        ret[rate<0] <- FALSE
    }
    ret
}

##' @export
##' @rdname means
mean_gamma <- function(shape, rate=1) {shape / rate}

var.gamma <- function(shape, rate=1) {shape / rate^2}

##' @export
##' @rdname means
rmst_gamma = function(t, shape, rate=1, start=0){
  rmst_generic(pgamma, t, start=start, shape=shape, rate=rate)
}
