### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

##' @export
##' @rdname hazard
hlnorm <- function(x, meanlog=0, sdlog=1, log=FALSE){
    h <- dbase("lnorm", log=log, x=x, meanlog=meanlog, sdlog=sdlog)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- dlnorm(x, meanlog, sdlog) / plnorm(x, meanlog, sdlog, lower.tail=FALSE)
    ret
}

##' @export
##' @rdname hazard
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

##' @export
##' @rdname means
mean_lnorm <- function(meanlog=0, sdlog=1){
    exp(meanlog + 0.5*sdlog^2)
}

##' @export
##' @rdname means
rmst_lnorm = function(t, meanlog=0, sdlog=1, start=0){
  rmst_generic(plnorm, t, start=start, meanlog=meanlog, sdlog=sdlog)
}

var.lnorm <- function(meanlog=0, sdlog=1){
    exp(2*meanlog + sdlog^2)*(exp(sdlog^2) - 1)
}
