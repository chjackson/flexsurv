### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

##' @export
##' @rdname hazard
hlnorm <- function(x, meanlog=0, sdlog=1, log=FALSE){
    h <- dbase("lnorm", log=log, x=x, meanlog=meanlog, sdlog=sdlog)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    logdens <- dlnorm(x = x, meanlog = meanlog, sdlog = sdlog, log=TRUE)
    logsurv <- plnorm(q = x, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE, log.p=TRUE)
    loghaz <- logdens - logsurv
    ret[ind] <- if (log) loghaz else exp(loghaz)
    ret
}

##' @export
##' @rdname hazard
Hlnorm <- function(x, meanlog=0, sdlog=1, log=FALSE){
    h <- dbase("lnorm", log=log, x=x, meanlog=meanlog, sdlog=sdlog)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- - plnorm(x, meanlog, sdlog, lower.tail=FALSE, log.p=TRUE)
    if (log) ret[ind] <- log(ret[ind])
    ret
}

check.lnorm <- function(meanlog=0, sdlog=1){
    ret <- rep(TRUE, length(meanlog))
    if (any(!is.na(sdlog) & sdlog<0)) {
        warning("Negative SD parameter")
        ret[!is.na(sdlog) & sdlog<0] <- FALSE
    }
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
