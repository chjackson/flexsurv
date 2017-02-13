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
