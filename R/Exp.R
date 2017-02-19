### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)

##' @export
##' @rdname hazard
hexp <- function(x, rate=1, log=FALSE){
    h <- dbase("exp", log=log, x=x, rate=rate)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    ret[ind] <- if (log) log(rate) else rate
    ret
}

##' @export
##' @rdname hazard
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

##' @export
##' @rdname means
mean_exp <- function(rate=1) { 1/rate }

##' @export
##' @rdname means
rmst_exp <- function(t, rate=1,start=0){
  p_start = pexp(start, rate=rate)
  auc_full = mean_exp(rate=rate)
  auc_right = auc_full * (1-pexp(t, rate=rate))
  auc_left = auc_full * p_start
  auc_uncond = auc_full - auc_right - auc_left
  auc_uncond/(1-p_start)
}

var.exp <- function(rate=1) { 1/rate^2 }
