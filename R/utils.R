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

qgeneric <- function(pdist, p, matargs=NULL, ...)
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
    ## args containing vector params of the distribution (e.g. gamma and knots in dsurvspline)
    args.mat <- args[matargs]
    args[c(matargs,"lower.tail","log.p","lbound","ubound")] <- NULL
    ## Other args assumed to contain scalar params of the distribution.
    ## Replicate all to their maximum length, along with p 
    matlen <- if(is.null(matargs)) NULL else sapply(args.mat, function(x){if(is.matrix(x))nrow(x) else 1})
    veclen <- sapply(args, length)
    maxlen <- max(c(length(p), veclen, matlen))
    for (i in seq(along=args))
        args[[i]] <- rep(args[[i]], length.out=maxlen)
    for (i in seq(along=args.mat)){
        if (is.matrix(args.mat[[i]]))
            args.mat[[i]] <- apply(args.mat[[i]], 2, function(x)rep(x, length=maxlen))
        else args.mat[[i]] <- matrix(args.mat[[i]], nrow=maxlen, ncol=length(args.mat[[i]]), byrow=TRUE)
    }
    p <- rep(p, length.out=maxlen)

    ret[p < 0 | p > 1] <- NaN
    ind <- (p > 0 & p < 1)
    if (any(ind)) {
        hind <- seq(along=p)[ind]
        h <- function(y) {
            args <- lapply(args, function(x)x[hind[i]])
            args.mat <- lapply(args.mat, function(x)x[hind[i],])
            p <- p[hind[i]]
            args$q <- y
            args <- c(args, args.mat)
            (do.call(pdist, args) - p)
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


## suppresses NOTE from checker about variables created with "assign"
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ind"))
