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




##' Generic function to find restricted mean survival of a distribution
##' 
##' Generic function to find the restricted mean of a distribution, given the
##' equivalent probability distribution function using numeric intergration.
##' 
##' This function is used by default for custom distributions for which an
##' rmst function is not provided.
##' 
##' This assumes a suitably smooth, continuous distribution.
##' 
##' @param pdist Probability distribution function, for example,
##' \code{\link{pnorm}} for the normal distribution, which must be defined in
##' the current workspace.  This should accept and return vectorised parameters
##' and values.  It should also return the correct values for the entire real
##' line, for example a positive distribution should have \code{pdist(x)==0}
##' for \eqn{x<0}.
##' @param t Vector of times to which rmst is evaluated
##' @param start Optional left-truncation time or times.  The returned
##' restricted mean survival will be conditioned on survival up to
##' this time.
##' @param matargs Character vector giving the elements of \code{...} which
##' represent vector parameters of the distribution.  Empty by default.  When
##' vectorised, these will become matrices.  This is used for the arguments
##' \code{gamma} and \code{knots} in \code{\link{qsurvspline}}.
##' @param ...  The remaining arguments define parameters of the distribution
##' \code{pdist}.  These MUST be named explicitly.
##' @return Vector of restricted means survival times of the distribution at
##' \code{p}.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @keywords distribution
##' @examples
##' 
##' rmst_lnorm(500, start=250, meanlog=7.4225, sdlog = 1.1138)
##' rmst_generic(plnorm, 500, start=250, c(0.025, 0.975), meanlog=7.4225, sdlog = 1.1138)
##' # must name the arguments
##' 
##' @export
rmst_generic <- function(pdist, t, start=0, matargs=NULL, ...)
{
  args <- list(...)
  t_len <- length(t)
  if(length(start) == 1) start <- rep(start, length(t))
  ret <- numeric(t_len)
  start_p = 1 - pdist(start,...)
  for(i in seq_len(t_len)){
    ret[i] <- integrate(
      function(end) (1 - pdist(end,...))/ start_p[i], start[i], t[i]
    )$value
  }
  ret[t<start] <- 0
  if (any(is.nan(ret))) warning("NaNs produced")
  ret
}

##' Generic function to find quantiles of a distribution
##' 
##' Generic function to find the quantiles of a distribution, given the
##' equivalent probability distribution function.
##' 
##' This function is used by default for custom distributions for which a
##' quantile function is not provided.
##' 
##' It works by finding the root of the equation \eqn{h(q) = pdist(q) - p = 0}.
##' Starting from the interval \eqn{(-1, 1)}, the interval width is expanded by
##' 50\% until \eqn{h()} is of opposite sign at either end.  The root is then
##' found using \code{\link{uniroot}}.
##' 
##' This assumes a suitably smooth, continuous distribution.
##' 
##' An identical function is provided in the \pkg{msm} package.
##' 
##' @param pdist Probability distribution function, for example,
##' \code{\link{pnorm}} for the normal distribution, which must be defined in
##' the current workspace.  This should accept and return vectorised parameters
##' and values.  It should also return the correct values for the entire real
##' line, for example a positive distribution should have \code{pdist(x)==0}
##' for \eqn{x<0}.
##' @param p Vector of probabilities to find the quantiles for.
##' @param matargs Character vector giving the elements of \code{...} which
##' represent vector parameters of the distribution.  Empty by default.  When
##' vectorised, these will become matrices.  This is used for the arguments
##' \code{gamma} and \code{knots} in \code{\link{qsurvspline}}.
##' @param ...  The remaining arguments define parameters of the distribution
##' \code{pdist}.  These MUST be named explicitly.
##' 
##' This may also contain the standard arguments \code{log.p} (logical; default
##' \code{FALSE}, if \code{TRUE}, probabilities p are given as log(p)), and
##' \code{lower.tail} (logical; if \code{TRUE} (default), probabilities are P[X
##' <= x] otherwise, P[X > x].).
##' 
##' If the distribution is bounded above or below, then this should contain
##' arguments \code{lbound} and \code{ubound} respectively, and these will be
##' returned if \code{p} is 0 or 1 respectively.  Defaults to \code{-Inf} and
##' \code{Inf} respectively.
##' @return Vector of quantiles of the distribution at \code{p}.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @keywords distribution
##' @examples
##' 
##' qnorm(c(0.025, 0.975), 0, 1)
##' qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
##' @export
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
        if (is.matrix(args.mat[[i]])){
            args.mat[[i]] <- matrix(
              apply(args.mat[[i]], 2, function(x)rep(x, length=maxlen)),
              ncol=ncol(args.mat[[i]]),
              byrow=F
            )
        }
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
