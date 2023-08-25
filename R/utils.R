### Standardised procedure for defining density, cumulative
### distribution, hazard and cumulative hazard functions for
### time-to-event distributions

dbase <- function(dname, lower.tail=TRUE, log=FALSE, ...){
    args <- list(...)
    ## Vectorise all arguments, replicating to length of longest argument
    n <- max(sapply(args, length))
    for (i in seq_along(args)) {
        args[[i]] <- rep(args[[i]], length.out=n)
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
        args[[i]] <- rep(args[[i]], length.out=n)
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


##' Generic function to find restricted mean survival time for some distribution
##'
##' Generic function to find the restricted mean of a distribution, given the
##' equivalent probability distribution function, using numeric integration.
##'
##' This function is used by default for custom distributions for which an
##' \code{rmst} function is not provided.
##'
##' This assumes a suitably smooth, continuous distribution.
##'
##' @param pdist Probability distribution function, for example,
##' \code{\link{pnorm}} for the normal distribution, which must be defined in
##' the current workspace.  This should accept and return vectorised parameters
##' and values.  It should also return the correct values for the entire real
##' line, for example a positive distribution should have \code{pdist(x)==0}
##' for \eqn{x<0}.
##'
##' @param t Vector of times at which rmst is evaluated
##'
##' @param start Optional left-truncation time or times.  The returned
##' restricted mean survival will be conditioned on survival up to
##' this time.
##'
##' @param matargs Character vector giving the elements of \code{...} which
##' represent vector parameters of the distribution.  Empty by default.  When
##' vectorised, these will become matrices.  This is used for the arguments
##' \code{gamma} and \code{knots} in \code{\link{psurvspline}}.
##'
##' @param scalarargs Character vector naming scalar arguments of the distribution function that cannot be vectorised.  This is used for the arguments \code{scale} and \code{timescale} in \code{\link{psurvspline}}.
##'
##' @param ...  The remaining arguments define parameters of the distribution
##' \code{pdist}.  These MUST be named explicitly.
##'
##' @return Vector of restricted mean survival times of the distribution at
##' \code{p}.
##'
##' @author Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
##'
##' @keywords distribution
##'
##' @examples
##'
##' rmst_lnorm(500, start=250, meanlog=7.4225, sdlog = 1.1138)
##' rmst_generic(plnorm, 500, start=250, meanlog=7.4225, sdlog = 1.1138)
##' # must name the arguments
##'
##' @export
rmst_generic <- function(pdist, t, start=0, matargs=NULL, scalarargs=NULL, ...)
{
  args <- list(...)
  args_mat <- args[matargs]
  args_scalar <- args[scalarargs]
  args[c(matargs,scalarargs)] <- NULL
  matlen <- if(is.null(matargs)) NULL else sapply(args_mat, function(x){if(is.matrix(x))nrow(x) else 1})
  veclen <- if (length(args) == 0) NULL else sapply(args, length)
  t_len <- length(t)
  maxlen <- max(c(t_len, veclen, matlen))  
  if(length(start) == 1) start <- rep(start, length.out=maxlen)
  na_inds <- rep(FALSE, maxlen)
  for (i in seq_along(args)){
      args[[i]] <- rep(args[[i]], length.out=maxlen)
      na_inds <- na_inds | is.na(args[[i]])
  }
  t <- rep(t, length.out=maxlen)
  for (i in seq_along(args_mat)){
      if (is.matrix(args_mat[[i]])){
          args_mat[[i]] <- matrix(
              apply(args_mat[[i]], 2, function(x)rep(x, length.out=maxlen)),
              ncol=ncol(args_mat[[i]]),
              byrow=F
          )
      }
      else args_mat[[i]] <- matrix(args_mat[[i]], nrow=maxlen, ncol=length(args_mat[[i]]), byrow=TRUE)
      na_inds <- na_inds | apply(args_mat[[i]], 1, anyNA)
  }
  ret <- numeric(maxlen)
  ret[na_inds] <- NA
  for (i in seq_len(maxlen)[!na_inds]){
      fargs_vec <- lapply(args, function(x)x[i])
      fargs_mat <- lapply(args_mat, function(x)x[i,,drop=FALSE])
      pdargs <- c(list(start[i]), fargs_vec, fargs_mat, args_scalar)
      start_p <- 1 - do.call(pdist, pdargs)
      fn <- function(end){
          pdargs <- c(list(end), fargs_vec, fargs_mat, args_scalar)
          pd <- do.call(pdist, pdargs)
          (1 - pd) / start_p
      }
      ret[i] <- integrate(fn, start[i], t[i])$value
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
##' @param pdist Probability distribution function, for example,
##' \code{\link{pnorm}} for the normal distribution, which must be defined in
##' the current workspace.  This should accept and return vectorised parameters
##' and values.  It should also return the correct values for the entire real
##' line, for example a positive distribution should have \code{pdist(x)==0}
##' for \eqn{x<0}.
##'
##' @param p Vector of probabilities to find the quantiles for.
##'
##' @param matargs Character vector giving the elements of \code{...} which
##' represent vector parameters of the distribution.  Empty by default.  When
##' vectorised, these will become matrices.  This is used for the arguments
##' \code{gamma} and \code{knots} in \code{\link{qsurvspline}}.
##'
##' @param scalarargs Character vector naming scalar arguments of the distribution function that cannot be vectorised.  This is used for the arguments \code{scale} and \code{timescale} in \code{\link{qsurvspline}}.
##'
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
##'
##' @return Vector of quantiles of the distribution at \code{p}.
##'
##' @author Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk>
##'
##' @keywords distribution
##'
##' @examples
##'
##' qnorm(c(0.025, 0.975), 0, 1)
##' qgeneric(pnorm, c(0.025, 0.975), mean=0, sd=1) # must name the arguments
##' @export
qgeneric <- function(pdist, p, matargs=NULL, scalarargs=NULL, ...)
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
    ## Arguments that cannot be vectorised
    args.scalar <- args[scalarargs]
    args[c(matargs,scalarargs,"lower.tail","log.p","lbound","ubound")] <- NULL
    ## Other args assumed to contain vectorisable parameters of the distribution.
    ## Replicate all to their maximum length, along with p
    matlen <- if(is.null(matargs)) NULL else sapply(args.mat, function(x){if(is.matrix(x))nrow(x) else 1})
    veclen <- if (length(args) == 0) NULL else sapply(args, length)
    maxlen <- max(c(length(p), veclen, matlen))
    na_inds <- rep(FALSE, length(ret))
    for (i in seq_along(args)){
        args[[i]] <- rep(args[[i]], length.out=maxlen)
        na_inds <- na_inds | is.na(args[[i]])
    }
    for (i in seq_along(args.mat)){
        if (is.matrix(args.mat[[i]])){
            args.mat[[i]] <- matrix(
              apply(args.mat[[i]], 2, function(x)rep(x, length.out=maxlen)),
              ncol=ncol(args.mat[[i]]),
              byrow=F
            )
        }
        else args.mat[[i]] <- matrix(args.mat[[i]], nrow=maxlen, ncol=length(args.mat[[i]]), byrow=TRUE)
        na_inds <- na_inds | apply(args.mat[[i]], 1, anyNA)
    }
    p <- rep(p, length.out=maxlen)
    ret[p < 0 | p > 1] <- NaN
    ret[na_inds] <- NA
    ind <- (p > 0 & p < 1 & !na_inds)
    if (any(ind)) {
        hind <- seq_along(p)[ind]
        n <- length(p[ind])
        ptmp <- numeric(n)
        interval <- matrix(rep(c(-1, 1), n), ncol=2, byrow=TRUE)
        h <- function(y) {
            args <- lapply(args, function(x)x[hind])
            args.mat <- lapply(args.mat, function(x)x[hind,])
            p <- p[hind]
            args$q <- y
            args <- c(args, args.mat, args.scalar)
            (do.call(pdist, args) - p)
        }
        ptmp <- rstpm2::vuniroot(h, interval, tol=.Machine$double.eps, extendInt="yes", maxiter=10000)$root
        ret[ind] <- ptmp
    }
    if (any(is.nan(ret))) warning("NaNs produced")
    ret
}


## suppresses NOTE from checker about variables created with "assign"
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ind"))

##' helper function to safely convert a Hessian matrix to covariance matrix
##'
##' @param hessian hessian matrix to convert to covariance matrix (must be evaluated at MLE)
##' @param tol.solve tolerance used for solve()
##' @param tol.evalues accepted tolerance for negative eigenvalues of the covariance matrix
##' @param ... arguments passed to Matrix::nearPD
##'
##' @importFrom Matrix nearPD
##' @keywords internal
.hess_to_cov <- function(hessian, tol.solve = 1e-9, tol.evalues = 1e-5, ...) {
  if(is.null(tol.solve)) tol.solve <- .Machine$double.eps
  if(is.null(tol.evalues)) tol.evalues <- 1e-5 
  # use solve(.) over chol2inv(chol(.)) to get an inverse even if not PD
  # less efficient but more stable
  inv_hessian <- solve(hessian, tol = tol.solve)
  if (any(is.infinite(inv_hessian)))
    stop("Inverse Hessian has infinite values.  This might indicate that the model is too complex to be identifiable from the data")
  evalues <- eigen(inv_hessian, symmetric = TRUE, only.values = TRUE)$values
  if (min(evalues) < -tol.evalues)
    warning(sprintf(
      "Hessian not positive definite: smallest eigenvalue is %.1e (threshold: %.1e). This might indicate that the optimization did not converge to the maximum likelihood, so that the results are invalid. Continuing with the nearest positive definite approximation of the covariance matrix.",
      min(evalues), -tol.evalues
    ))
  # make sure we return a plain positive definite symmetric matrix
  as.matrix(Matrix::nearPD(inv_hessian, ensureSymmetry = TRUE, ...)$mat)
}


#' Numerical evaluation of the hessian of a function using numDeriv::hessian
#'
#' We perform a quick check about the expected runtime and adjust the
#' precision accordingly.
#'
#' @param f function to compute Hessian for
#' @param x location to evaluate Hessian at
#' @param seconds.warning time threshold in seconds to trigger message and
#'    reduce the number of iterations for Richardson extrapolation of
#'    numDeriv::hessian
#' @param default.r default number of iterations (high-precision recommendation
#'    of numDeriv)
#' @param min.r minial number of iteration, must be at least 2,
#' @param ... further arguments passed to method.args of numDeriv::hessian
#'
#' @importFrom numDeriv hessian
#' @keywords internal
.hessian <- function(f, x, seconds.warning = 60, default.r = 6, min.r = 2, ...) {
  # dimensionality of the problem
  k <- length(x)
  # estimate evaluation time of f(x)
  start_time <- Sys.time()
  for (i in 1:3) f(x)
    mean_eval_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))/3
  default_runtime <- mean_eval_sec * (1 + default.r*(k^2 + k))
  # reduce iterations for Richardson extrapolation
  if (default_runtime > seconds.warning) {
          r <- default.r
    runtime <- default_runtime
    while ((r > min.r) & (runtime > seconds.warning)) {
            r <- r - 1
      runtime <- mean_eval_sec * (1 + r*(k^2 + k))
    }
    message(sprintf(
      "estimated runtime for evaluating the hessian with r=%i is %i minutes, reducing r to %i, estimated runtime is %i minutes",
      default.r, round(default_runtime/60), r, round(runtime/60)
    ))
  } else {
    r <- default.r
  }
  numDeriv::hessian(f, x, method = "Richardson", method.args = list(r = r, ...))
}

check_numeric <- function(...){
  args <- list(...)
  nm <- names(args)
  for (i in seq_along(args)){
    nm <- names(args)[i]
    nmstr <- if (is.null(nm) || nm=="") "" else sprintf(" for `%s`", nm)
    if (!is.numeric(args[[i]]))
        stop(sprintf("Non-numeric value supplied%s", nmstr))
  }
}
