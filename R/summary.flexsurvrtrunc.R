##' Summarise quantities of interest from fitted flexsurvrtrunc models
##' 
##' This function extracts quantities of interest from the untruncated 
##' version of a model with individual-specific right truncation points 
##' fitted by \code{\link{flexsurvrtrunc}}.  Note that covariates are
##' currently not supported by \code{\link{flexsurvrtrunc}}.
##' 
##' 
##' @param type \code{"survival"} for survival probabilities.
##' 
##' \code{"cumhaz"} for cumulative hazards.
##' 
##' \code{"hazard"} for hazards.
##' 
##' \code{"rmst"} for restricted mean survival.
##' 
##' \code{"mean"} for mean survival.
##' 
##' \code{"median"} for median survival (alternative to \code{type="quantile"} with \code{quantiles=0.5}).
##' 
##' \code{"quantile"} for quantiles of the survival time distribution.
##' 
##'  Ignored if \code{"fn"} is specified.
##'  
##' @param fn Custom function of the parameters to summarise against time.
##' This has optional first argument \code{t} representing time, and any remaining
##' arguments must be parameters of the distribution.  It should return a
##' vector of the same length as \code{t}.
##' 
##' 
##' @inheritParams summary.flexsurvreg
##' 
##' @export
##' 
summary.flexsurvrtrunc <- function(object, type="survival", fn=NULL,
                                t=NULL, quantiles=0.5, ci=TRUE, se=FALSE,
                                B=1000, cl=0.95,
                                ...)
{
    x <- object
    dat <- x$data
    type <- match.arg(type, c("survival","cumhaz","hazard","rmst","mean","median", "quantile","link"))
    start <- 0
    
    if(type == "mean"){
      if(!is.null(t)) warning("Mean selected, but time specified.  For restricted mean, set type to 'rmst'.")
      # Type = mean same as RMST w/ time = Inf
      t <- rep(Inf,length(start))
    }
    else if(type == "median"){
      if(!is.null(t)) warning("Median selected, but time specified.")
      t <- rep(0.5,length(start))
    }
    else if(type == "quantile"){
      t <- quantiles
      if((any(t<0) | any(t>1))){
        stop("Quantiles should not be less than 0 or greater than 1")
      }
      t <- rep(t,length(start))
    }
    else if(type == "rmst"){
        if (is.null(t))
            t <- max(dat$t)
    }
    else if (is.null(t))
        t <- sort(unique(dat$t))

    if (is.null(fn)) {
        fn <- summary_fns(x, type)
    }

    fn <- expand.summfn.args(fn)
    fncall <- list(t,start)
    dlist <- x$dlist

    basepars <- x$res[dlist$pars,"est"]
    names(basepars) <- dlist$pars
    basepars <- as.list(basepars)
    fncall[dlist$pars] <- basepars
    y <- do.call(fn, fncall)
    if (ci){
        res.ci <- cisumm.flexsurvreg.old(x, t, start, X=NULL, fn=fn, B=B, cl=cl)
        ly <- res.ci[,1]
        uy <-  res.ci[,2]
    }
    if (se){
        res.se <- sesumm.flexsurvreg.old(x, t, start, X=NULL, fn=fn, B=B)
    }
    if (type %in% c("median","mean"))
        ret <- data.frame(est=y, row.names=NULL)
    else if (type == "quantile")
        ret <- data.frame(quantile=t, est=y, row.names=NULL)
    else ret <- data.frame(time=t, est=y, row.names=NULL)
    if (ci) { ret$lcl <- ly; ret$ucl <- uy}
    if (se) { ret$se <- res.se }
    class(ret) <- c("summary.flexsurvrtrunc", class(ret))
    ret
}

cisumm.flexsurvreg.old <- function(x, t, start, X, fn, B=1000, cl=0.95) {
    if (all(is.na(x$cov)) || (B==0))
        ret <- array(NA, dim=c(length(t), 2))
    else {
        ret <- normbootfn.flexsurvreg(x=x, t=t, start=start, X=X, fn=fn, B=B)
        ret <- apply(ret, c(1,3), function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        ret <- t(ret[,1,])
    }
    ret
}

sesumm.flexsurvreg.old <- function(x, t, start, X, fn, B=1000) {
  if (all(is.na(x$cov)) || (B==0))
        ret <- numeric(length(t))
    else {
        ret <- normbootfn.flexsurvreg(x=x, t=t, start=start, X=X, fn=fn, B=B)
        ret <- apply(ret, c(1,3), sd, na.rm=TRUE)
        ret <- ret[1,]
    }
    ret
}
