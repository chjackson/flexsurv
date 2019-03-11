##' Summaries of fitted flexible survival models
##' 
##' Return fitted survival, cumulative hazard or hazard at a series of times
##' from a fitted \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}
##' model.
##' 
##' Time-dependent covariates are not currently supported.  The covariate
##' values are assumed to be constant through time for each fitted curve.
##' 
##' @param object Output from \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvspline}}, representing a fitted survival model object.
##' @param newdata Data frame containing covariate values to produce fitted
##' values for.  Or a list that can be coerced to such a data frame.  There
##' must be a column for every covariate in the model formula, and one row for
##' every combination of covariates the fitted values are wanted for.  These
##' are in the same format as the original data, with factors as a single
##' variable, not 0/1 contrasts.
##' 
##' If this is omitted, if there are any continuous covariates, then a single
##' summary is provided with all covariates set to their mean values in the
##' data - for categorical covariates, the means of the 0/1 indicator variables
##' are taken.  If there are only factor covariates in the model, then all
##' distinct groups are used by default.
##' @param X Alternative way of defining covariate values to produce fitted
##' values for.  Since version 0.4, \code{newdata} is an easier way that
##' doesn't require the user to create factor contrasts, but \code{X} has been
##' kept for backwards compatibility.
##' 
##' Columns of \code{X} represent different covariates, and rows represent
##' multiple combinations of covariate values.  For example
##' \code{matrix(c(1,2),nrow=2)} if there is only one covariate in the model,
##' and we want survival for covariate values of 1 and 2.  A vector can also be
##' supplied if just one combination of covariates is needed.
##' 
##' For ``factor'' (categorical) covariates, the values of the contrasts
##' representing factor levels (as returned by the \code{\link{contrasts}}
##' function) should be used.  For example, for a covariate \code{agegroup}
##' specified as an unordered factor with levels \code{20-29, 30-39, 40-49,
##' 50-59}, and baseline level \code{20-29}, there are three contrasts.  To
##' return summaries for groups \code{20-29} and \code{40-49}, supply \code{X =
##' rbind(c(0,0,0), c(0,1,0))}, since all contrasts are zero for the baseline
##' level, and the second contrast is ``turned on'' for the third level
##' \code{40-49}.
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
##' \code{"link"} for the fitted value of the location parameter (i.e. the "linear predictor")
##' 
##' Ignored if \code{"fn"} is specified.
##' @param fn Custom function of the parameters to summarise against time.
##' This has optional first two arguments \code{t} representing time, and
##' \code{start} representing left-truncation points, and any remaining
##' arguments must be parameters of the distribution.  It should return a
##' vector of the same length as \code{t}.
##' @param t Times to calculate fitted values for. By default, these are the
##' sorted unique observation (including censoring) times in the data - for
##' left-truncated datasets these are the "stop" times.
##' @param quantiles If \code{type="quantile"}, this specifies the quantiles of the survival time distribution to return estimates for.
##' @param start Optional left-truncation time or times.  The returned
##' survival, hazard or cumulative hazard will be conditioned on survival up to
##' this time.
##' 
##' A vector of the same length as \code{t} can be supplied to allow different
##' truncation times for each prediction time, though this doesn't make sense
##' in the usual case where this function is used to calculate a predicted
##' trajectory for a single individual.  This is why the default \code{start}
##' time was changed for version 0.4 of \pkg{flexsurv} - this was previously a
##' vector of the start times observed in the data.
##' @param ci Set to \code{FALSE} to omit confidence intervals.
##' @param se Set to \code{TRUE} to include standard errors.
##' @param B Number of simulations from the normal asymptotic distribution of
##' the estimates used to calculate confidence intervals or standard errors.
##' Decrease for greater
##' speed at the expense of accuracy, or set \code{B=0} to turn off calculation
##' of CIs and SEs.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @param tidy If \code{TRUE}, then the results are returned as a tidy data
##' frame instead of a list.  This can help with using the \pkg{ggplot2}
##' package to compare summaries for different covariate values.
##' @param ... Further arguments passed to or from other methods.  Currently
##' unused.
##' @return If \code{tidy=FALSE}, a list with one component for each unique
##' covariate value (if there are only categorical covariates) or one component
##' (if there are no covariates or any continuous covariates).  Each of these
##' components is a matrix with one row for each time in \code{t}, giving the
##' estimated survival (or cumulative hazard, or hazard) and 95\% confidence
##' limits.  These list components are named with the covariate names and
##' values which define them.
##' 
##' If \code{tidy=TRUE}, a data frame is returned instead.  This is formed by
##' stacking the above list components, with additional columns to identify the
##' covariate values that each block corresponds to.
##' 
##' If there are multiple summaries, an additional list component named
##' \code{X} contains a matrix with the exact values of contrasts (dummy
##' covariates) defining each summary.
##' 
##' The \code{\link{plot.flexsurvreg}} function can be used to quickly plot
##' these model-based summaries against empirical summaries such as
##' Kaplan-Meier curves, to diagnose model fit.
##' 
##' Confidence intervals are obtained by sampling randomly from the asymptotic
##' normal distribution of the maximum likelihood estimates and then taking quantiles
##' (see, e.g. Mandel (2013)).
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{flexsurvreg}}, \code{\link{flexsurvspline}}.
##' @references Mandel, M. (2013). "Simulation based confidence intervals for
##' functions with complicated derivatives." The American Statistician (in
##' press).
##' @keywords models
##' @export
summary.flexsurvreg <- function(object, newdata=NULL, X=NULL, type="survival", fn=NULL,
                                t=NULL, quantiles=0.5, start=0, ci=TRUE, se=FALSE,
                                B=1000, cl=0.95, tidy=FALSE,
                                ...)
{
    x <- object
    dat <- x$data
    Xraw <- model.frame(x)[,unique(attr(model.frame(x),"covnames.orig")),drop=FALSE]
    isfac <- sapply(Xraw, function(x){is.factor(x) || is.character(x)})
    type <- match.arg(type, c("survival","cumhaz","hazard","rmst","mean","median", "quantile","link"))
    if (is.null(newdata)){
        if (is.vector(X)) X <- matrix(X, nrow=1)
        if (x$ncovs > 0 && is.null(X)) {
            ## if any continuous covariates, calculate fitted survival for "average" covariate value
            if (!all(isfac)){
                nd <- colMeans(model.matrix(x))
                X <- matrix(nd ,nrow=1, dimnames=list(NULL,names(nd)))
                attr(X, "newdata") <- as.data.frame(X)
            }
            ## else calculate for all different factor groupings
            else {
                X <- unique(model.matrix(x))
                ## build names like "COVA=value1,COVB=value2"
                nam <- as.matrix(unique(Xraw))
                for (i in 1:ncol(nam)) nam[,i] <- paste(colnames(nam)[i], nam[,i], sep="=")
                rownames(X) <- apply(nam, 1, paste, collapse=",")
                attr(X, "newdata") <- unique(Xraw)
            }
        }
        else if (is.null(X)) X <- as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1))
        else if (!is.matrix(X) || (is.matrix(X) && ncol(X) != x$ncoveffs)) {
            plural <- if (x$ncoveffs > 1) "s" else ""
            stop("expected X to be a matrix with ", x$ncoveffs, " column", plural, " or a vector with ", x$ncoveffs, " element", plural)
        }
        else {
            attr(X, "newdata") <- X
            colnames(attr(X, "newdata")) <- colnames(model.matrix(x))
        }
    } else
        X <- form.model.matrix(object, as.data.frame(newdata))
    
    if(type == "mean"){
      if(!is.null(t)) warning("Mean selected, but time specified.  For restricted mean, set type to 'rmst'.")
      # Type = mean same as RMST w/ time = Inf
      t <- rep(Inf,length(start))
    }
    else if(type == "median"){
      if(!is.null(t)) warning("Median selected, but time specified.")
      t <- rep(0.5,length(start))
    }
    else if(type == "link"){
      if(!is.null(t)) warning("`link` selected, but time specified.")
      t <- rep(0,length(start))
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
            t <- max(dat$Y[,"time1"])
    }
    else if (is.null(t))
        t <- sort(unique(dat$Y[,"stop"]))
    if (length(start)==1)
        start <- rep(start, length(t))
    else if (length(start) != length(t))
        stop("length of \"start\" is ",length(start),". Should be 1, or length of \"t\" which is ",length(t))

    if (is.null(fn)) {
        fn <- summary.fns(x, type)
    }
    fn <- expand.summfn.args(fn)
    fncall <- list(t,start)
    beta <- if (x$ncovs==0) 0 else x$res[x$covpars,"est"]
    dlist <- x$dlist
    ret <- vector(nrow(X), mode="list")
    if(!is.null(newdata)){
        nd <- attr(X, "newdata")
        covnames <- apply(as.data.frame(nd), 1, function(x)paste0(names(nd), "=", x, collapse=", "))
    }
    else covnames <- rownames(X)
    names(ret) <- covnames
    for (i in 1:nrow(X)) {
        basepars.mat <- add.covs(x, x$res.t[dlist$pars,"est"], beta, X[i,,drop=FALSE], transform=FALSE)
        basepars <- as.list(as.data.frame(basepars.mat))
        fncall[dlist$pars] <- basepars
        if (type=="link")
            x$aux$location <- x$dlist$location
        for (j in seq_along(x$aux)){
            fncall[[names(x$aux)[j]]] <- x$aux[[j]]
        }
        y <- do.call(fn, fncall)
        if (ci){
            res.ci <- cisumm.flexsurvreg(x, t, start, X[i,,drop=FALSE], fn=fn, B=B, cl=cl)
            ly <- res.ci[,1]
            uy <-  res.ci[,2]
        }
        if (se){
            res.se <- sesumm.flexsurvreg(x, t, start, X[i,,drop=FALSE], fn=fn, B=B)
        }
        if (type %in% c("median","mean"))
            ret[[i]] <- data.frame(est=y, row.names=NULL)
        else if (type == "quantile")
            ret[[i]] <- data.frame(quantile=t, est=y, row.names=NULL)
        else ret[[i]] <- data.frame(time=t, est=y, row.names=NULL)
        if (ci) { ret[[i]]$lcl <- ly; ret[[i]]$ucl <- uy}
        if (se) { ret[[i]]$se <- res.se }
    }
    if (x$ncovs>0) attr(ret,"X") <- X
    if (tidy) {
        ret <- do.call("rbind", ret)
        if (x$ncovs>0) {
            nd <- attr(X, "newdata")
            covdf <- nd[rep(seq_len(nrow(nd)), each=length(t)), , drop=FALSE]
            rownames(ret) <- rownames(covdf) <- NULL
            ret <- cbind(ret, covdf)
        }
    }
    class(ret) <- c("summary.flexsurvreg",class(ret))
    ret
}

summary.fns <- function(x, type){
    switch(type,   # TODO warn for clashing arguments in dfns
           "survival" = function(t,start,...) {
               ret <- (1 - x$dfns$p(t,...))/(1 - x$dfns$p(start,...))
               ret[t<start] <- 1 # prob[t<start] was previously 0
               ret
           },
           "median" = function(start,...) {
             start_p = 1 - x$dfns$p(start,...)
             med_from_start = start_p/2
             ret = x$dfns$q(med_from_start,...)
           },
           "quantile" = function(t=0.5, start,...) {
             start_p = 1 - x$dfns$p(start,...)
             med_from_start = start_p * t
             ret = x$dfns$q(med_from_start,...)
           },
           "hazard" = function(t,start,...) {
               ret <- x$dfns$h(t,...) * (1 - x$dfns$p(start,...))
               ret[t<start] <- 0
               ret
           },
           "cumhaz" = function(t,start,...) {
               ret <- x$dfns$H(t,...) - x$dfns$H(start,...)
               ret[t<start] <- 0
               ret
           },
           "rmst" = function(t,start,...) x$dfns$rmst(t,start=start, ...),
           "mean" = function(t,start,...) x$dfns$mean(...),
           "link" = function(...){
               args <- list(...)
               args[[args$location]]
           }
    )
}

##' @export
print.summary.flexsurvreg <- function(x, ...){
    if (!inherits(x, "data.frame")){
        for (i in seq_along(x)){
            cat(names(x)[i], "\n")
            print(x[[i]])
            if (i<length(x)) cat("\n")
        }
    } else print.data.frame(x)
}

## TODO would converting newdata to X be better handled in this function

add.covs <- function(x, pars, beta, X, transform=FALSE){  ## TODO option to transform on input
    nres <- nrow(X)
    if (!is.matrix(pars)) pars <- matrix(pars, nrow=nres, ncol=length(pars), byrow=TRUE)
    if (!is.matrix(beta)) beta <- matrix(beta, nrow=1)
    for (j in seq(along=x$dlist$pars)){
        covinds <- x$mx[[x$dlist$pars[j]]]
        if (length(covinds) > 0){
            pars[,j] <- pars[,j] + beta[,covinds] %*% t(X[,covinds,drop=FALSE])
        }
        if (!transform)
            pars[,j] <- x$dlist$inv.transforms[[j]](pars[,j])
    }
    colnames(pars) <- x$dlist$pars
    pars
}

## Draw B samples from multivariate normal distribution of baseline
## parameter estimators, for given covariate values



##' Simulate from the asymptotic normal distribution of parameter estimates.
##' 
##' Produce a matrix of alternative parameter estimates under sampling
##' uncertainty, at covariate values supplied by the user.  Used by
##' \code{\link{summary.flexsurvreg}} for obtaining confidence intervals around
##' functions of parameters.
##' 
##' 
##' @param x A fitted model from \code{\link{flexsurvreg}} (or \code{\link{flexsurvspline}}).
##' @param B Number of samples.
##' @param newdata Data frame or list containing the covariate values to
##' evaluate the parameters at.  If there are covariates in the model, at least
##' one of \code{newdata} or \code{X} must be supplied, unless \code{raw=TRUE}.
##' @param X Alternative (less convenient) format for covariate values: a
##' matrix with one row, with one column for each covariate or factor contrast.
##' Formed from all the "model matrices", one for each named parameter of the
##' distribution, with intercepts excluded, \code{cbind}ed together.
##' @param transform \code{TRUE} if the results should be transformed to the
##' real-line scale, typically by log if the parameter is defined as positive.
##' The default \code{FALSE} returns parameters on the natural scale.
##' @param raw Return samples of the baseline parameters and the covariate
##' effects, rather than the default of adjusting the baseline parameters for
##' covariates.
##' @return If \code{newdata} includes only one covariate combination, a matrix
##' will be returned with \code{B} rows, and one column for each named
##' parameter of the survival distribution.
##' 
##' If more than one covariate combination is requested (e.g. \code{newdata} is
##' a data frame with more than one row), then a list of matrices will be
##' returned, one for each covariate combination.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{summary.flexsurvreg}}
##' @references Mandel, M. (2013). "Simulation based confidence intervals for
##' functions with complicated derivatives." The American Statistician (in
##' press).
##' @keywords models
##' @examples
##' 
##'     fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
##'     normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
##'     normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
##'     normboot.flexsurvreg(fite, B=10, newdata=list(age=0))  ## closer to...
##'     fite$res
##' @export
normboot.flexsurvreg <- function(x, B, newdata=NULL, X=NULL, transform=FALSE, raw=FALSE){
    if (x$ncovs > 0 && !raw) {
        if (is.null(X)) {
            if (is.null(newdata)) stop("neither \"newdata\" nor \"X\" supplied")
            X <- form.model.matrix(x, as.data.frame(newdata))
        }
    } else X <- as.matrix(0, nrow=1, ncol=1)
    sim <- matrix(nrow=B, ncol=nrow(x$res))
    colnames(sim) <- rownames(x$res)
    if (is.na(x$cov[1])) stop("Covariance matrix not available from non-converged model")
    sim[,x$optpars] <- rmvnorm(B, x$opt$par, x$cov)
    sim[,x$fixedpars] <- rep(x$res.t[x$fixedpars,"est"],each=B)
    if (x$ncovs > 0 && !raw){
        beta <- sim[, x$covpars, drop=FALSE]
        if (nrow(X)==1){
            res <- sim[,x$dlist$pars,drop=FALSE]
            res <- add.covs(x=x, pars=res, beta=beta, X=X, transform=transform)
        }  else {
            res <- vector(nrow(X), mode="list")
            for (i in 1:nrow(X)) {
                res[[i]] <- sim[,x$dlist$pars,drop=FALSE]
                res[[i]] <- add.covs(x=x, pars=res[[i]], beta=beta, X=X[i,,drop=FALSE], transform=transform)
            }
        }
    } else {
        res <- sim
        if (!transform){
            for (j in seq(along=x$dlist$pars)){
                res[,j] <- x$dlist$inv.transforms[[j]](res[,j])
            }
        }
    }
    attr(res, "X") <- X
    res
}

### Compute CIs for survival, cumulative hazard, hazard, or user
### defined function, at supplied times t and covariates X, using
### random sample of size B from the assumed MVN distribution of MLEs.

normbootfn.flexsurvreg <- function(x, t, start, newdata=NULL, X=NULL, fn, B){
    sim <- normboot.flexsurvreg(x, B, newdata=newdata, X=X)
    X <- attr(sim, "X")
    if (!is.list(sim)) sim <- list(sim)
    ret <- array(NA_real_, dim=c(nrow(X), B, length(t)))
    for (k in 1:nrow(X)){
        for (i in seq(length=B)) {
            fncall <- list(t,start)
            for (j in seq(along=x$dlist$pars))
                fncall[[x$dlist$pars[j]]] <- sim[[k]][i,j]
            for (j in seq_along(x$aux))
                fncall[[names(x$aux)[j]]] <- x$aux[[j]]
            ret[k,i,] <- do.call(fn, fncall)
        }
    }
    if (nrow(X)==1) ret[1,,,drop=FALSE] else ret
}

cisumm.flexsurvreg <- function(x, t, start, X, fn, B=1000, cl=0.95) {
    if (any(is.na(x$res[,2])) || (B==0))
        ret <- array(NA, dim=c(length(t), 2))
    else {
        ret <- normbootfn.flexsurvreg(x=x, t=t, start=start, X=X, fn=fn, B=B)
        ret <- apply(ret, c(1,3), function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        ret <- t(ret[,1,])
    }
    ret
}

sesumm.flexsurvreg <- function(x, t, start, X, fn, B=1000) {
    if (any(is.na(x$res[,2])) || (B==0))
        ret <- numeric(length(t))
    else {
        ret <- normbootfn.flexsurvreg(x=x, t=t, start=start, X=X, fn=fn, B=B)
        ret <- apply(ret, c(1,3), sd, na.rm=TRUE)
        ret <- ret[1,]
    }
    ret
}
