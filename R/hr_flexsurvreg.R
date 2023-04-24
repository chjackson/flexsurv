##' Hazard ratio as a function of time from a parametric survival model
##' 
##' @inheritParams summary.flexsurvreg
##' 
##' @param x Object returned by \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}.
##' 
##' @param newdata A data frame with two rows, each specifying a set of covariate values.
##' The hazard ratio is calculated as hazard(z2)/hazard(z1), where z1 is the first row 
##' of \code{newdata} and z2 is the second row.
##'
##' \code{newdata} must be supplied unless the model \code{x} includes just one covariate.
##' With one covariate, a default is constructed, which defines the hazard ratio between
##' the second and first level of the factor (if the covariate is a factor), or between
##' a value of 1 and a value of 0 (if the covariate is numeric).
##' 
##' @return A data frame with estimate and confidence limits for the hazard ratio, and
##' one row for each of the times requested in \code{t}.
##'
##' @export
hr_flexsurvreg <- function(x, newdata=NULL, t=NULL, start=0, ci=TRUE, B=1000, cl=0.95, na.action=na.pass){
    if (is.null(newdata)){
        if (is.null(x$data))
            stop("`newdata` must be specified explicitly if the data have been removed from the model object")
        else 
            newdata <- hr_default_newdata(x)
    }
    if (is.null(newdata))
        stop("`newdata` must be specified unless the model has a single covariate")
  newdata <- as.data.frame(newdata)
  X <- newdata_to_X(x, newdata, na.action=na.action)
  args1 <- xt_to_fnargs(x, X[2,,drop=FALSE], t, start=start, type="hazard")
  args0 <- xt_to_fnargs(x, X[1,,drop=FALSE], t, start=start, type="hazard")
  fn <- expand.summfn.args(summary.fns(x, type="hazard"))
  est1 <- do.call(fn, args1)
  est0 <- do.call(fn, args0)
  
  tstart <- summfn_to_tstart(x=x, type="hr", t=t)
  t <- tstart$t; start <- tstart$start
  
  res <- data.frame(t = t) 
  res$est <- est1 / est0
  if (ci){
    argsboth <- xt_to_fnargs(x, X, t, start=start, type="hazard")
    sim <- normboot.flexsurvreg(x, B, X=attr(argsboth, "X"), tidy=TRUE)
    sim$"(time)" <- rep(rep(t, each=B), 2)
    sim$"(cov)" <- rep(c(1,2), each=B*length(t))
    sim0 <- sim[sim$"(cov)" == 1,,drop=FALSE]
    sim1 <- sim[sim$"(cov)" == 2,,drop=FALSE]
    argst <- list(t = rep(t, each=B), start = start)
    for (j in seq_along(x$aux))
      argst[[names(x$aux)[j]]] <- x$aux[[j]]
    args1 <- c(argst, as.list(as.data.frame(sim1))[x$dlist$pars])
    args0 <- c(argst, as.list(as.data.frame(sim0))[x$dlist$pars])
    ret <- do.call(fn, args1) / do.call(fn, args0)
    ret <- matrix(ret, nrow=B)
    retci <- apply(ret, 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
    res$lcl <- retci[1,]
    res$ucl <- retci[2,]
    # note it doesn't draw new simulations for each covariate value. OK 
  }
  res
}

hr_default_newdata <- function(x){
  mf <- model.frame(x)
  Xraw <- mf[,unique(attr(mf,"covnames.orig")),drop=FALSE]
  isfac <- sapply(Xraw, function(x){is.factor(x) || is.character(x)})
  isnum <- sapply(Xraw, function(x){is.numeric(x)})
  ncovs <- ncol(Xraw)
  if (ncovs == 1) {
      if (isfac){
          levs <- levels(Xraw[,1])
          if (length(levs)==1) stop("the covariate has only one one factor level")
      newdata <- list(levs[1:2])
    } else if (isnum) {
      newdata <- list(c(0,1))
    } else newdata <- NULL
    if (is.list(newdata))
       names(newdata) <- attr(mf, "covnames")
  } else newdata <- NULL
  newdata
}
