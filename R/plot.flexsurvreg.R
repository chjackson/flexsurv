##' Plots of fitted flexible survival models
##' 
##' Plot fitted survival, cumulative hazard or hazard from a parametric model
##' against nonparametric estimates to diagnose goodness-of-fit.  Alternatively
##' plot a user-defined function of the model parameters against time.
##' 
##' 
##' @param x Output from \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvspline}}, representing a fitted survival model object.
##' @param newdata Data frame containing covariate values to produce fitted
##' values for.  See \code{\link{summary.flexsurvreg}}.
##' 
##' If there are only factor covariates in the model, then Kaplan-Meier (or
##' nonparametric hazard...)  curves are plotted for all distinct groups, and
##' by default, fitted curves are also plotted for these groups.  To plot
##' Kaplan-Meier and fitted curves for only a subset of groups, use
##' \code{plot(survfit())} followed by \code{lines.flexsurvreg()}.
##' 
##' If there are any continuous covariates, then a single population
##' Kaplan-Meier curve is drawn. By default, a single fitted curve is drawn
##' with the covariates set to their mean values in the data - for categorical
##' covariates, the means of the 0/1 indicator variables are taken.
##' @param X Alternative way to supply covariate values, as a model matrix.
##' See \code{\link{summary.flexsurvreg}}.  \code{newdata} is an easier way.
##' @param type \code{"survival"} for survival, to be plotted against
##' Kaplan-Meier estimates from \code{\link[survival]{plot.survfit}}.
##' 
##' \code{"cumhaz"} for cumulative hazard, plotted against transformed
##' Kaplan-Meier estimates from \code{\link[survival]{plot.survfit}}.
##' 
##' \code{"hazard"} for hazard, to be plotted against smooth nonparametric
##' estimates from \code{\link[muhaz]{muhaz}}.  The nonparametric estimates
##' tend to be unstable, and these plots are intended just to roughly indicate
##' the shape of the hazards through time.  The \code{min.time} and
##' \code{max.time} options to \code{\link[muhaz]{muhaz}} may sometimes need to
##' be passed as arguments to \code{\link{plot.flexsurvreg}} to avoid an error
##' here.
##' 
##' Ignored if \code{"fn"} is specified.
##' @param fn Custom function of the parameters to summarise against time.  The
##' first two arguments of the function must be \code{t} representing time, and
##' \code{start} representing left-truncation points, and any remaining
##' arguments must be parameters of the distribution.  It should return a
##' vector of the same length as \code{t}.
##' @param t Vector of times to plot fitted values for, see
##' \code{\link{summary.flexsurvreg}}.
##' @param start Left-truncation points, see \code{\link{summary.flexsurvreg}}.
##' @param est Plot fitted curves (\code{TRUE} or \code{FALSE}.)
##' @param ci Plot confidence intervals for fitted curves. By default, this is
##' \code{TRUE} if one observed/fitted curve is plotted, and \code{FALSE} if
##' multiple curves are plotted.
##' @param B Number of simulations controlling accuracy of confidence
##' intervals, as used in \code{\link[=summary.flexsurvreg]{summary}}.
##' Decrease for greater speed at the expense of accuracy, or set \code{B=0} to
##' turn off calculation of CIs.
##' @param cl Width of confidence intervals, by default 0.95 for 95\%
##' intervals.
##' @param col.obs Colour of the nonparametric curve.
##' @param lty.obs Line type of the nonparametric curve.
##' @param lwd.obs Line width of the nonparametric curve.
##' @param col Colour of the fitted parametric curve(s).
##' @param lty Line type of the fitted parametric curve(s).
##' @param lwd Line width of the fitted parametric curve(s).
##' @param col.ci Colour of the fitted confidence limits, defaulting to the
##' same as for the fitted curve.
##' @param lty.ci Line type of the fitted confidence limits.
##' @param lwd.ci Line width of the fitted confidence limits.
##' @param ylim y-axis limits: vector of two elements.
##' @param add If \code{TRUE}, add lines to an existing plot, otherwise new
##' axes are drawn.
##' @param ... Other options to be passed to \code{\link{plot.survfit}} or
##' \code{\link[muhaz]{muhaz}}, for example, to control the smoothness of the
##' nonparametric hazard estimates.  The \code{min.time} and \code{max.time}
##' options to \code{\link[muhaz]{muhaz}} may sometimes need to be changed from
##' the defaults.
##' @note Some standard plot arguments such as \code{"xlim","xlab"} may not
##' work.  This function was designed as a quick check of model fit.  Users
##' wanting publication-quality graphs are advised to set up an empty plot with
##' the desired axes first (e.g. with \code{plot(...,type="n",...)}), then use
##' suitable \code{\link{lines}} functions to add lines.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{flexsurvreg}}
##' @keywords models hplot
##' @export
plot.flexsurvreg <- function(x, newdata=NULL, X=NULL, type="survival", fn=NULL, t=NULL, start=0,
                             est=TRUE, ci=NULL, B=1000, cl=0.95,
                             col.obs="black", lty.obs=1, lwd.obs=1,
                             col="red",lty=1,lwd=2,
                             col.ci=NULL,lty.ci=2,lwd.ci=1,ylim=NULL,
                             add=FALSE,...)
{
    ## don't calculate or plot CIs by default if all covs are categorical -> multiple curves
    mf <- model.frame(x)
    Xraw <- mf[,attr(mf, "covnames.orig"), drop=FALSE]
    if (is.null(ci))
        ci <- ((x$ncovs == 0) || (!(sapply(Xraw,is.factor))))
    if (!ci) B <- 0
    summ <- summary(x, newdata=newdata, X=X, type=type, fn=fn, t=t, start=start, ci=ci, B=B, cl=cl)
    t <- summ[[1]]$time
    X <- if (is.null(attr(summ,"X"))) as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1)) else attr(summ,"X")
    if (is.null(col.ci)) col.ci <- col
    if (is.null(lwd.ci)) lwd.ci <- lwd
    dat <- x$data
    isfac <- sapply(Xraw,is.factor)
    if (!is.null(fn)) type <- ""
    if (!add) {
        mm <- as.data.frame(model.matrix(x))
        form <- "Surv(dat$Y[,\"start\"],dat$Y[,\"stop\"],dat$Y[,\"status\"]) ~ "
        form <- paste(form, if (x$ncovs > 0 && all(isfac)) paste("mm[,",1:x$ncoveffs,"]", collapse=" + ") else 1)
        form <- as.formula(form)
        ## If any continuous covariates, it is hard to define subgroups
        ## so just plot the population survival
        if (type=="survival") {
            plot(survfit(form, data=mm), col=col.obs, lty=lty.obs, lwd=lwd.obs, ylim=ylim, ...)
        }
        else if (type=="cumhaz") {
            plot(survfit(form, data=mm), fun="cumhaz", col=col.obs, lty=lty.obs, lwd=lwd.obs, ylim=ylim, ...)
        }
        else if (type=="hazard") {
            muhaz.args <- list(...)[names(list(...)) %in% names(formals(muhaz))]
            if (is.null(muhaz.args$min.time)) muhaz.args$min.time <- 0
            if (is.null(muhaz.args$max.time)) muhaz.args$max.time <- with(as.data.frame(dat$Y), max(time[status==1]))
            plot.args <- list(...)[!names(list(...)) %in% names(formals(muhaz))]
            if (!all(dat$Y[,"start"]==0)) warning("Left-truncated data not supported by muhaz: ignoring truncation point when plotting observed hazard")
            if (any(dat$Y[,"status"] > 1)) stop("Interval-censored data not supported by muhaz")
            if (!all(isfac)){
                haz <- do.call("muhaz", c(list(times=dat$Y[,"stop"], delta=dat$Y[,"status"]), muhaz.args))
                do.call("plot", c(list(haz), list(col=col.obs, lty=lty.obs, lwd=lwd.obs), plot.args))
            }
            else {
                ## plot hazard for all groups defined by unique combinations of covariates
                group <- if(x$ncovs>0) do.call("interaction", mm) else factor(rep(0,nrow(dat$Y)))
                Xgroup <- factor(do.call("interaction", as.data.frame(X)), levels=levels(group))
                haz <- list()
                for (i in 1:nrow(X)) {
                    haz[[i]] <- do.call("muhaz", c(list(times=dat$Y[,"time"], delta=dat$Y[,"status"], subset=(group==Xgroup[i])), muhaz.args))
                }
                if (missing(ylim))
                    ylim <- range(sapply(haz, function(x)range(x$haz.est)))
                do.call("plot", c(list(haz[[1]]), list(col=col.obs, lty=lty.obs, lwd=lwd.obs, ylim=ylim), plot.args))
                if (nrow(X)>1) {
                    for (i in 1:nrow(X)) {
                        lines(haz[[i]], col=col.obs, lty=lty.obs, lwd=lwd.obs)
                    }
                }
            }
        }
    }
    col <- rep(col, length=nrow(X)); lty=rep(lty, length=nrow(X)); lwd=rep(lwd, length=nrow(X))
    col.ci <- rep(col.ci, length=nrow(X)); lty.ci=rep(lty.ci, length=nrow(X)); lwd.ci=rep(lwd.ci, length=nrow(X))
    for (i in 1:nrow(X)) {
        if (est) lines(summ[[i]]$time, summ[[i]]$est, col=col[i], lty=lty[i], lwd=lwd[i])
        if (ci) {
            lines(summ[[i]]$time, summ[[i]]$lcl, col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i])
            lines(summ[[i]]$time, summ[[i]]$ucl, col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i])
        }
    }
}



##' Add fitted flexible survival curves to a plot
##' 
##' Add fitted survival (or hazard or cumulative hazard) curves from a
##' \code{\link{flexsurvreg}} model fit to an existing plot.
##' 
##' Equivalent to \code{\link{plot.flexsurvreg}(...,add=TRUE)}.
##' 
##' 
##' @param x Output from \code{\link{flexsurvreg}}, representing a fitted
##' survival model object.
##' @param newdata Covariate values to produce fitted curves for, as a data
##' frame, as described in \code{\link{plot.flexsurvreg}}.
##' @param X Covariate values to produce fitted curves for, as a matrix, as
##' described in \code{\link{plot.flexsurvreg}}.
##' @param type \code{"survival"} for survival, \code{"cumhaz"} for cumulative
##' hazard, or \code{"hazard"} for hazard, as in
##' \code{\link{plot.flexsurvreg}}.
##' @param t Vector of times to plot fitted values for.
##' @param est Plot fitted curves (\code{TRUE} or \code{FALSE}.)
##' @param ci Plot confidence intervals for fitted curves.
##' @param B Number of simulations controlling accuracy of confidence
##' intervals, as used in \code{\link[=summary.flexsurvreg]{summary}}.
##' @param cl Width of confidence intervals, by default 0.95 for 95\%
##' intervals.
##' @param col Colour of the fitted curve(s).
##' @param lty Line type of the fitted curve(s).
##' @param lwd Line width of the fitted curve(s).
##' @param col.ci Colour of the confidence limits, defaulting to the same as
##' for the fitted curve.
##' @param lty.ci Line type of the confidence limits.
##' @param lwd.ci Line width of the confidence limits, defaulting to the same
##' as for the fitted curve.
##' @param ... Other arguments to be passed to the generic \code{\link{plot}}
##' and \code{\link{lines}} functions.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{flexsurvreg}}
##' @keywords models aplot
##' @export
lines.flexsurvreg <- function(x, newdata=NULL, X=NULL, type="survival", t=NULL,
                              est=TRUE, ci=NULL, B=1000, cl=0.95,
                              col="red",lty=1,lwd=2,
                              col.ci=NULL,lty.ci=2,lwd.ci=1, ...)
{
    plot.flexsurvreg(x, newdata=newdata, X=X, type=type, t=t, est=est, ci=ci, B=B, cl=cl,
                     col=col, lty=lty, lwd=lwd, col.ci=col.ci,lty.ci=lty.ci,lwd.ci=lwd.ci, add=TRUE, ...)
}
