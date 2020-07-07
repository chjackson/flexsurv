##' Nonparametric estimator of survival from right-truncated, uncensored data
##'
##' Estimates the survivor function from right-truncated, uncensored data by reversing time, interpreting the data as left-truncated, applying the Kaplan-Meier / Lynden-Bell estimator and transforming back.
##'
##' Note that this does not estimate the untruncated survivor function - instead it estimates the survivor function truncated above at a time defined by the maximum possible time that might have been observed in the data. 
##'
##' @param t Vector of observed times from an initial event to a final event.
##'
##' @param rtrunc Individual-specific right truncation points, so that each individual's survival time \code{t} would not have been observed if it was greater than the corresponding element of \code{rtrunc}.
##' @param tmax Maximum possible time to event that could have been observed.
##'
##' @param data Data frame to find \code{t} and \code{rtrunc} in.  If not supplied, these should be in the working environment.
##'
##' @param conf.type Transformation used to compute a confidence interval from estimates and standard errors, as in \code{\link{survfit.formula}}. 
##'
##' @param conf.int Confidence level, defaulting to 0.95.
##'
##' @return A list with components:
##'
##' \code{time} Time points where the estimated survival changes.
##'
##' \code{surv} Estimated survival at \code{time}, truncated above at \code{tmax}.
##'
##' \code{std.err} Standard error of -log(survival).
##'
##' \code{lower} Lower confidence limits for survival.
##'
##' \code{upper} Upper confidence limits for survival.
##'
##' @details
##'
##' Define \eqn{X} as the time of the initial event, \eqn{Y} as the time of the final event, then we wish to determine the distribution of \eqn{T = Y- X}.
##'
##' Observations are only recorded if \eqn{Y \leq t_{max}}.  Then the distribution of \eqn{T} in the resulting sample is right-truncated by \code{rtrunc} \eqn{ = t_{max} - X}. 
##'
##' Equivalently, the distribution of \eqn{t_{max} - T} is left-truncated, since it is only observed if \eqn{t_{max} - T \geq X}.  Then the standard Kaplan-Meier type estimator as implemented in \code{\link[survival]{survfit}} is used (as described by Lynden-Bell, 1971) and the results transformed back.
##'
##' This situation might happen in a disease epidemic, where \eqn{X} is the date of disease onset for an individual, \eqn{Y} is the date of death, and we wish to estimate the distribution of the time \eqn{T} from onset to death, given we have only observed people who have died by the date \eqn{t_{max}}.
##'
##' If the estimated survival is unstable at the highest times, then consider replacing \code{tmax} by a slightly lower value, then if necessary, removing individuals with \code{t > tmax}, so that the estimand is changed to the survivor function truncated over a slightly narrower interval.
##'
##' @examples
##' 
##' ## simulate some event time data
##' set.seed(1) 
##' X <- rweibull(100, 2, 10)
##' T <- rweibull(100, 2, 10) 
##'
##' ## truncate above
##' tmax <- 20
##' obs <- X + T < tmax 
##' rtrunc <- tmax - X
##' dat <- data.frame(X, T, rtrunc)[obs,]
##' sf <-    survrtrunc(T, rtrunc, data=dat, tmax=tmax)
##' plot(sf, conf.int=TRUE)
##' ## Kaplan-Meier estimate ignoring truncation is biased 
##' sfnaive <- survfit(Surv(T) ~ 1, data=dat)
##' lines(sfnaive, conf.int=TRUE, lty=2, col="red")
##'
##' ## truncate above the maximum observed time 
##' tmax <- max(X + T) + 10
##' obs <- X + T < tmax 
##' rtrunc <- tmax - X
##' dat <- data.frame(X, T, rtrunc)[obs,]
##' sf <-    survrtrunc(T, rtrunc, data=dat, tmax=tmax)
##' plot(sf, conf.int=TRUE)
##' ## estimates identical to the standard Kaplan-Meier
##' sfnaive <- survfit(Surv(T) ~ 1, data=dat)
##' lines(sfnaive, conf.int=TRUE, lty=2, col="red")
##'
##' @references
##'
##' D. Lynden-Bell (1971)  A method of allowing for known observational selection in small samples applied to 3CR quasars. Monthly Notices of the Royal Astronomical Society, 155:95–118.
##' 
##' Seaman, S., Presanis, A. and Jackson, C. (2020) Review of methods for estimating distribution of time to event from right-truncated data.
##' 
##' @export
survrtrunc <- function(t, rtrunc, tmax, data=NULL, conf.type="log", conf.int=0.95){
    t <- eval(substitute(t), data, parent.frame())
    rtrunc <- eval(substitute(rtrunc), data, parent.frame())
    check_survrtrunc(t, rtrunc, tmax)
    X <- tmax - rtrunc
    trev <- tmax - t 
    event <- rep(1, length(X))
    sfleft <- survival::survfit(Surv(time=X, time2=trev, event=event, type="counting") ~ 1)
    sf <- list(time = rev(tmax - sfleft$time),
               surv = c(rev(1 - sfleft$surv)[-1], 0),
               std.err = sfleft$std.err)
    ci <- survfit_confint(p = sf$surv, se = sf$std.err,
                                     conf.type=conf.type, conf.int=conf.int)
    sf <- c(sf, ci)
    class(sf) <- "survrtrunc"
    sf
}

check_survrtrunc <- function(t, rtrunc, tmax) {
    if (!all(t <= rtrunc)) stop("Not all `t` are <= `rtrunc`")
    if (!all(rtrunc <= tmax)) stop("Not all `rtrunc` are <= `tmax`")
    if (!all(t <= tmax)) stop("Not all `t` are <= `tmax`")
}

##' Plot nonparametric estimates of survival from right-truncated data.
##'
##' \code{plot.survrtrunc} creates a new plot, while \code{lines.survrtrunc} adds lines to an exising plot.  
##' 
##' @param x Object of class \code{"survrtrunc"} as returned by \code{\link{survrtrunc}}. 
##' 
##' @param ... Other arguments to be passed to \code{\link{plot.survfit}} or \code{\link{lines.survfit}}. 
##' 
##' @rdname plot_survtrunc
##' @export
plot.survrtrunc <- function(x, ...){
    class(x) <- "survfit"
    plot(x, ...)
}

##' 
##' @rdname plot_survtrunc
##' @export
lines.survrtrunc <- function(x, ...){
    class(x) <- "survfit"
    lines(x, ...)
}


