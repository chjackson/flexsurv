#' Calculate residuals for flexible survival models
#'
#' Calculates residuals for \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}} model fits.
#'
#' @param object Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted survival model object.
#' 
#' @param type Character string for the type of residual desired. Currently only \code{"response"} and \code{"coxsnell"} are supported. More residual types may become available in future versions.
#'
#' @param ... Not currently used.
#'
#' @details Residuals of \code{type = "response"} are calculated as the naive difference between the observed survival and the covariate-specific predicted mean survival from \code{\link{predict.flexsurvreg}}, ignoring whether the event time is observed or censored.
#'
#' \code{type="coxsnell"} returns the Cox-Snell residual, defined as the estimated cumulative hazard at each data point.  To check the fit of the
#' A more fully featured utility for this is provided in the function \code{\link{coxsnell_flexsurvreg}}. 
#'
#' @return Numeric vector with the same length as \code{nobs(object)}.
#'
#' @seealso \code{\link{predict.flexsurvreg}}
#'
#' @importFrom stats residuals
#'
#' @export
#'
#' @examples
#'
#' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
#' residuals(fitg, type="response")
#'
#' 
#'
residuals.flexsurvreg <- function(object, type = "response", ...)
{
  type <- match.arg(type, c("response","coxsnell"))

  if (type=="response"){ 
  obs_surv <- unname(object$data$Y[, 1])
  fit_surv <- predict(object, type = type)$.pred
  res <- obs_surv - fit_surv
  } else if (type=="coxsnell") {
      cx <- coxsnell_flexsurvreg(object)
      res <- cx$est
  }

  res
}

##' Cox-Snell residuals from a parametric survival model
##' 
##' @param x Object returned by \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}} representing a fitted survival model
##'
##' @return A data frame with a column called \code{est} giving the Cox-Snell residual, defined as the fitted cumulative hazard at each data point.
##'  fitted cumulative hazard at the given observed data point, and other columns indicating the observation time
##'  and covariate values defining the data at this point.   
##'
##' An extra column \code{"(qexp)"} gives the equally-spaced quantiles of a standard 
##' exponential distribution in the same order as \code{est}.   To check the fit of the model, 
##' \code{"(qexp)"} is plotted against \code{est}, and the points should form a straight line 
##' through the origin with slope 1. 
##' 
##' @examples
##'
##'   fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
##'   cs <- coxsnell_flexsurvreg(fitg)
##'   
##'   ## Model doesn't appear to fit well since the cumulative hazards are underestimated.
##'   ## In this example, this is probably because the dataset is small, 
##'   ## hence the point estimate is noisy.
##'   plot(cs$"(qexp)", cs$est, pch=19, xlab="Theoretical quantiles", ylab="Cumulative hazard")
##'   abline(a=0,b=1,col="red",lwd=2)
##'   
##'   ## Alternative way to produce the same plot using "qqplot"
##'   qy <- qexp(ppoints(nrow(cs),0))
##'   qqplot(qy, cs$est)
##'   abline(a=0,b=1, col="red", lwd=2)
##'   
##'   ## A log transform may or may not bring out the pattern more clearly
##'   plot(log(cs$"(qexp)"), log(cs$est), pch=19)
##'   abline(a=0,b=1, col="red", lwd=2)
##'   
##'   ## In the model `fitg`, the fitted cumulative hazard is lower than the true cumulative hazard
##'   ## Another way to show this is to compare parametric vs nonparametric estimates of 
##'   ## the cumulative hazard 
##'   plot(fitg, type="cumhaz", ci=FALSE)
##'   
##'   ## Alternative example: where the true model is fitted to simulated data
##'   ## The model fits well
##'   y <- rweibull(10000, 2, 2)
##'   fite <- flexsurvreg(Surv(y) ~ 1, dist="weibull")
##'   cs <- coxsnell_flexsurvreg(fite)
##'   plot(cs$"(qexp)", cs$est, pch=19, xlab="Theoretical quantiles", ylab="Cumulative hazard")
##'   abline(a=0,b=1,col="red",lwd=2)
##'   
##' @export
coxsnell_flexsurvreg <- function(x){
    mf <- model.frame(x, orig=TRUE)
    t <- mf[,1][,"time"]
    covnames <- attr(model.frame(x), "covnames")
    nd <- mf[,covnames,drop=FALSE]
    res <- summary(x, type="cumhaz", t=t, newdata=nd, cross=FALSE, 
                   ci=FALSE, se=FALSE, tidy=TRUE)
    res$"(qexp)" <- NA
    res$"(qexp)"[order(res$est)] <- qexp(ppoints(nrow(res), 0))
    res    
}
