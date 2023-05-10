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
  fit_surv <- predict(object, type = type)$.pred_time
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
##'  fitted cumulative hazard at the given observed data point, and other columns indicating the observation time,
##'  observed event status, and covariate values defining the data at this point.   
##'
##' The cumulative hazards \code{est} should form a censored sample from an Exponential(1).  
##' Therefore to check the fit of the model, plot a nonparametric estimate of the cumulative
##' hazard curve against a diagonal line through the origin, which is the theoretical cumulative
##' hazard trajectory of the Exponential(1).   
##' 
##' @examples
##'
##'   fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
##'   cs <- coxsnell_flexsurvreg(fitg)
##'   
##'   ## Model appears to fit well, with some small sample noise 
##'   surv <- survfit(Surv(cs$est, ovarian$fustat) ~ 1)
##'   plot(surv, fun="cumhaz")
##'   abline(0, 1, col="red")
##'   
##' @export
coxsnell_flexsurvreg <- function(x){
  mf <- model.frame(x, orig=TRUE)
  startstop <- "start" %in% colnames(mf[,1])
  tind <- if (startstop) "stop" else "time"
  t <- mf[,1][,tind]
  start <- if (startstop) mf[,1][,"start"] else 0  
  covnames <- attr(model.frame(x), "covnames")
  nd <- mf[,covnames,drop=FALSE]
  res <- summary(x, type="cumhaz", t=t, start=start, newdata=nd, cross=FALSE, 
                 ci=FALSE, se=FALSE, tidy=TRUE)
  res$status <- mf[,1][,"status"]
  res <- res[c("time","status", setdiff(names(res), c("time","status","est")), "est")]
  res    
}
