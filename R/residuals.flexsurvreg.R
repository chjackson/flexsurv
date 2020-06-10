c#' Calculate residuals for a flexible survival models
#'
#' Calculates residuals for \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}} model fits.
#'
#' @param object Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted survival model object.
#' @param type Character string for the type of residual desired. Currently only \code{"response"} is supported. More residual types may become available in future versions.
#'
#' @param ... Not currently used.
#'
#' @details Residuals of \code{type = "response"} are calculated as the naive difference between the observed survival and the covariate-specific predicted mean survival from \code{\link{predict.flexsurvreg}}, ignoring whether the event time is observed or censored. 
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
#' residuals(fitg)
#'
residuals.flexsurvreg <- function(object, type = "response", ...)
{
  type <- match.arg(type, "response")

  obs_surv <- unname(object$data$Y[, 1])
  fit_surv <- predict(object, type = type)$.pred

  res <- obs_surv - fit_surv
  res
}
