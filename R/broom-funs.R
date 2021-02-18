#' Tidy a flexsurv model object
#'
#' Tidy summarizes information about the components of the model into a tidy data frame.
#'
#' @param x Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted survival model object.
#'
#' @param conf.int Logical. Should confidence intervals be returned? Default is \code{FALSE}.
#'
#' @param conf.level The confidence level to use for the confidence interval if \code{conf.int = TRUE}. Default is \code{0.95}.
#'
#' @param pars Character vector for one of \code{"all"}, \code{"coefs"}, or \code{"baseline"} for all parameters, covariate effects (i.e. regression betas), or baseline distribution paramaters, respectively. Default is \code{"all"}.
#'
#' @param transform Character vector of transformations to apply to requested \code{pars}. Default is \code{"none"}, which returns \code{pars} as-is.
#'
#' Users can specify one or both types of transformations:
#'
#' * \code{"baseline.real"} which transforms the baseline distribution parameters to the real number line used for estimation.
#'
#' * \code{"coefs.exp"} which exponentiates the covariate effects.
#'
#' See \code{Details} for a more complete explanation.
#'
#' @param ... Not currently used.
#'
#' @details \code{flexsurvreg} models estimate two types of coefficients, baseline distribution parameters, and covariate effects which act on the baseline distribution. By design, \code{flexsurvreg} returns distribution parameters on the same scale as is found in the relevant \code{d/p/q/r} functions. Covariate effects are returned on the log-scale, which represents either log-time ratios (accelerated failure time models) or log-hazard ratios for proportional hazard models. By default, \code{tidy()} will return baseline distribution parameters on their natural scale and covariate effects on the log-scale.
#'
#' To transform the baseline distribution parameters to the real-value number line (the scale used for estimation), pass the character argument \code{"baseline.real"} to \code{transform}. To get time ratios or hazard ratios, pass \code{"coefs.exp"} to \code{transform}. These transformations may be done together by submitting both arguments as a character vector.
#'
#' @return A \code{\link{tibble}} containing the columns: \code{term}, \code{estimate}, \code{std.error}, \code{statistic}, \code{p.value}, \code{conf.low}, and \code{conf.high}, by default.
#'
#' \code{statistic} and \code{p.value} are only provided for covariate effects (\code{NA} for baseline distribution parameters). These are computed as Wald-type test statistics with p-values from a standard normal distribution.
#'
#' @importFrom purrr map2_dbl
#'
#' @md
#'
#' @export
#'
#' @examples
#'
#' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
#' tidy(fitg)
#' tidy(fitg, pars = "coefs", transform = "coefs.exp")
#'
tidy.flexsurvreg <- function(x, conf.int = FALSE, conf.level = 0.95,
                             pars = "all", transform = "none", ...)
{
  assertthat::assert_that(is.logical(conf.int))

  if (conf.int) assertthat::assert_that(is.numeric(conf.level),
                                        conf.level > 0, conf.level < 1,
                                        length(conf.level) == 1,
                                        msg = "`conf.level` must be length one and between 0 and 1")

  pars <- match.arg(pars, c("all", "coefs", "baseline"))
  transform <- match.arg(transform, c("none", "baseline.real", "coefs.exp"),
                         several.ok = TRUE)

  dist_pars <- x$dlist$pars

  if ("baseline.real" %in% transform) vals <- x$res.t else vals <- x$res

  coefs <- vals[, "est"]
  terms <- names(coefs)
  ses <- vals[, "se"]
  stats <- coefs / ses; stats[names(stats) %in% dist_pars] <- NA
  pvals <- pnorm(abs(stats), lower.tail = FALSE)

  coef_pars <- terms[which(!terms %in% dist_pars)]

  if (conf.int) cis <- confint(x, level = conf.level) else cis <- NULL

  if (!"baseline.real" %in% transform) {
    cis[rownames(cis) %in% dist_pars, 1] <-
      purrr::map2_dbl(cis[rownames(cis) %in% dist_pars, 1], x$dlist$inv.transforms, ~.y(.x))

    cis[rownames(cis) %in% dist_pars, 2] <-
      purrr::map2_dbl(cis[rownames(cis) %in% dist_pars, 2], x$dlist$inv.transforms, ~.y(.x))
  }

  if ("coefs.exp" %in% transform & conf.int) {
    coefs[terms %in% coef_pars] <- exp(coefs[terms %in% coef_pars])
    ses[terms %in% coef_pars] <- exp(ses[terms %in% coef_pars])
    cis[terms %in% coef_pars] <- exp(cis[terms %in% coef_pars])
  } else if ("coefs.exp" %in% transform) {
    coefs[terms %in% coef_pars] <- exp(coefs[terms %in% coef_pars])
    ses[terms %in% coef_pars] <- exp(ses[terms %in% coef_pars])
  }

  lcl <- cis[, 1]
  ucl <- cis[, 2]

  res <- tibble::tibble(
    term = terms,
    estimate = coefs,
    std.error = ses,
    statistic = stats,
    p.value = pvals,
    conf.low = !!lcl,
    conf.high = !!ucl
  )

  if (pars == "coefs") {
    res[res$term %in% coef_pars, ]
  } else if (pars == "bdist") {
    res[res$term %in% dist_pars, ]
  } else {
    res
  }
}

#' @importFrom generics tidy
#' @export
generics::tidy

#' Glance at a flexsurv model object
#'
#' Glance accepts a model object and returns a tibble with exactly one row of model summaries.
#'
#' @param x Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted survival model object.
#'
#' @param ... Not currently used.
#'
#' @return A one-row \code{\link{tibble}} containing columns:
#'
#' * \code{N} Number of observations used in fitting
#'
#' * \code{events} Number of events
#'
#' * \code{censored} Number of censored events
#'
#' * \code{trisk} Total length of time-at-risk (i.e. follow-up)
#'
#' * \code{df} Degrees of freedom (i.e. number of estimated parameters)
#'
#' * \code{logLik} Log-likelihood
#'
#' * \code{AIC} Akaike's "An Information Criteria"
#'
#' * \code{BIC} Bayesian Information Criteria
#'
#' @export
#'
#' @md
#'
#' @examples
#' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
#' glance(fitg)
#'
glance.flexsurvreg <- function(x, ...)
{
  tibble::tibble(
    N = x$N,
    events = x$events,
    censored = x$N - x$events,
    trisk = x$trisk,
    df = x$npars,
    # global test statistic?
    # global p.value?
    logLik = x$loglik,
    AIC = x$AIC,
    BIC = BIC(x)
  )
}

#' @importFrom generics glance
#' @export
generics::glance

#' Augment data with information from a flexsurv model object
#'
#' Augment accepts a model object and a dataset and adds information about each observation in the dataset. Most commonly, this includes predicted values in the \code{.fitted} column, residuals in the \code{.resid} column, and standard errors for the fitted values in a \code{.se.fit} column. New columns always begin with a . prefix to avoid overwriting columns in the original dataset.
#'
#' @param x Output from \code{\link{flexsurvreg}} or \code{\link{flexsurvspline}}, representing a fitted survival model object.
#'
#' @param data A \code{\link{data.frame}} or \code{\link{tibble}} containing the original data that was used to produce the object \code{x}.
#'
#' @param newdata A \code{\link{data.frame}} or \code{\link{tibble}} containing all the original predictors used to create \code{x}. Defaults to \code{NULL}, indicating that nothing has been passed to \code{newdata}. If \code{newdata} is specified, the \code{data} argument will be ignored.
#'
#' @param type.predict Character indicating type of prediction to use. Passed to the \code{type} argument of the \code{\link{predict}} generic. Allowed arguments vary with model class, so be sure to read the \code{predict.my_class} documentation.
#'
#' @param type.residuals Character indicating type of residuals to use. Passed to the type argument of \code{\link{residuals}} generic. Allowed arguments vary with model class, so be sure to read the \code{residuals.my_class} documentation.
#'
#' @param ... Additional arguments. Not currently used.
#'
#' @details If neither of \code{data} or \code{newdata} are specified, then \code{model.frame(x)} will be used. It is worth noting that \code{model.frame(x)} will include a \code{\link{Surv}} object and not the original time-to-event variables used when fitting the \code{flexsurvreg} object. If the original data is desired, specify \code{data}.
#'
#' @return A \code{\link{tibble}} containing \code{data} or \code{newdata} and possible additional columns:
#'
#' * \code{.fitted} Fitted values of model
#'
#' * \code{.se.fit} Standard errors of fitted values
#'
#' * \code{.resid} Residuals (not present if \code{newdata} specified)
#'
#' @importFrom tidyr unnest
#'
#' @md
#'
#' @export
#'
#' @examples
#' fit <- flexsurvreg(formula = Surv(time, status) ~ age, data = lung, dist = "exp")
#' augment(fit, data = lung)
#'
augment.flexsurvreg <- function(x, data = NULL, newdata = NULL,
                                type.predict = "response",
                                type.residuals = "response", ...)
{
  if (is.null(data) && is.null(newdata)) {
    data <- model.frame(x)
  }

  type.predict <- match.arg(type.predict, c("response"))
  type.residuals <- match.arg(type.residuals, c("response"))

  if (is.null(newdata)) {
    predictions <- tidyr::unnest(predict(x, type = type.predict, se.fit = TRUE),
                                 .pred)
    res <- tibble::as_tibble(data)
    res <- tibble::add_column(res,
                              .fitted = predictions$.pred,
                              .se.fit = predictions$.std_error,
                              .resid = residuals(x, type = type.residuals))
  } else {
    predictions <- tidyr::unnest(predict(x, type = type.predict, se.fit = TRUE),
                                 .pred)
    res <- tibble::as_tibble(newdata)
    res <- tibble::add_column(res,
                              .fitted = predictions$.pred,
                              .se.fit = predictions$.std_error,)
  }
  res
}

utils::globalVariables(".pred")

#' @importFrom generics augment
#' @export
generics::augment
