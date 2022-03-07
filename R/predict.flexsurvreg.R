#' Predictions from flexible survival models
#'
#' Predict outcomes from flexible survival models at the covariate values
#'   specified in \code{newdata}.
#'
#' @param object Output from \code{\link{flexsurvreg}} or
#'   \code{\link{flexsurvspline}}, representing a fitted survival model object.
#'
#' @param newdata Data frame containing covariate values at which to produce
#'   fitted values. There must be a column for every covariate in the model
#'   formula used to fit \code{object}, and one row for every combination of
#'   covariate values at which to obtain the fitted predictions.
#'
#'   If \code{newdata} is omitted, then the original data used to fit the model
#'   are used, as extracted by \code{model.frame(object)}. However this will
#'   currently not work if the model formula contains functions, e.g.
#'   \code{~ factor(x)}. The names of the model frame must correspond to
#'   variables in the original data.
#'
#' @param type Character vector for the type of predictions desired.
#'
#' * \code{"response"} for mean survival time (the default). \code{"mean"} is
#'   an acceptable synonym
#'
#' * \code{"quantile"} for quantiles of the survival distribution as specified
#'   by \code{p}
#'
#' * \code{"rmst"} for restricted mean survival time
#'
#' * \code{"survival"} for survival probabilities
#'
#' * \code{"cumhaz"} for cumulative hazards
#'
#' * \code{"hazard"} for hazards
#'
#' * \code{"link"} for fitted values of the location parameter, analogous to
#'   the linear predictor in generalized linear models (\code{type = "lp"} and
#'   \code{type = "linear"} are acceptable synonyms)
#'
#' @param times Vector of time horizons at which to compute fitted values.
#'   Only applies when \code{type} is \code{"survival"}, \code{"cumhaz"},
#'   \code{"hazard"}, or \code{"rmst"}. Will be silently ignored for all other
#'   types.
#'
#'   If not specified, predictions for \code{"survival"}, \code{"cumhaz"}, and
#'   \code{"hazard"} will be made at each observed event time in
#'   \code{model.frame(object)}.
#'
#'   For \code{"rmst"}, when \code{times} is not specified predictions will be
#'   made at the maximum observed event time from the data used to fit
#'   \code{object}. Specifying \code{times = Inf} is valid, and will return
#'   mean survival (equal to \code{type = "response"}).
#'
#' @param start Optional left-truncation time or times. The returned
#'   survival, hazard, or cumulative hazard will be conditioned on survival up
#'   to this time. `start` must be length 1 or the same length as `times`.
#'
#' @param conf.int Logical. Should confidence intervals be returned?
#'   Default is \code{FALSE}.
#'
#' @param conf.level Width of symmetric confidence intervals, relative to 1.
#'
#' @param se.fit Logical. Should standard errors of fitted values be returned?
#'   Default is \code{FALSE}.
#'
#' @param p Vector of quantiles at which to return fitted values when
#'   \code{type = "quantile"}. Default is \code{c(0.1, 0.9)}.
#'
#' @param ... Not currently used.
#'
#' @return A \code{\link{tibble}} with same number of rows as \code{newdata}
#'   and in the same order. If multiple predictions are requested, a
#'   \code{\link{tibble}} containing a single list-column of data frames.
#'
#'   For the list-column of data frames - the dimensions of each data frame
#'   will be identical. Rows are added for each value of \code{times} or
#'   \code{p} requested.
#'
#' @seealso \code{\link{summary.flexsurvreg}},
#'   \code{\link{residuals.flexsurvreg}}
#'
#' @importFrom tibble tibble
#' @importFrom stats predict
#'
#' @md
#'
#' @export
#'
#' @examples
#'
#' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist = "gengamma")
#'
#' ## Simplest prediction: mean or median, for covariates defined by original dataset
#' predict(fitg)
#' predict(fitg, type = "quantile", p = 0.5)
#'
#' ## Simple prediction for user-defined covariate values
#' predict(fitg, newdata = data.frame(age = c(40, 50, 60)))
#' predict(fitg, type = "quantile", p = 0.5, newdata = data.frame(age = c(40,50,60)))
#'
#' ## Predict multiple quantiles and unnest
#' require(tidyr)
#' pr <- predict(fitg, type = "survival", times = c(600, 800))
#' tidyr::unnest(pr, .pred)
#'
predict.flexsurvreg <- function(object,
                                newdata,
                                type = "response",
                                times,
                                start = 0,
                                conf.int = FALSE,
                                conf.level = 0.95,
                                se.fit = FALSE,
                                p = c(0.1, 0.9),
                                ...
                                )
{
    if (missing(newdata)) newdata <- model.frame(object)

    assertthat::assert_that(inherits(newdata, "data.frame"),
                            msg = "`newdata` must inherit class `data.frame`")
    assertthat::assert_that(is.logical(conf.int), is.logical(se.fit))
    assertthat::assert_that(all(is.numeric(p), p <= 1, p >=0),
                            msg = "`p` should be a vector of quantiles between 0 and 1")

    assertthat::assert_that(
        is.numeric(start),
        msg = "`start` must be a numeric vector of left-truncation times"
    )

    if (conf.int) assertthat::assert_that(is.numeric(conf.level),
                                          conf.level > 0, conf.level < 1,
                                          length(conf.level) == 1,
                                          msg = "`conf.level` must be length one and between 0 and 1")

    type <- match.arg(type, c("mean", "response", "quantile", "link", "lp", "linear",
                              "survival", "cumhaz", "hazard", "rmst"))

    stype <- switch(
        type,
        response = "mean",
        lp = "link",
        linear = "link",
        type # all others keep their type
    )

    if (stype %in% c("survival", "cumhaz", "hazard")) {
        if (missing(times)) times <- object$data$Y[, 1][order(object$data$Y[, 1])]
        assertthat::assert_that(all(is.numeric(times), times > 0),
                                msg = "`times` must be a vector of positive real-valued numbers.")
    } else if (stype == "rmst" && !missing(times)) {
        assertthat::assert_that(all(is.numeric(times), times > 0),
                                msg = "`times` must be a vector of positive real-valued numbers.")
    } else {
        times <- NULL
    }

    assertthat::assert_that(
        length(start) == 1 | length(start) == length(times),
        msg =
            paste0(
                "Length of `start` is ", length(start), ". Length should be 1, or the same length as `times`, which is ", length(times)
            )
    )

    nest_output <- ((stype == "quantile" && length(p) > 1) |
                        (stype %in% c("survival", "cumhaz", "hazard", "rmst") &&
                        length(times) > 1))

    res <- if (stype %in% c("survival", "hazard", "cumhaz", "rmst")) {
      summary(
        object = object,
        newdata = newdata,
        type = stype,
        quantiles = p,
        t = times,
        start = start,
        ci = conf.int,
        cl = conf.level,
        se = se.fit,
        tidy = TRUE,
        na.action = na.pass,
        cross = TRUE
      )
    } else {
      # Avoid passing `t = times` for non-time based predictions
      # to avoid noisy warnings
      summary(
        object = object,
        newdata = newdata,
        type = stype,
        quantiles = p,
        start = start,
        ci = conf.int,
        cl = conf.level,
        se = se.fit,
        tidy = TRUE,
        na.action = na.pass,
        cross = TRUE
      )
    }

    res <- tidy_rename(res, stype)

    if (nest_output) {
      if (stype == 'quantile') {
        num_reps <- length(p)
      } else {
        num_reps <- length(times)
      }
      orig_nrow <- nrow(newdata)
      res <- dplyr::mutate(res, .id = rep(1:orig_nrow, each = num_reps))
      res <- dplyr::group_nest(res, .id, .key = '.pred')
      res <- dplyr::select(res, .pred)
    }
    res
}

utils::globalVariables('.id')

tidy_names <- function() {
    tibble::tibble(
        old = c("time", "quantile", "est", "se", "lcl", "ucl"),
        new = c(".time", ".quantile", ".pred", ".std_error",
                ".pred_lower", ".pred_upper")
    )
}

tidy_rename <- function(x, type) {
  names_tbl <- tidy_names()
  names_to_change <- names_tbl[names_tbl$old %in% names(x), ]
  out <- dplyr::select(x, names_to_change$old)
  colnames(out) <- names_to_change$new
  out <- tibble::as_tibble(out)
  if (type == "mean") type <- "time"
  name_with_type <- rlang::sym(paste0('.pred_', type))
  dplyr::rename(out, !!name_with_type := .pred)
}
