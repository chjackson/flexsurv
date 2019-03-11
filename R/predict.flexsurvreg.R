##' Prediction from flexible survival models
##'
##' Predict outcomes from flexible survival models at the covariate values given in \code{"newdata"}.
##'
##' @inheritParams summary.flexsurvreg
##'
##' @param type \code{"response"} for mean survival
##'
##' \code{"quantile"} for quantiles of the survival distribution specified by \code{p}
##'
##' \code{"link"} for fitted values of the location parameter, analogous to the linear predictor in generalized linear models.   \code{type="lp"} and \code{type="linear"} are synonyms.
##'
##' @param interval Return a confidence interval \code{interval="confidence"}.
##' 
##' @param se.fit Return standard errors
##'
##' @param p vector of quantiles to return when \code{type="quantile"}
##'
##' @param tidy If \code{TRUE} then a tidyverse-friendly format is returned.  This is a tibble with one row per individual in \code{newdata}.  If multiple predictions are requested per individual (e.g. multiple quantiles) then list columns are returned... TODO finalize.    If \code{FALSE} then a data frame with one row per individual and prediction is returned. 
##'
##' Otherwise a data frame is returned... 
##'
##' TODO remove covariates?
##'
##' If \code{FALSE} then a simple data frame is returned 
##' 
##' @seealso \code{\link{summary.flexsurvreg}}
##'
##' @return TODO document format
##'
##' @importFrom tidyr unnest
##' @importFrom tibble tibble
##'
##' @examples ## TODO document extraction with tidyverse tools
##'
##' \dontrun{
##' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist="gengamma")
##' predict(fitg)
##' predict(fitg, type="quantile", p=0.5)
##' library(tidyverse) # TODO sort out deps
##' predict(fitg, type="quantile", p=c(0.1, 0.9))
##' }
predict.flexsurvreg <- function(object,
                                newdata,
                                type="response",
                                interval="none",
                                se.fit=FALSE,
                                p = c(0.1, 0.9),
                                tidy=TRUE
                                )
{
    if (missing(newdata))
        newdata <- model.frame(object)
    type <- match.arg(type, c("response", "quantile", "link", "lp", "linear"))
    interval <- match.arg(interval, c("none", "confidence"))
    stype <- type
    if (type=="response") stype <- "mean"
    if (type %in% c("lp", "linear")) stype <- "link"
    multi_output <- (type=="quantile") && (length(p) > 1)
    if (stype %in% c("mean", "quantile", "link")){ # may not need this line 
        res <- summary.flexsurvreg(object=object, newdata=newdata, type=stype,
                                   quantiles=p, ci=(interval=="confidence"), se=se.fit,
                                   tidy=!multi_output)
    }
    else res <- NULL
    rename <- function(x){
        ## if depend on tidyverse stuff could use that for this
        oldnames <- c("est", "se",       "lcl",       "ucl")
        newnames  <- c("pred","std_error","pred_lower","pred_upper")
        for (i in seq_along(x))
            names(x)[names(x)==oldnames[i]] <- newnames[i]
        x
    }
    if (multi_output){ ## TODO nicer way in purrr I expect, learning exercise for me...
        res <- lapply(res, rename)
        res2 <- vector(ncol(res[[1]]), mode = "list")
        for (i in seq_along(res2))
            res2[[i]] <- lapply(res, function(x)x[,i])
        res2 <- do.call(tibble, res2)
        names(res2) <- names(res[[1]])
        res <- if (tidy) res2 else tidyr::unnest(res2)
    } else {
        res <- rename(res)
    }
    res
}
