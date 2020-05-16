##' Prediction from flexible survival models
##'
##' Predict outcomes from flexible survival models at the covariate values given in \code{"newdata"}.
##'
##' @inheritParams summary.flexsurvreg
##'
##' @param newdata Data frame containing covariate values to produce fitted
##' values for, or a list that can be coerced to such a data frame.  There
##' must be a column for every covariate in the model formula, and one row for
##' every combination of covariates the fitted values are wanted for.  If this is
##' omitted, then the original data used to fit the model are used, extracted by
##' \code{model.frame(object)}. 
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
##' @param tidy If \code{TRUE} then a tidyverse-friendly format is returned.  This is a tibble with one row per individual in \code{newdata}.  If multiple predictions are requested per individual (e.g. multiple quantiles) then list columns are returned.    If \code{FALSE} then a data frame with one row per individual and prediction is returned. 
##'
##' @param return_covs If \code{FALSE} (the default) then covariates are removed from the returned data, to match the behaviour of standard predict functions such as \code{\link{predict.lm}}.  If \code{TRUE}, then columns for the covariates are included in the returned data, which may be helpful for post-processing the results.
##'
##' @seealso \code{\link{summary.flexsurvreg}}
##'
##' @return If \code{tidy=TRUE} then a tibble with one row per individual in \code{newdatata}.  If multiple quantities are predicted per individual, then list columns are returned, see, e.g. \url{https://r4ds.had.co.nz/many-models.html#list-columns-1}.   If \code{tidy=FALSE}, then a data frame with one row per individual and prediction is returned. 
##'
##' @importFrom tidyr unnest
##' 
##' @importFrom tibble tibble
##'
##' @examples ## TODO as well as just examples of running the functions, document extraction with tidyverse tools for the examples that return list columns. 
##'
##' \dontrun{
##' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist="gengamma")
##'
##' ## Simplest prediction: mean or median, for covariates defined by original dataset
##' predict(fitg)
##' predict(fitg, type="quantile", p=0.5)
##' predict(fitg, return_covs=TRUE)
##'
##' ## Simple prediction for user-defined covariate values
##' predict(fitg, newdata=list(age=c(40, 50, 60)))
##' predict(fitg, type="quantile", p=0.5, newdata=list(age=c(40,50,60)))
##'
##' ## Predict multiple quantiles 
##' require(tibble) 
##' require(tidyr) # TODO declare dependencies in DESCRIPTION, decide if Imports or Suggests
##' predict(fitg, type="quantile", p=c(0.1, 0.9))
##'
##' ## Multiple quantiles and multiple covariate values
##' 
##' }
predict.flexsurvreg <- function(object,
                                newdata,
                                type="response",
                                interval="none",
                                se.fit=FALSE,
                                p = c(0.1, 0.9),
                                tidy=TRUE,
                                return_covs=FALSE
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
        ## If depend on tidyverse stuff could use dplyr::rename. Should only add dependencies where necessary, though if we are making list columns, we could assume the user is willing to load tidyverse to deal with them.
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
        res2 <- do.call(tibble, res2) # TODO This bit's currently broken
        names(res2) <- names(res[[1]])
        res <- if (tidy) res2 else tidyr::unnest(res2)
    } else {
        res <- rename(res)
    }
    if (!return_covs){
        covnames <- attr(model.frame(object),"covnames.orig")
        res[,covnames] <- NULL
    }
    res
}
