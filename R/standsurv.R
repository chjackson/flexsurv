#' Marginal survival and hazards of fitted flexsurvreg models
#'
#' Returns a tidy data.frame of marginal survival probabilities, or hazards, 
#' restricted mean survival, or quantiles of the marginal survival function
#' at user-defined time points and covariate patterns.
#' Standardization is performed over any undefined covariates in the model. 
#' The user provides the data to standardize over. Contrasts can be calculated 
#' resulting in estimates of the average treatment effect or the average 
#' treatment effect in the treated if a treated subset of the data are supplied.
#'
#' The syntax of \code{standsurv} follows closely that of Stata's 
#' \code{standsurv} command written by Paul Lambert and Michael Crowther. The 
#' function calculates standardized (marginal) measures including standardized
#' survival functions, standardized restricted mean survival times, quantiles
#' and the hazard of standardized survival. The standardized survival is defined as
#' \deqn{S_s(t|X=x) = E(S(t|X=x,Z)) = \frac{1}{N} \sum_{i=1}^N S(t|X=x,Z=z_i)}{S(t|X=x) = E[S(t|X=x,Z)] = 1/N * sum(S(t|X=x,Z=z_i))}
#' The hazard of the standardized survival is a weighted average of 
#' individual hazard functions at time t, weighted by the survival
#' function at this time:
#' \deqn{h_s(t|X=x) = \frac{\sum_{i=1}^N S(t|X=x,Z=z_i)h(t|X=x,Z=z_i)}{\sum_{i=1}^N S(t|X=x,Z=z_i)}}{h(t|X=x) = sum(S(t|X=x,Z=z_i) * h(t|X=x,Z=z_i)) / sum(S(t|X=x,Z=z_i))}
#' Marginal expected survival and hazards can be calculated by providing a 
#' population-based lifetable of class ratetable in \code{ratetable} and a 
#' mapping between stratification factors in the lifetable and the user dataset
#' using \code{rmap}. If these stratification factors are not in the fitted
#' survival model then the user must specify them in \code{newdata} along with
#' the covariates of the model. The marginal expected survival is calculated 
#' using the "Ederer" method that assumes no censoring as this is most relevant 
#' approach for forecasting (see 
#' \code{\link[survival]{survexp}}). A worked example is given below.
#' 
#' Marginal all-cause survival and hazards can be calculated after fitting a
#' relative survival model, which utilise the expected survival from a population
#' ratetable. See Rutherford et al. (Chapter 6) for further details.
#' 
#'
#' @param object Output from \code{\link{flexsurvreg}} or
#' \code{\link{flexsurvspline}}, representing a fitted survival model object.
#' @param newdata Data frame containing covariate values to produce marginal
##' values for. If not specified then the fitted model data.frame is used.
##' There must be a column for every covariate in the model formula
##' for which the user wishes to standardize over.  These are in the same format
##' as the original data, with factors as a single variable, not 0/1 contrasts.
##' Any covariates that are to be fixed should be specified in \code{at}. 
##' There should be one row for every combination of covariates in which to 
##' standardize over. If newdata contains a variable named '(weights)' then a 
##' weighted mean will be used to create the standardized estimates. This is the
##' default behaviour if the fitted model contains case weights, which are stored 
##' in the fitted model data.frame.
#' @param at A list of scenarios in which specific covariates are fixed to 
#' certain values. Each element of \code{at} must itself be a list. For example,
#' for a covariate \code{group} with levels "Good", "Medium" and "Poor", the 
#' standardized survival plots for each group averaging over all other 
#' covariates is specified using 
#' \code{at=list(list(group="Good"), list(group="Medium"), list(group="Poor"))}.
#' @param atreference The reference scenario for making contrasts. Default is 1
#' (i.e. the first element of \code{at}).
#' @param type \code{"survival"} for marginal survival probabilities. In a 
#' relative survival framework this returns the marginal all-cause survival 
#' (see details).
##' 
##' \code{"hazard"} for the hazard of the marginal survival probability. In a 
#' relative survival framework this returns the marginal all-cause hazard 
#' (see details).
##' 
##' \code{"rmst"} for standardized restricted mean survival.
##' 
##' \code{"relsurvival"} for marginal relative survival (can only be specified
##' if a relative survival model has been fitted in flexsurv).
##' 
##' \code{"excesshazard"} for marginal excess hazards (can only be specified
##' if a relative survival model has been fitted in flexsurv).
##' 
##' \code{"quantile"} for quantiles of the marginal all-cause survival 
##' distribution. The \code{quantiles} option also needs to be provided.
#' @param t Times to calculate marginal values at.
#' @param ci Should confidence intervals be calculated? 
#' Defaults to FALSE
#' @param se Should standard errors be calculated? 
#' Defaults to FALSE
#' @param boot Should bootstrapping be used to calculate standard error and 
#' confidence intervals? Defaults to FALSE, in which case the delta method is 
#' used
#' @param B Number of bootstrap simulations from the normal asymptotic 
#' distribution of the estimates used to calculate confidence intervals or 
#' standard errors. Decrease for greater speed at the expense of accuracy. Only 
#' specify if \code{boot = TRUE}
#' @param cl Width of symmetric confidence intervals, relative to 1.
#' @param trans Transformation to apply when calculating standard errors via the
#' delta method to obtain confidence intervals. The default transformation is 
#' "log". Other possible names are "none", "loglog", "logit".
#' @param contrast Contrasts between standardized measures defined by \code{at}
#' scenarios. Options are \code{"difference"} and \code{"ratio"}. There will be
#' n-1 new columns created where n is the number of \code{at} scenarios. Default
#' is NULL (i.e. no contrasts are calculated).
#' @param trans.contrast Transformation to apply when calculating standard errors
#' for contrasts via the delta method to obtain confidence intervals. The default
#' transformation is "none" for differences in survival, hazard, quantiles, or RMST, 
#' and "log" for ratios of survival, hazard, quantiles or RMST.
#' @param seed The random seed to use (for bootstrapping confidence intervals)
#' @param rmap An list that maps data set names to expected ratetable names. 
#' This must be specified if all-cause survival and hazards are required after
#' fitting a relative survival model. This can also be specified if expected
#' rates are required for plotting purposes. See the details section below.
#' @param ratetable A table of expected event rates 
#' (see \code{\link[survival]{ratetable}})
#' @param scale.ratetable Transformation from the time scale of the fitted 
#' flexsurv model to the time scale in \code{ratetable}. For example, if the 
#' analysis time of the fitted model is in years and the ratetable is in 
#' units/day then we should use \code{scale.ratetable = 365.25}. This is the 
#' default as often the ratetable will be in units/day (see example).
#' @param n.gauss.quad Number of Gaussian quadrature points used for integrating 
#' the all-cause survival function when calculating RMST in a relative survival 
#' framework (default = 100)
#' @param quantiles If \code{type="quantile"}, this specifies the quantiles of 
#' the survival time distribution to return estimates for.
#' @param interval Interval of survival times for quantile root finding. 
#' Default is c(1e-08, 500). 
#'
#' @return A \code{tibble} containing one row for each 
#' time-point. The column naming convention is \code{at{i}} for the ith scenario
#' with corresponding confidence intervals (if specified) named \code{at{i}_lci}
#' and \code{at{i}_uci}. Contrasts are named \code{contrast{k}_{j}} for the 
#' comparison of the kth versus the jth \code{at} scenario.
#' 
#' In addition tidy long-format data.frames are returned in the attributes
#' \code{standsurv_at} and \code{standsurv_contrast}. These can be passed to 
#' \code{ggplot} for plotting purposes (see \code{\link{plot.standsurv}}).
#' @importFrom tibble as_tibble
#' @importFrom rlang :=
#' @importFrom dplyr bind_cols
#' @importFrom dplyr inner_join
#' @export
#' @author Michael Sweeting <mikesweeting79@@gmail.com>
#' @references Paul Lambert, 2021. "STANDSURV: Stata module to compute 
#' standardized (marginal) survival and related functions," 
#' Statistical Software Components S458991, Boston College Department of 
#' Economics. https://ideas.repec.org/c/boc/bocode/s458991.html
#' 
#' Rutherford, MJ, Lambert PC, Sweeting MJ, Pennington B, Crowther MJ, Abrams KR,
#' Latimer NR. 2020. "NICE DSU Technical Support Document 21: Flexible Methods 
#' for Survival Analysis" 
#' https://nicedsu.sites.sheffield.ac.uk/tsds/flexible-methods-for-survival-analysis-tsd
#' 
#' @examples
#'## mean age is higher in those with smaller observed survival times 
#' newbc <- bc
#' set.seed(1)
#' newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE),
#'  sd = 5)
#' 
#' ## Fit a Weibull flexsurv model with group and age as covariates
#' weib_age <- flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc, 
#'                        dist="weibull")
#'                        
#'## Calculate standardized survival and the difference in standardized survival
#'## for the three levels of group across a grid of survival times                        
#'standsurv_weib_age <- standsurv(weib_age, 
#'                                            at = list(list(group="Good"), 
#'                                                      list(group="Medium"), 
#'                                                      list(group="Poor")), 
#'                                            t=seq(0,7, length.out=100),
#'                                            contrast = "difference", ci=FALSE)
#'standsurv_weib_age
#'
#'## Calculate hazard of standardized survival and the marginal hazard ratio
#'## for the three levels of group across a grid of survival times
#'## 10 bootstraps for confidence intervals (this should be larger)
#'\dontrun{          
#'haz_standsurv_weib_age <- standsurv(weib_age, 
#'                                            at = list(list(group="Good"), 
#'                                                      list(group="Medium"), 
#'                                                      list(group="Poor")), 
#'                                            t=seq(0,7, length.out=100),
#'                                            type="hazard",
#'                                            contrast = "ratio", boot = TRUE,
#'                                            B=10, ci=TRUE)
#'haz_standsurv_weib_age                                            
#'plot(haz_standsurv_weib_age, ci=TRUE)
#'## Hazard ratio plot shows a decreasing marginal HR 
#'## Whereas the conditional HR is constant (model is a PH model)
#'plot(haz_standsurv_weib_age, contrast=TRUE, ci=TRUE)
#'
#'## Calculate standardized survival from a Weibull model together with expected
#'## survival matching to US lifetables
#'
#'# age at diagnosis in days. This is required to match to US ratetable, whose
#'# timescale is measured in days
#'newbc$agedays <- floor(newbc$age * 365.25)  
#'## Create some random diagnosis dates centred on 01/01/2010 with SD=1 year
#'## These will be used to match to expected rates in the lifetable
#'newbc$diag <- as.Date(floor(rnorm(dim(newbc)[1], 
#'                      mean = as.Date("01/01/2010", "%d/%m/%Y"), sd=365)), 
#'                      origin="1970-01-01")
#'## Create sex (assume all are female)
#'newbc$sex <- factor("female")
#'standsurv_weib_expected <- standsurv(weib_age, 
#'                                            at = list(list(group="Good"), 
#'                                                      list(group="Medium"), 
#'                                                      list(group="Poor")), 
#'                                            t=seq(0,7, length.out=100),
#'                                            rmap=list(sex = sex,
#'                                                      year = diag,
#'                                                      age = agedays),
#'                                            ratetable = survival::survexp.us,
#'                                            scale.ratetable = 365.25,
#'                                            newdata = newbc)
#'## Plot marginal survival with expected survival superimposed                                            
#'plot(standsurv_weib_expected, expected=TRUE)
#'}
standsurv <- function(object, newdata = NULL, at = list(list()), atreference = 1, 
                                  type = "survival", t = NULL, 
                                  ci = FALSE, se = FALSE, boot = FALSE, B = NULL, 
                                  cl =0.95, trans = "log", contrast = NULL, 
                                  trans.contrast = NULL, seed = NULL, 
                                  rmap, ratetable, scale.ratetable = 365.25, 
                                  n.gauss.quad = 100, quantiles = 0.5,
                                  interval = c(1e-08, 500)) {
  x <- object
  
  if(!is.null(seed)) set.seed(seed)
  
  ## Add checks
  ## Currently type is restricted to survival, hazard, rmst, quantile, or relsurvival or excesshazard (for relative survival models)
  type <- match.arg(type, c("survival", "hazard", "rmst", "quantile", "relsurvival", "excesshazard"))
  # Check that models is a relative survival model if relsurv or excesshazard are specified
  if(type %in% c("relsurvival", "excesshazard") & !("bhazard" %in% names(x$call))){
    stop(paste0(type, " can only be specified for a relative survival flexsurv model"))
  }
  type2 <- type
  
  ## Checks for relative survival models
  if("bhazard" %in% names(x$call)){
    if(type %in% c("quantile", "survival", "hazard") & (missing(rmap) | missing(ratetable))) 
      stop("'rmap' and 'ratetable' must be specified to calculate all-cause survival/hazard in a relative survival model")
    # Swap type to acsurvival (all-cause survival) or achazard (all-cause hazard) if survival or hazard is specified
    # Swap type to survival or hazard if relsurv or excesshazard are specified
    type2 <- switch(type, 
           "survival"= "acsurvival",
           "hazard"= "achazard",
           "quantile" = "acquantile",
           "relsurvival" = "survival",
           "excesshazard" = "hazard",
           "rmst" = "acrmst"
    )
    
    if(type2 == "acsurvival") message("Marginal all-cause survival will be calculated")
    if(type2 == "achazard") message("Marginal all-cause hazard will be calculated")
    if(type2 == "acquantile") message("Quantiles of marginal all-cause survival will be calculated")
    if(type2 == "acrmst"){
      message("Marginal restricted mean survival will be calculated")
      if(length(t)>2){
        message("Marginal RMST is currently slow. We suggest using only one or two time points")
      }
    }
    
    if(!missing(rmap) & is.null(newdata)) 
      stop("Must provide a 'newdata' data.frame containing all covariates and matching variables with 'rmap'")
    
  }
  
  # checks for type="quantile"
  if(type == "quantile" & !is.numeric(quantiles)) 
    stop("'quantiles' argument must be provided as a numeric vector for type='quantile'")
  
  if(!is.null(contrast)) {
    contrast <- match.arg(contrast, c("difference", "ratio"))
    if(is.null(trans.contrast)){
      trans.contrast <- ifelse(contrast=="difference", "none", "log")
      trans.contrast <- match.arg(trans.contrast, c("log", "none", "loglog", "logit"))
    }
  }
  trans <- match.arg(trans, c("log", "none", "loglog", "logit"))
  
  ## Check that at is a list and that all elements of at are lists
  if(!is.list(at)){
    stop("'at' must be at list")
  }
  if (!all(sapply(at, is.list))) {
    stop("All elements of 'at' must be lists")
  }
  
  ## Check sensible transformations have been specified for type
  if(boot == F & type %in% c("hazard", "rmst", "quantile") & trans %in% c("loglog", "logit")){
    warning(paste0("type ",type, " with transformation ",trans, " may not be sensible"))
  }
  if(boot == F & !is.null(B)){
    stop(paste0("'boot'=FALSE but 'B' is non-null"))
  }
  if(boot == T & is.null(B)){
    stop(paste("'B' must be specified if 'boot'=TRUE"))
  }
  
  ## Contrast numbers
  cnums <- (1:length(at))[-atreference]

  ## If no contrasts (i.e. only one at() list) then 'contrast' should be NULL
  if(length(cnums)==0 & !is.null(contrast)){
    stop("'contrast' cannot be specified if length of 'at' < 2")
  }
  
  ## Standardize over fitted dataset by default
  if(is.null(newdata)){
    data <- model.frame(x)
  } else{
    data <- newdata
  }
  ## Was weighted regression used, and do these weights feature in newdata?
  weighted <- FALSE
  if("(weights)" %in% names(data)){
    if(var(data$`(weights)`)!=0){
      weighted <- TRUE
      message("Weighted regression was used, standardization will be weighted accordingly")
    }
  }

  ## Set t to be unique event times if NULL  
  if(is.null(t)) t <- sort(unique(x$data$Y[,"stop"]))
  
  ## Calculate individual expected survival and expected hazard 
  if((!missing(rmap) & !missing(ratetable))){
    message("Calculating marginal expected survival and hazard")
    expsurv <- expsurv.fn(t, substitute(rmap), ratetable, data, weighted, scale.ratetable)
  } else expsurv <- NULL
  
  ## Loop over at()
  stand.pred.list <- dat.list <- list()
  for(i in 1:length(at)){
    dat <- data
    covs <- at[[i]]
    covnames <- names(covs)
    ## If all covariates have been specified in 'at' then we have no further covariates
    ## to standardize over, so just use 1 row of data
    ## do not use this shortcut for type2 = "acsurvival", type2 = "achazard", type = "acrmst"
    ## or type = "acquantile" as we require individual expected survivals for these methods
    allcovs <- all.vars(formula(x)[-2])
    if(all(allcovs %in% covnames) & !weighted & !(type2 %in% c("acsurvival", "achazard", "acrmst", "acquantile"))){
      dat <- dat[1,,drop=F]
    }
    ## If at is not specified then no further manipulation of data is required, 
    ## we standardize over original or passes dataset
    if(!is.null(covnames)){
      for(j in 1:length(covnames)) dat[, covnames[j]] <- covs[j]
    } 
    dat.list[[i]] <- dat
    predsum <- standsurv.fn(object, type = type2, newdata=dat, t=t, i=i, weighted=weighted, expsurv=expsurv,
                            rmap=substitute(rmap), ratetable=ratetable, scale.ratetable=scale.ratetable,
                            quantiles=quantiles, interval=interval, n.gauss.quad=n.gauss.quad)
    
    if(ci == TRUE | se == TRUE){
      
      if(boot == TRUE){
        if(i==1)
          message("Calculating bootstrap standard errors / confidence intervals")      
        rawsim <- NULL
        bootresults <- boot.standsurv(object, B, dat, i, t, type, type2, weighted, 
                                      se, ci, cl, rawsim, predsum, expsurv,
                                      rmap=substitute(rmap), ratetable=ratetable, scale.ratetable=scale.ratetable,
                                      quantiles=quantiles, interval=interval, n.gauss.quad=n.gauss.quad)  
        stand.pred.list[[i]] <- bootresults$stand.pred
        predsum <- bootresults$predsum
        if(is.null(rawsim)) 
          rawsim <- bootresults$rawsim
      } else{
        if(i==1) message("Calculating standard errors / confidence intervals using delta method")
        predsum <- deltamethod.standsurv(object, newdata=dat, type2, t, i, se, ci, predsum, 
                                         trans, cl, weighted, expsurv,
                                         rmap=substitute(rmap), ratetable=ratetable, 
                                         scale.ratetable=scale.ratetable,
                                         quantiles=quantiles, interval=interval, 
                                         n.gauss.quad=n.gauss.quad)
      }
    }
    
    if(i == 1) {
      standpred <- predsum
    } else {
      by <- ifelse(type2 %in% c("quantile","acquantile"), "probability", "time")
      standpred <- standpred %>% inner_join(predsum, by = by)
    }
    
  }
  if (anyNA(standpred)) warning("Missing values present in newdata")
  
  if(!is.null(contrast)){
    if(boot == TRUE){
      if(ci==TRUE | se==TRUE) message("Calculating bootstrap standard errors / confidence intervals for contrasts")
      if(contrast == "difference"){
        for(i in cnums){
          standpred <- standpred %>% mutate("contrast{i}_{atreference}" := .data[[paste0("at", i)]] - .data[[paste0("at", atreference)]])
          if(ci == TRUE){
            stand.pred.quant <- apply(stand.pred.list[[i]] - stand.pred.list[[atreference]], 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
            stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% rename("contrast{i}_{atreference}_lci" := "2.5%", "contrast{i}_{atreference}_uci" := "97.5%")
            standpred <- standpred %>% bind_cols(stand.pred.quant)
          }
          if(se == TRUE){
            stand.pred.se <- tibble("contrast{i}_{atreference}_se" := apply(stand.pred.list[[i]] - stand.pred.list[[atreference]], 2, sd, na.rm=TRUE))
            standpred <- standpred %>% bind_cols(stand.pred.se)
          }
        }
      }
      if(contrast == "ratio"){
        for(i in cnums){
          standpred <- standpred %>% mutate("contrast{i}_{atreference}" := .data[[paste0("at", i)]] / .data[[paste0("at", atreference)]])
          if(ci == TRUE){
            stand.pred.quant <- apply(stand.pred.list[[i]] / stand.pred.list[[atreference]], 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
            stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% rename("contrast{i}_{atreference}_lci" := "2.5%", "contrast{i}_{atreference}_uci" := "97.5%")
            standpred <- standpred %>% bind_cols(stand.pred.quant)
          }
          if(se == TRUE){
            stand.pred.se <- tibble("contrast{i}_{atreference}_se" := apply(stand.pred.list[[i]] / stand.pred.list[[atreference]], 2, sd, na.rm=TRUE))
            standpred <- standpred %>% bind_cols(stand.pred.se)
          }
        }
      }
    }
    else {
      if(ci==TRUE | se==TRUE) message("Calculating standard errors / confidence intervals for contrasts using delta method")
      for(i in cnums){
        standpred <- deltamethod.contrast.standsurv(object, dat=dat.list[[i]], dat.ref=dat.list[[atreference]],
                                       type2, t, i, atreference, se, ci, standpred, trans.contrast, 
                                       cl, contrast, weighted, expsurv,
                                       rmap=substitute(rmap), ratetable=ratetable, scale.ratetable=scale.ratetable,
                                       quantiles=quantiles, interval=interval, n.gauss.quad=n.gauss.quad)
      }
    }
  }
  label <- unlist(lapply(at,function(k){paste(names(k),k,sep="=",collapse=", ")}))
  attr(standpred, "label") <- label
  attr(standpred, "type") <- type
  attr(standpred, "contrast") <- contrast
  attr(standpred, "at") <- at
  attr(standpred, "atreference") <- atreference
  if(!is.null(expsurv)){
    attr(standpred, "expected") <- expsurv$marginal
  }
  class(standpred) <- c("standsurv", class(standpred))
  
  ## Create tidy versions of the data.frame and store as attributes
  standpred <- tidy(standpred)
  
  standpred
}

#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr left_join
#' @importFrom dplyr slice
#' @importFrom dplyr n
#' @importFrom dplyr ungroup
#' @importFrom statmod gauss.quad
#' @import rlang
standsurv.fn <- function(object, type, newdata, t, i, trans="none", weighted, expsurv, rmap, ratetable, 
                         scale.ratetable, quantiles, interval, n.gauss.quad){
  tr.fun <- tr(trans)
  if(! type %in% c("hazard", "acsurvival", "achazard", "acrmst", "quantile", "acquantile")){
    predsum <- standsurv.fn.generic(t, object, type, newdata, i, weighted, tr.fun)
  } else if(type == "hazard"){
    predsum <- standsurv.fn.hazard(t, object, newdata, i, weighted, tr.fun)
  } else if(type == "quantile"){
    predsum <- standsurv.fn.quantile(object, newdata, i, weighted, tr.fun, 
                                     quantiles, interval)
  } else if(type=="acsurvival"){
    predsum <- standsurv.fn.acsurvival(t, object, newdata, i, weighted, tr.fun, expsurv)
  } else if(type=="achazard"){
    predsum <- standsurv.fn.achazard(t, object, newdata, i, weighted, tr.fun, expsurv)
  } else if(type=="acrmst"){
    pred <- standsurv.fn.acrmst(t, object, newdata, i, rmap, ratetable, scale.ratetable, weighted, tr.fun, n.gauss.quad)
    predsum <- tibble(time=t, "at{i}":=pred)
  } else if(type=="acquantile"){
    predsum <- standsurv.fn.acquantile(t, object, newdata, i,  rmap, ratetable, scale.ratetable,
                                       weighted, tr.fun, quantiles, interval)
  }
  predsum
}

# generic standardisation function, for use with types "survival", "relsurv", "excesshazard"
standsurv.fn.generic <- function(t, object, type, newdata, i, weighted, tr.fun){
  pred <- summary(object, type = type, tidy = T, newdata=newdata, t=t, ci=F)
  pred$levels.fct <- factor(seq_along(t))
  if(weighted){
    pred$weights <- rep(newdata$`(weights)`, each=length(t))
    predsum <- pred %>% group_by(levels.fct, time) %>%
      summarise("at{i}" := tr.fun(weighted.mean(.data$est, .data$weights))) %>%
      ungroup() %>% select(-levels.fct)
  } else{
    predsum <- pred %>% group_by(levels.fct, time) %>% summarise("at{i}" := tr.fun(mean(.data$est))) %>%
      ungroup() %>% select(-levels.fct)
  }
  predsum
}

# hazard standardisation function, for use with type "hazard"
standsurv.fn.hazard <- function(t, object, newdata, i, weighted, tr.fun){
  pred <- summary(object, type = "hazard", tidy = T, newdata=newdata, t=t, ci=F)
  names(pred)[names(pred)=="est"] <- "h"
  pred <- cbind(pred, S = summary(object, type = "survival", tidy = T, newdata=newdata, t=t, ci=F)[,"est"])
  pred$levels.fct <- factor(seq_along(t))
  
  if(weighted){
    pred$weights <- rep(newdata$`(weights)`, each=length(t))
    predsum <- pred %>% group_by(levels.fct, time) %>% 
      summarise("at{i}" := tr.fun(weighted.mean(.data$h,.data$S * .data$weights))) %>%
      ungroup() %>% select(-levels.fct)
  } else{
    predsum <- pred %>% group_by(levels.fct, time) %>% summarise("at{i}" := tr.fun(weighted.mean(.data$h,.data$S))) %>% 
      ungroup() %>% select(-levels.fct)
  }
}

# quantile standardisation function, for use with type "quantile"
standsurv.fn.quantile <- function(object, newdata, i, weighted, tr.fun, quantiles, interval){
  quantile.root.fn <- function(t, q, object, newdata, i, weighted){
    as.numeric(standsurv.fn.generic(t, object, type = "survival", newdata, i, 
                                    weighted, tr.fun=tr("none"))[1,2]) - (1-q)
  }
  predsum <- tibble()
  for(q in quantiles){
    root <- uniroot(quantile.root.fn, interval = interval,
                    q=q, object=object, newdata=newdata, i=i, weighted=weighted)$root
    predsum <- predsum %>% bind_rows(tibble(probability = q) %>% mutate("at{i}" := tr.fun(root)))
  }
  predsum
}

# acsurvival standardisation function, for use with type "acsurvival"
standsurv.fn.acsurvival <- function(t, object, newdata, i, weighted, tr.fun, expsurv){
  rs <- summary(object, type = "survival", tidy = T, newdata=newdata, t=t, ci=F) ## this gives predictions of relative survival for each individual
  rs$id <- rep(1:dim(newdata)[1],each=length(t))
  names(rs)[names(rs)=="est"] <- "rs"
  pred <- rs %>% left_join(expsurv$expsurv, by=c("id","time"))
  pred$levels.fct <- factor(seq_along(t))
  if(weighted){
    pred$weights <- rep(newdata$`(weights)`, each=length(t))
    predsum <- pred %>% group_by(levels.fct, time) %>%
      summarise("at{i}" := tr.fun(weighted.mean(.data$rs*.data$es, .data$weights))) %>%
      ungroup() %>% select(-levels.fct)
  } else {
    predsum <- pred %>% group_by(levels.fct, time) %>% summarise("at{i}" := tr.fun(mean(.data$rs*.data$es))) %>%
      ungroup() %>% select(-levels.fct)
  }
  predsum
}

# achazard standardisation function, for use with type "achazard"
standsurv.fn.achazard <- function(t, object, newdata, i, weighted, tr.fun, expsurv){
  rs <- summary(object, type = "survival", tidy = T, newdata=newdata, t=t, ci=F) ## this gives predictions of relative survival for each individual
  rs$id <- rep(1:dim(newdata)[1],each=length(t))
  names(rs)[names(rs)=="est"] <- "rs"
  excessh <- summary(object, type = "hazard", tidy = T, newdata=newdata, t=t, ci=F) ## this gives excess hazard
  excessh$id <- rep(1:dim(newdata)[1],each=length(t))
  names(excessh)[names(excessh)=="est"] <- "excessh"
  pred <- rs %>% left_join(excessh, by=c("id", "time")) %>%
    left_join(expsurv$expsurv, by=c("id","time"))
  pred$levels.fct <- factor(seq_along(t))
  if(weighted){
    pred$weights <- rep(newdata$`(weights)`, each=length(t))
    predsum <- pred %>% group_by(levels.fct, time) %>% 
      summarise("at{i}" := tr.fun(weighted.mean(.data$excessh*.data$eh, .data$rs*.data$es * .data$weights))) %>%
      ungroup() %>% select(-levels.fct)
  } else {
    predsum <- pred %>% group_by(levels.fct, time) %>% 
      summarise("at{i}" := tr.fun(weighted.mean(.data$excessh + .data$eh, .data$rs*.data$es))) %>%
      ungroup() %>% select(-levels.fct)
  }
}

standsurv.fn.acrmst <- function(t, object, newdata, i, rmap, ratetable, scale.ratetable, weighted, 
                   tr.fun = tr("none"), n.gauss.quad = 100){
  # Using Gauss-Legendre quadrature
  newdata$id <- 1:dim(newdata)[1]
  gaussxw <- gauss.quad(n.gauss.quad)
  pred <- vector()
  for(j in seq_along(t)){
    if(t[j]==0){
      pred[j] <- 0
    } else {
      scale <- t[j]/2
      points <- scale*(gaussxw$nodes + 1)
      eval_fn <- acsurv.int.fn(points, object, newdata, rmap, ratetable, scale.ratetable, weighted)
      pred[j] <- tr.fun(scale*sum(gaussxw$weights * eval_fn))
    }
  }
  pred
}

# Function to marginal all-cause survival at any given time points average(S*(t)R*(t))
acsurv.int.fn <- function(t1, object, newdata, rmap, ratetable, scale.ratetable, weighted){
  # Relative survival
  rs <- summary(object, type = "survival", tidy = T, newdata=newdata, t=t1, ci=F)
  rs$id <- rep(1:dim(newdata)[1],each=length(t1))
  rs <- rs %>% left_join(newdata, by="id")
  ## Expected survival
  rs$t1.scale <- rs$time * scale.ratetable
  rs$es <- do.call("survexp", list(formula = t1.scale~1, 
                                   rmap = rmap, method="individual.s", 
                                   ratetable = ratetable, data=rs
  ))
  rs$levels.fct <- factor(seq_along(t))
  
  if(weighted){
    rssum <- rs %>% group_by(levels.fct, time) %>% 
      summarise(acsurv = weighted.mean(.data$est*.data$es, .data$weights)) %>%
      ungroup() %>% select(-levels.fct)
  } else {
    rssum <- rs %>% group_by(levels.fct, time) %>% summarise(acsurv = mean(.data$est*.data$es)) %>%
      ungroup() %>% select(-levels.fct)
  }
  rssum$acsurv
}


# acquantile standardisation function, for use with type "acquantile"
standsurv.fn.acquantile <- function(t, object, newdata, i,  rmap, ratetable, scale.ratetable,
                                    weighted, tr.fun, quantiles, interval){
  newdata$id <- 1:dim(newdata)[1]
  quantile.root.fn <- function(t, q, object, newdata, i, weighted){
      acsurv.int.fn(t, object, newdata, rmap, ratetable, scale.ratetable, weighted) - (1-q)
  }
  predsum <- tibble()
  for(q in quantiles){
    root <- uniroot(quantile.root.fn, interval = interval,
                    q=q, object=object, newdata=newdata, i=i, weighted=weighted)$root
    predsum <- predsum %>% bind_rows(tibble(probability = q) %>% mutate("at{i}" := tr.fun(root)))
  }
  predsum
}


tr <- function(trans){
  switch(trans,
         "log"= log,
         "none"= function(x) x,
         "loglog" = function(x) log(-log(1-x)),
         "logit" = qlogis
  )
}

inv.tr <- function(trans){
  switch(trans,
         "log"= exp,
         "none"= function(x) x,
         "loglog" = function(x) 1-exp(-exp(x)),
         "logit" = plogis
  )
}

#' @importFrom dplyr pull
boot.standsurv <- function(object, B, dat, i, t, type, type2, weighted, se, ci, cl, 
                           rawsim, predsum, expsurv, rmap, ratetable, scale.ratetable, 
                           quantiles, interval, n.gauss.quad){
  if(is.null(rawsim))
    rawsim <- attributes(normboot.flexsurvreg(object, B=B, raw=T))$rawsim ## only run this once, not for every specified _at
  
  X <- form.model.matrix(object, as.data.frame(dat), na.action=na.pass)
  if(!(type2 %in% c("acrmst","quantile","acquantile"))){
    if(!(type2 %in% c("acsurvival", "achazard"))){
      sim.pred <- normbootfn.flexsurvreg(object, t=t, start=0, X=X, fn=summary.fns(object, type2), B=B, rawsim=rawsim) # pts, sims, times
    } else {
      sim.pred <- normbootfn.flexsurvreg(object, t=t, start=0, X=X, fn=summary.fns(object, type), B=B, rawsim=rawsim) # pts, sims, times
    }
    
    if(weighted){
      weights <- array(dat$`(weights)`, dim(sim.pred))
    } else {
      weights <- array(1, dim(sim.pred))
    }
    if(type2=="hazard"){
      # Weight individual hazards by survival function to get hazard of the standardized survival
      haz <- sim.pred
      surv <- normbootfn.flexsurvreg(object, t=t, start=0, X=X, fn=summary.fns(object, "survival"), B=B, rawsim=rawsim) # pts, sims, times
      stand.pred <- apply(haz*surv*weights, c(2,3),sum) / apply(surv*weights,c(2,3),sum)
    } else if(type2=="acsurvival"){
      rs <- sim.pred
      es <- array(expsurv$expsurv$es, dim= dim(rs)[c(1,3,2)])
      es <- aperm(es, c(1, 3, 2))
      stand.pred <- apply(es*rs*weights, c(2,3),sum) / apply(weights,c(2,3),sum)
    } else if(type2=="achazard"){
      excessh <- sim.pred
      rs <- normbootfn.flexsurvreg(object, t=t, start=0, X=X, fn=summary.fns(object, "survival"), B=B, rawsim=rawsim) # pts, sims, times
      es <- array(expsurv$expsurv$es, dim= dim(rs)[c(1,3,2)])
      es <- aperm(es, c(1, 3, 2))
      eh <- array(expsurv$expsurv$eh, dim= dim(rs)[c(1,3,2)])
      eh <- aperm(eh, c(1, 3, 2))
      stand.pred <- apply(rs*es*weights*(excessh+eh), c(2,3),sum) / apply(rs*es*weights,c(2,3),sum)
    } else {
      stand.pred <- apply(sim.pred*weights, c(2,3), sum) / apply(weights,c(2,3),sum)
    }
  } 
  if(type2=="acrmst") { ## if type2=="acrmst"
    ## This for now manipulates rawsim within standsurv.fn.acrmst. This code could be made more slick.
    stand.pred <- matrix(nrow=B, ncol=length(t))
    for(b in seq(length.out=B)) {
      newobject <- object
      newobject$res.t[,"est"] <- rawsim[b,]
      newobject$res[newobject$covpars,"est"] <- rawsim[b, newobject$covpars]
      newobject$res[newobject$dlist$pars,"est"] <- NA  ## setting to NA to be safe
      newobject$res.t[,-1] <- newobject$res[,-1] <- NA ## setting to NA to be safe
      ## standardisation (averaging across patients) is done within acrmst itself 
      stand.pred[b,] <- standsurv.fn.acrmst(t, newobject, newdata=dat, i, rmap, ratetable, scale.ratetable, weighted, 
                                            tr.fun=tr("none"), n.gauss.quad)
    }
  }
  if(type2=="quantile"){
    ## This for now manipulates rawsim within standsurv.fn.quantile. This code could be made more slick.
    stand.pred <- matrix(nrow=B, ncol=length(quantiles))
    for(b in seq(length.out=B)) {
      newobject <- object
      newobject$res.t[,"est"] <- rawsim[b,]
      newobject$res[newobject$covpars,"est"] <- rawsim[b, newobject$covpars]
      newobject$res[newobject$dlist$pars,"est"] <- NA  ## setting to NA to be safe
      newobject$res.t[,-1] <- newobject$res[,-1] <- NA ## setting to NA to be safe
      ## standardisation (averaging across patients) is done within standsurv.fn.quantile itself 
      stand.pred[b,] <- pull(standsurv.fn.quantile(newobject, newdata=dat, i, weighted, tr.fun=tr("none"), 
                                              quantiles, interval)[,2])
    }
  }
  if(type2=="acquantile"){
    stop("bootstrap not yet applied with acquantile")
  }
  if(se == TRUE){
    stand.pred.se <- tibble("at{i}_se" := apply(stand.pred, 2, sd, na.rm=TRUE))
    predsum <- predsum %>% bind_cols(stand.pred.se)
  }
  if(ci == TRUE){
    stand.pred.quant <- apply(stand.pred, 2, function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE) )
    stand.pred.quant <- as_tibble(t(stand.pred.quant)) %>% 
      rename("at{i}_lci" := "2.5%", "at{i}_uci" := "97.5%")
    predsum <- predsum %>% bind_cols(stand.pred.quant)
  }
  
  return(list(predsum = predsum, stand.pred=stand.pred, rawsim=rawsim))
}

#' @importFrom numDeriv grad
deltamethod.standsurv <- function(object, newdata, type2, t, i, se, ci, 
                                  predsum, trans, cl, weighted, expsurv,
                                  rmap, ratetable, scale.ratetable, quantiles, interval,
                                  n.gauss.quad){
  if(!(type2 %in% c("quantile","acquantile"))){
    g1 <- function(coef, t, trans) {
      object$res[,"est"] <- object$res.t[,"est"] <- coef
      standsurv.fn(object, type=type2, newdata=newdata, t=t, i=i, trans, 
                   weighted=weighted, expsurv=expsurv, rmap, ratetable, scale.ratetable, 
                   quantiles, interval, n.gauss.quad)[,2,drop=T]
    }
  } else {
    g2 <- function(coef, quantiles, trans) {
      object$res[,"est"] <- object$res.t[,"est"] <- coef
      standsurv.fn(object, type=type2, newdata=newdata, t=t, i=i, trans, 
                   weighted=weighted, expsurv=expsurv, rmap, ratetable, scale.ratetable, 
                   quantiles=quantiles, interval, n.gauss.quad)[,2,drop=T]
    }
  }
  est <- standsurv.fn(object, type=type2, newdata=newdata, t=t, i=i, trans="none", 
                      weighted=weighted, expsurv=expsurv, rmap, ratetable, scale.ratetable, 
                      quantiles, interval, n.gauss.quad)[,2,drop=T]
  
  var.none <- NULL
  if(se==TRUE){
    # Calculate for each value of t the untransformed standardized measure
    if(!(type2 %in% c("quantile","acquantile"))){
      var.none <- sapply(t, function(ti){
        gd <- grad(g1, coef(object), method="simple" ,t=ti, 
                   trans="none")
        gd %*% vcov(object) %*% gd
      })
    } else {
      var.none <- sapply(quantiles, function(q){
        gd <- grad(g2, coef(object), method="simple" ,quantiles=q, 
                   trans="none")
        gd %*% vcov(object) %*% gd
      })
    }
    stand.pred.se <- as_tibble(sqrt(var.none)) %>% rename("at{i}_se" := "value")
    predsum <- predsum %>% bind_cols(stand.pred.se)   
  }  
  if(ci==TRUE){
    # Calculate for each value of t the transformed standardized measure
    if(trans=="none" & !is.null(var.none)){
      var.trans <- var.none ## use already calculated variances
    } else {
      if(!(type2 %in% c("quantile","acquantile"))){
        var.trans <- sapply(t, function(ti){
          gd <- grad(g1, coef(object), method="simple" ,t=ti, 
                     trans=trans)
          gd %*% vcov(object) %*% gd
        })
      } else {
        var.trans <- sapply(quantiles, function(q){
          gd <- grad(g2, coef(object), method="simple" ,quantiles=q, 
                     trans=trans)
          gd %*% vcov(object) %*% gd
        })
      }
    }
    tr.fun <- tr(trans)
    inv.tr.fun <- inv.tr(trans)
    stand.pred.lcl <- as_tibble(inv.tr.fun(tr.fun(est)+qnorm((1-cl)/2, lower.tail=T)*sqrt(var.trans))) %>%
      rename("at{i}_lci" := "value")
    stand.pred.ucl <- as_tibble(inv.tr.fun(tr.fun(est)+qnorm((1-cl)/2, lower.tail=F)*sqrt(var.trans))) %>%
      rename("at{i}_uci" := "value")
    predsum <- predsum %>% bind_cols(stand.pred.lcl, stand.pred.ucl)   
  }
  predsum
}

#' @importFrom numDeriv grad
deltamethod.contrast.standsurv <- function(object, dat, dat.ref,
                               type2, t, i, atreference, se, ci, predsum, 
                               trans.contrast, cl, contrast, weighted, expsurv,
                               rmap, ratetable, scale.ratetable, quantiles, interval,
                               n.gauss.quad){
  tr.fun <- tr(trans.contrast)
  inv.tr.fun <- inv.tr(trans.contrast)
  contrast.fn <- switch(contrast, "difference"= `-`, "ratio"= `/` )
  
  if(!(type2 %in% c("quantile","acquantile"))){
    g1 <- function(coef, t, tr.fun, contrast.fn) {
      object$res[,"est"] <- object$res.t[,"est"] <- coef
      tr.fun(contrast.fn(standsurv.fn(object, type=type2, newdata=dat, t=t, i=i, 
                                      trans="none", weighted=weighted, expsurv=expsurv,
                                      rmap, ratetable, scale.ratetable, quantiles, interval,
                                      n.gauss.quad)[,2,drop=T], 
                         standsurv.fn(object, type=type2, newdata=dat.ref, t=t, i=i, 
                                      trans="none", weighted=weighted, expsurv=expsurv,
                                      rmap, ratetable, scale.ratetable, quantiles, interval,
                                      n.gauss.quad)[,2,drop=T]))
    }
  } else {
    g2 <- function(coef, quantiles, tr.fun, contrast.fn) {
      object$res[,"est"] <- object$res.t[,"est"] <- coef
      tr.fun(contrast.fn(standsurv.fn(object, type=type2, newdata=dat, t=t, i=i, 
                                      trans="none", weighted=weighted, expsurv=expsurv,
                                      rmap, ratetable, scale.ratetable, quantiles=quantiles, 
                                      interval, n.gauss.quad)[,2,drop=T], 
                         standsurv.fn(object, type=type2, newdata=dat.ref, t=t, i=i, 
                                      trans="none", weighted=weighted, expsurv=expsurv,
                                      rmap, ratetable, scale.ratetable, quantiles=quantiles, 
                                      interval, n.gauss.quad)[,2,drop=T]))
    }
  }
  
  est <- contrast.fn(standsurv.fn(object, type=type2, newdata=dat, t=t, i=i, 
                                  trans="none", weighted=weighted, expsurv=expsurv,
                                  rmap, ratetable, scale.ratetable, quantiles, interval,
                                  n.gauss.quad)[,2,drop=T],
                     standsurv.fn(object, type=type2, newdata=dat.ref, t=t, i=i, 
                                  trans="none", weighted=weighted, expsurv=expsurv,
                                  rmap, ratetable, scale.ratetable, quantiles, interval,
                                  n.gauss.quad)[,2,drop=T])
  stand.pred <- as_tibble(est) %>%
    rename("contrast{i}_{atreference}" := "value")
  predsum <- predsum %>% bind_cols(stand.pred)   
  
  var.none <- NULL
  if(se==TRUE){
    # Calculate for each value of t the untransformed standardized measure
    if(!(type2 %in% c("quantile","acquantile"))){
      var.none <- sapply(t, function(ti){
        gd <- grad(g1, coef(object), method="simple" ,t=ti, 
                   tr.fun=function(x) x, contrast.fn=contrast.fn)
        gd %*% vcov(object) %*% gd
      })
    } else {
      var.none <- sapply(quantiles, function(q){
        gd <- grad(g2, coef(object), method="simple" ,quantiles=q, 
                   tr.fun=function(x) x, contrast.fn=contrast.fn)
        gd %*% vcov(object) %*% gd
      })
    }
    stand.pred.se <- as_tibble(sqrt(var.none)) %>% rename("contrast{i}_{atreference}_se" := "value") 
    predsum <- predsum %>% bind_cols(stand.pred.se)   
  }  
  if(ci==TRUE){
    # Calculate for each value of t the transformed standardized measure
    if(trans.contrast=="none" & !is.null(var.none)){
      var.trans <- var.none ## use already calculated variances
    } else {
      if(!(type2 %in% c("quantile","acquantile"))){
        var.trans <- sapply(t, function(ti){
          gd <- grad(g1, coef(object), method="simple" ,t=ti, 
                     tr.fun=tr.fun, contrast.fn=contrast.fn)
          gd %*% vcov(object) %*% gd
        })
      } else {
        var.trans <- sapply(quantiles, function(q){
          gd <- grad(g2, coef(object), method="simple" ,quantiles=q, 
                     tr.fun=tr.fun, contrast.fn=contrast.fn)
          gd %*% vcov(object) %*% gd
        })
      }
    }
    stand.pred.lcl <- as_tibble(inv.tr.fun(tr.fun(est)+qnorm((1-cl)/2, lower.tail=T)*sqrt(var.trans))) %>%
      rename("contrast{i}_{atreference}_lci" := "value")
    stand.pred.ucl <- as_tibble(inv.tr.fun(tr.fun(est)+qnorm((1-cl)/2, lower.tail=F)*sqrt(var.trans))) %>%
      rename("contrast{i}_{atreference}_uci" := "value")
    predsum <- predsum %>% bind_cols(stand.pred.lcl, stand.pred.ucl)   
  }
  predsum
}
    


#' Tidy a standsurv object. 
#' 
#' This function is used internally by \code{standsurv} and tidy
#' data.frames are automatically returned by the function.
#'
#' @param x A standsurv object.
#' @param ... Not currently used.
#'
#' @return Returns additional tidy data.frames (tibbles)
#' stored as attributes named standpred_at and standpred_contrast.
#' @importFrom dplyr select
#' @importFrom dplyr inner_join
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect matches
#' @export
#'
tidy.standsurv <- function(x, ...){
  standpred <- x
  at <-attributes(standpred)$at
  atreference <- attributes(standpred)$atreference
  type <- attributes(standpred)$type
  label <- attributes(standpred)$label
  contrast <- attributes(standpred)$contrast
  ci <- any(grepl("_lci",names(standpred)))
  by <- ifelse(type == "quantile", "probability", "time")
  
  class(standpred) <- class(standpred)[class(standpred)!="standsurv"]

  standpred_at <- standpred %>% 
    select(c(all_of(by), matches("at[0-9]+$"))) %>%
    pivot_longer(cols=matches("at[0-9]+$"),
                 names_to = "at",
                 values_to = type,
                 names_prefix = "at")
  if(ci){
    standpred_at_lci <- standpred %>% 
      select(c(all_of(by), matches("at[0-9]+_lci"))) %>%
      pivot_longer(cols=matches("at[0-9]+_lci"),
                   names_to = "at",
                   names_pattern = "at(.+)_lci",
                   values_to = paste0(type,"_lci")) %>%
      distinct()
    standpred_at_uci <- standpred %>% 
      select(c(all_of(by), matches("at[0-9]+_uci"))) %>%
      pivot_longer(cols=matches("at[0-9]+_uci"),
                   names_to = "at",
                   names_pattern = "at(.+)_uci",
                   values_to = paste0(type,"_uci")) %>%
      distinct()
    standpred_at <- standpred_at %>% inner_join(standpred_at_lci, by=c(by,"at")) %>%
      inner_join(standpred_at_uci, by=c(by, "at"))
  }
  for(i in 1:length(at)){
    standpred_at <- standpred_at %>% 
      mutate(at = replace(at, at==i,
                          label[i]))
  }  
  attr(standpred,"standpred_at") <- standpred_at
  if(!is.null(contrast)){
    ## Contrast numbers
    cnums <- (1:length(at))[-atreference]
    standpred_contrast <- standpred %>% 
      select(c(all_of(by), matches("contrast[0-9]+_[0-9]+$"))) %>%
      pivot_longer(cols=matches("contrast[0-9]+_[0-9]+$"),
                   names_to = "contrast",
                   values_to = contrast,
                   names_prefix = "contrast")
    if(ci){
      standpred_contrast_lci <- standpred %>% 
        select(c(all_of(by), matches("contrast[0-9]+_[0-9]+_lci"))) %>%
        pivot_longer(cols=matches("contrast[0-9]+_[0-9]+_lci"),
                     names_to = "contrast",
                     names_pattern = "contrast(.+)_lci",
                     values_to = paste0(contrast,"_lci"))
      standpred_contrast_uci <- standpred %>% 
        select(c(all_of(by), matches("contrast[0-9]+_[0-9]+_uci"))) %>%
        pivot_longer(cols=matches("contrast[0-9]+_[0-9]+_uci"),
                     names_to = "contrast",
                     names_pattern = "contrast(.+)_uci",
                     values_to = paste0(contrast,"_uci"))
      standpred_contrast <- standpred_contrast %>% inner_join(standpred_contrast_lci, by=c(by, "contrast")) %>%
        inner_join(standpred_contrast_uci, by=c(by, "contrast"))
    }
    for(i in cnums){
      standpred_contrast <- standpred_contrast %>% 
        mutate(contrast = replace(contrast, contrast==paste0(i,"_",atreference),
                                  paste0(label[i]," vs ",label[atreference])))
    }  
    attr(standpred,"standpred_contrast") <- standpred_contrast
  }
  class(standpred) <- c("standsurv",class(standpred))
  standpred
}


#' Plot standardized metrics from a fitted flexsurv model
#'
#' Plot standardized metrics such as the marginal survival, restricted mean 
#' survival and hazard, based on a fitted flexsurv model.
#' 
#' @param x A standsurv object returned by \code{standsurv}
#' @param contrast Should contrasts of standardized metrics be plotted. Defaults
#' to FALSE
#' @param ci Should confidence intervals be plotted (if calculated in 
#' \code{standsurv})? 
#' @param expected Should the marginal expected survival / hazard also be 
#' plotted? This can only be invoked if \code{rmap} and \code{ratetable} have 
#' been passed to \code{standsurv}
#' @param ... Not currently used
#'
#' @return A ggplot showing the standardized metric calculated by 
#' \code{standsurv} over time. Modification of the plot is
#' possible by adding further ggplot objects, see Examples.
#' @import ggplot2
#' @export
#'
#' @examples
#'## Use bc dataset, with an age variable appended
#'## mean age is higher in those with smaller observed survival times 
#'newbc <- bc
#'newbc$age <- rnorm(dim(bc)[1], mean = 65-scale(newbc$recyrs, scale=FALSE), 
#' sd = 5)
#'
#'## Fit a Weibull flexsurv model with group and age as covariates
#'weib_age <- flexsurvreg(Surv(recyrs, censrec) ~ group+age, data=newbc,
#'                        dist="weibull")
#'## Calculate standardized survival and the difference in standardized survival
#'## for the three levels of group across a grid of survival times
#'standsurv_weib_age <- standsurv(weib_age,
#'                                            at = list(list(group="Good"),
#'                                                      list(group="Medium"),
#'                                                      list(group="Poor")),
#'                                            t=seq(0,7, length=100),
#'                                            contrast = "difference", ci=TRUE,
#'                                            boot = TRUE, B=10, seed=123)
#'plot(standsurv_weib_age)
#'plot(standsurv_weib_age) + ggplot2::theme_bw() + ggplot2::ylab("Survival") +
#'  ggplot2::xlab("Time (years)") + 
#'  ggplot2::guides(color=ggplot2::guide_legend(title="Prognosis"),
#'                                fill=ggplot2::guide_legend(title="Prognosis"))
#'plot(standsurv_weib_age, contrast=TRUE, ci=TRUE) + 
#'  ggplot2::ylab("Difference in survival") 
plot.standsurv <- function(x, contrast = FALSE, ci = FALSE, expected = FALSE, ...){
  if(!contrast){
    obj <- attributes(x)$standpred_at
    obj <- obj %>% mutate(Population = "Study")
    y <- attributes(x)$type
    if(y=="quantile") stop("plot method not configured for type='quantile'")
    group <- "at"
    if(expected){
      if(is.null(attributes(x)$expected)) 
        stop("Expected survival/hazards have not been calculated by standsurv")
      if(!(y %in% c("hazard", "achazard", "survival", "acsurvival")))
        stop(paste0("Expected survival/hazards cannot be plotted with type = ",y))
      if(y %in% c("hazard", "achazard"))
        y2 <- "exphaz"
      if(y %in% c("survival", "acsurvival"))
        y2 <- "expsurv"
      obj2 <- attributes(x)$expected[,c("time",y2)] %>%
        rename({{y}} := all_of(y2)) %>% 
        mutate({{group}} := "Expected", Population = "Expected")
      obj <- obj %>% bind_rows(obj2)
    }
  } else {
    if(is.null(attributes(x)$standpred_contrast)) 
      stop("Contrasts have not been calculated by standsurv")
    obj <- attributes(x)$standpred_contrast
    obj <- obj %>% mutate(Population = "Study")
    y <- attributes(x)$contrast
    group <- "contrast"
    if(expected) stop("Expected survival/hazard cannot be plotted with contrasts")
  }
  linetype <- c("Study" = "solid", "Expected" = "dashed")
  p <- ggplot() + 
    geom_line(aes(x=.data[["time"]], y=.data[[y]], color=.data[[group]],
                   linetype= .data[["Population"]]), data=obj) + 
    xlab("Time") +
    scale_linetype_manual(values = linetype, guide="none") 
    
  if(ci){
    if(any(grepl("_lci",names(obj)))){
      p <- p + geom_ribbon(aes(x= .data[["time"]],
                               ymin=.data[[paste0(y,"_lci")]],
                               ymax=.data[[paste0(y,"_uci")]],
                              fill= .data[[group]]), data = obj, alpha=0.2)
    } else warning("Confidence intervals have not been calculated in standsurv. None will be plotted")
  }
  p
}
              
#' @importFrom dplyr bind_rows
expsurv.fn <- function(t, rmap, ratetable, data, weighted, scale.ratetable){
  expsurv <- exphaz <- tibble()
  for(l in 1:length(t)){
    data$t.temp <- t[l] * scale.ratetable
    es <- do.call("survexp", list(formula = t.temp~1, 
                                        rmap = rmap, method="individual.s", 
                                        ratetable = ratetable, data = data  
                                        ))
    expsurv <- expsurv %>% bind_rows(tibble(time = t[l], es = es, id = 1:dim(data)[1]))
    # individual hazards from difference in cumulative hazards / epsilon
    epsilon <- 0.001
    data$t.temp.epsilon <- data$t.temp + epsilon
    eh.epsilon <- do.call("survexp", list(formula = t.temp.epsilon~1, 
                                  rmap = rmap, method="individual.h", 
                                  ratetable = ratetable, data = data  
    ))
    eh.epsilon0 <- do.call("survexp", list(formula = t.temp~1, 
                                           rmap = rmap, method="individual.h", 
                                           ratetable = ratetable, data=data
    ))
    eh <- scale.ratetable * (eh.epsilon -  eh.epsilon0) / epsilon
    exphaz <- exphaz %>% bind_rows(tibble(time = t[l], eh=eh, id=1:dim(data)[1]))
  }
  expsurv <- expsurv %>% inner_join(exphaz, by= c("id", "time"))
  if(weighted){
    expsurv$weights <- rep(data$`(weights)`, each=length(t))
    marginal <- expsurv %>% group_by(time) %>% 
      summarise("expsurv" = weighted.mean(.data$es, .data$weights),
                "exphaz" = weighted.mean(.data$eh, .data$es * .data$weights))
  } else {
    marginal <- expsurv %>% group_by(time) %>% 
      summarise("expsurv" = mean(.data$es), 
                "exphaz" = weighted.mean(.data$eh, .data$es))
  }
  return(list(marginal=marginal, expsurv=expsurv))
}

## to remove R CMD check NOTE
utils::globalVariables(c("levels.fct"))
