## sinh(log(y))
logh <- function(x) { 0.5 * (x - 1/x) }

buildTransformer <- function(inits, nbpars, dlist) {
    par.transform <-
        lapply(seq_len(nbpars),
               function(ind) {
                   xform  <- dlist$inv.transforms[[ind]]
                   function(pars) {
                       xform(pars[[ind]])
                   }
               })
    names(par.transform) <- names(inits)[seq_len(nbpars)]
    function(pars) {
        lapply(par.transform,
               function(item, par) { item(par) },
               pars)
    }
}

buildAuxParms <- function(aux, dlist) {
    aux.transform <- list()
    for (ind in seq_along(aux)) {
        name <- names(aux)[[ind]]
        if (!(name %in% dlist$pars)) {
            aux.transform[[name]] <- aux[[ind]]
        }
    }
    aux.transform
}

logLikFactory <- function(Y, X=0, weights, bhazard, dlist,
                          inits, dfns, aux, mx, fixedpars=NULL) {
    pars   <- inits
    npars  <- length(pars)
    nbpars <- length(dlist$pars)
    insert.locations <- setdiff(seq_len(npars),
                                fixedpars)
    
    ## which are the subjects with events
    event <- Y[,"status"] == 1
    event.times <- Y[event, "time1"]
    left.censor <- Y[!event, "time2"]
    right.censor <- Y[!event, "time1"]

    event.weights <- weights[event]
    no.event.weights <- weights[!event]
    
    par.transform <- buildTransformer(inits, nbpars, dlist)

    aux.pars <- buildAuxParms(aux, dlist)

    default.offset <- rep.int(0, length(event.times))
    do.hazard <- any(bhazard > 0)

    loglik <- rep.int(0, nrow(Y))
    ## the ... here is to work around optim
    function(optpars, ...) {

        pars[insert.locations] <- optpars
        raw.pars <- pars
        pars <- as.list(pars)

        pars.event <- pars.nevent <- pars
        if (npars > nbpars) {
            beta <- raw.pars[(nbpars+1):npars]
            for (i in dlist$pars){
                pars[[i]] <- pars[[i]] +
                    X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
                pars.event[[i]] <- pars[[i]][event]
                pars.nevent[[i]] <- pars[[i]][!event]
            }
        }

        fnargs <- c(par.transform(pars),
                    aux.pars)
        fnargs.event <- c(par.transform(pars.event),
                          aux.pars)
        fnargs.nevent <- c(par.transform(pars.nevent),
                           aux.pars)
        
        ## Generic survival model likelihood contributions
        ## Observed deaths
        dargs <- fnargs.event
        dargs$x <- event.times
        dargs$log <- TRUE
        logdens <- do.call(dfns$d, dargs)
        
        ## Left censoring times (upper bound for event time) 
        if (any(!event)){
            pmaxargs <- fnargs.nevent
            pmaxargs$q <- left.censor # Inf if right-censored, giving pmax=1
            pmax <- do.call(dfns$p, pmaxargs)
            pmax[pmaxargs$q==Inf] <- 1  # in case user-defined function doesn't already do this
            
        ## Right censoring times (lower bound for event time) 
            pargs <- fnargs.nevent
            pargs$q <- right.censor
            pmin <- do.call(dfns$p, pargs)
        }
        
        ## Left-truncation
        targs   <- fnargs
        targs$q <- Y[,"start"]
        pobs <- 1 - do.call(dfns$p, targs) # prob of being observed = 1 unless left-truncated

        ## Hazard offset for relative survival models
        if (do.hazard){
            pargs   <- fnargs.event
            pargs$q <- event.times
            pminb   <- do.call(dfns$p, pargs)
            loghaz  <- logdens - log(1 - pminb)
            offseti <- log(1 + bhazard[event] / exp(loghaz)*weights[event])
        } else {
            offseti <- default.offset
        }
        ## Express as vector of individual likelihood contributions
        loglik[event] <- (logdens*event.weights) + offseti
        if (any(!event))
            loglik[!event] <- (log(pmax - pmin)*no.event.weights)
        loglik <- loglik - log(pobs)*weights
        
        ret <- -sum(loglik)
        attr(ret, "indiv") <- loglik
        ret
    }
}

minusloglik.flexsurv <- function(optpars, Y, X=0, weights, bhazard,
                                 dlist, inits,
                                 dfns, aux, mx, fixedpars=NULL) {
    logLikFactory(Y, X, weights, bhazard, dlist, inits,
                  dfns, aux, mx, fixedpars=fixedpars)(optpars)
}

check.dlist <- function(dlist){
    ## put tests in testthat
    if (is.null(dlist$name)) stop("\"name\" element of custom distribution list not given")
    if (!is.character(dlist$name)) stop("\"name\" element of custom distribution list should be a string")
    if (is.null(dlist$pars)) stop("parameter names \"pars\" not given in custom distribution list")
    if (!is.character(dlist$pars)) stop("parameter names \"pars\" should be a character vector")
    npars <- length(dlist$pars)
    if (is.null(dlist$location)) {
        warning("location parameter not given, assuming it is the first one")
        dlist$location <- dlist$pars[1]
    }
    if (!(dlist$location %in% dlist$pars)) {
        stop(sprintf("location parameter \"%s\" not in list of parameters", dlist$location))
    }
    if (is.null(dlist$transforms)) stop("transforms not given in custom distribution list")
    if (is.null(dlist$inv.transforms)) stop("inverse transforms not given in custom distribution list")
    if (!is.list(dlist$transforms)) stop("\"transforms\" must be a list of functions")
    if (!is.list(dlist$inv.transforms)) stop("\"inv.transforms\" must be a list of functions")
    if (!all(sapply(dlist$transforms, is.function))) stop("some of \"transforms\" are not functions")
    if (!all(sapply(dlist$inv.transforms, is.function))) stop("some of \"inv.transforms\" are not functions")
    if (length(dlist$transforms) != npars) stop("transforms vector of length ",length(dlist$transforms),", parameter names of length ",npars)
    if (length(dlist$inv.transforms) != npars) stop("inverse transforms vector of length ",length(dlist$inv.transforms),", parameter names of length ",npars) #
    for (i in 1:npars){
        if (is.character(dlist$transforms[[i]])) dlist$transforms[[i]] <- get(dlist$transforms[[i]])
        if (is.character(dlist$inv.transforms[[i]])) dlist$inv.transforms[[i]] <- get(dlist$inv.transforms[[i]])
        if (!is.function(dlist$transforms[[i]])) stop("Transformation function for parameter ", i, " not defined")
        if (!is.function(dlist$inv.transforms[[i]])) stop("Inverse transformation function for parameter ", i, " not defined")
    }
    if (!is.null(dlist$inits) && !is.function(dlist$inits)) stop("\"inits\" element of custom distribution list must be a function")
    dlist
}
## Return formula for linear model on parameter called "par"
## Parameters should not have the same name as anything that might
## appear as a function in a formula (such as "I", "strata", or
## "factor").  If any parameters of the distribution being used are
## named like this, then such model functions will be interpreted as
## parameters and will not work

check.formula <- function(formula, dlist){
    if (!inherits(formula,"formula")) stop("\"formula\" must be a formula object")
    if (!("strata" %in% dlist$pars)){
        labs <- attr(terms(formula), "term.labels")
        strat <- grep("strata\\((.+)\\)",labs)
        if (any(strat)){
            cov <- gsub("strata\\((.+)\\)","\\1",labs[strat[1]])
            warning("Ignoring \"strata\" function: interpreting \"",cov, "\" as a covariate on \"", dlist$location, "\"")
        }
    }
}

ancpar.formula <- function(formula, par){
    labs <- attr(terms(formula), "term.labels")
    pattern <- paste0(par,"\\((.+)\\)")
    labs <- grep(pattern,labs,value=TRUE)
    if (length(labs)==0) return(NULL)
    labs <- gsub(pattern, "\\1", labs)
    f2 <- reformulate(labs)
    environment(f2) <- environment(formula)
    f2
}

## Omit formula terms containing ancillary parameters, leaving only
## the formula for the location parameter

get.locform <- function(formula, ancnames){
    labs <- attr(terms(formula), "term.labels")
    dropx <- unlist(lapply(ancnames, function(x){grep(paste0(x,"\\((.+)\\)"),labs)}))
    formula(terms(formula)[c(0,setdiff(seq_along(labs),dropx))])
}

## Concatenate location formula (that includes Surv response term)
## with list of ancillary formulae, giving a merged formula to obtain
## the model frame

concat.formulae <- function(formula,forms){
    covnames <- unlist(lapply(forms, function(x)attr(terms(x),"term.labels")))
    covform <- if(length(covnames)==0) "1" else paste(covnames, collapse=" + ")
    respname <- as.character(formula[2])
    form <- paste0(respname, " ~ ", covform)
    f2 <- eval(parse(text=form))
    environment(f2) <- environment(formula)
    ## names of variables in the data, not the formula, with functions such as factor() stripped
    ## used for error message with incomplete "newdata" in summary()
    covnames.bare <- unlist(lapply(forms, function(x)all.vars(delete.response(terms(x)))))
    attr(f2, "covnames") <- covnames.bare
    attr(f2, "covnames.orig") <- covnames
    f2
}

## User-supplied initial value functions don't have to include all
## four possible arguments: this expands them if they don't

expand.inits.args <- function(inits){
    inits2 <- inits
    formals(inits2) <- alist(t=,mf=,mml=,aux=)
    body(inits2) <- body(inits)
    inits2
}

## User-supplied summary output functions don't have to include all
## two possible arguments: this expands them if they don't

expand.summfn.args <- function(summfn){
    summfn2 <- summfn
    args <- c(alist(t=,start=), formals(summfn))
    formals(summfn2) <- args[!duplicated(names(args))]
    body(summfn2) <- body(summfn)
    summfn2
}

### On entry:
### event (status=1)            time1=event time
### right-censoring (status=0)  time1=lower bound
### left-censoring (status=2)   time1=upper bound 
### interval-censoring (status=3)  time1=lower, time2=upper

### On exit
### time1=lower bound or event time
### time2=upper bound
### start=left truncation time
### so meaning of time1,time2 reversed with left-censoring

check.flexsurv.response <- function(Y){
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
### convert Y from Surv object to numeric matrix
### though "time" only used for initial values, printed time at risk, empirical hazard
    if (attr(Y, "type") == "counting")
        Y <- cbind(Y, time=Y[,"stop"] - Y[,"start"], time1=Y[,"stop"], time2=Inf)
    else if (attr(Y, "type") == "interval"){
        Y[,"time2"][Y[,"status"]==0] <- Inf   # upper bound with right censoring 
        Y[,"time2"][Y[,"status"]==2] <- Y[,"time1"][Y[,"status"]==2]
        Y[,"time1"][Y[,"status"]==2] <- 0  # 
        Y <- cbind(Y, start=0, stop=Y[,"time1"], time=Y[,"time1"])
    }
    else if (attr(Y, "type") == "right")
        Y <- cbind(Y, start=0, stop=Y[,"time"], time1=Y[,"time"], time2=Inf)
    else stop("Survival object type \"", attr(Y, "type"), "\"", " not supported")
    if (any(Y[,"time1"]<0)){
        stop("Negative survival times in the data")
    }
    Y
}

compress.model.matrices <- function(mml){
    cbind.drop.intercept <- function(...)do.call("cbind", lapply(list(...), function(x)x[,-1,drop=FALSE]))
    X <- do.call("cbind.drop.intercept",mml)
    loc.cnames <- colnames(mml[[1]])[-1]
    anc.cnames <- unlist(mapply(function(x,y)sprintf("%s(%s)",x,y), names(mml[-1]), lapply(mml[-1], function(x)colnames(x)[-1])))
    cnames <- c(loc.cnames, anc.cnames)
    colnames(X) <- cnames
    X
}



##' Flexible parametric regression for time-to-event data
##' 
##' Parametric modelling or regression for time-to-event data.  Several
##' built-in distributions are available, and users may supply their own.
##' 
##' Parameters are estimated by maximum likelihood using the algorithms
##' available in the standard R \code{\link{optim}} function.  Parameters
##' defined to be positive are estimated on the log scale.  Confidence
##' intervals are estimated from the Hessian at the maximum, and transformed
##' back to the original scale of the parameters.
##' 
##' The usage of \code{\link{flexsurvreg}} is intended to be similar to
##' \code{\link[survival]{survreg}} in the \pkg{survival} package.
##' 
##' @aliases flexsurvreg flexsurv.dists
##' @param formula A formula expression in conventional R linear modelling
##' syntax. The response must be a survival object as returned by the
##' \code{\link{Surv}} function, and any covariates are given on the right-hand
##' side.  For example,
##' 
##' \code{Surv(time, dead) ~ age + sex}
##' 
##' \code{Surv} objects of \code{type="right"},\code{"counting"},
##' \code{"interval1"} or \code{"interval2"} are supported, corresponding to
##' right-censored, left-truncated or interval-censored observations.
##' 
##' If there are no covariates, specify \code{1} on the right hand side, for
##' example \code{Surv(time, dead) ~ 1}.
##' 
##' By default, covariates are placed on the ``location'' parameter of the
##' distribution, typically the "scale" or "rate" parameter, through a linear
##' model, or a log-linear model if this parameter must be positive.  This
##' gives an accelerated failure time model or a proportional hazards model
##' (see \code{dist} below) depending on how the distribution is parameterised.
##' 
##' Covariates can be placed on other (``ancillary'') parameters by using the
##' name of the parameter as a ``function'' in the formula.  For example, in a
##' Weibull model, the following expresses the scale parameter in terms of age
##' and a treatment variable \code{treat}, and the shape parameter in terms of
##' sex and treatment.
##' 
##' \code{Surv(time, dead) ~ age + treat + shape(sex) + shape(treat)}
##' 
##' However, if the names of the ancillary parameters clash with any real
##' functions that might be used in formulae (such as \code{I()}, or
##' \code{factor()}), then those functions will not work in the formula.  A
##' safer way to model covariates on ancillary parameters is through the
##' \code{anc} argument to \code{\link{flexsurvreg}}.
##' 
##' \code{\link{survreg}} users should also note that the function
##' \code{strata()} is ignored, so that any covariates surrounded by
##' \code{strata()} are applied to the location parameter.
##' @param anc An alternative and safer way to model covariates on ancillary
##' parameters, that is, parameters other than the main location parameter of
##' the distribution.  This is a named list of formulae, with the name of each
##' component giving the parameter to be modelled.  The model above can also be
##' defined as:
##' 
##' \code{Surv(time, dead) ~ age + treat, anc = list(shape = ~ sex + treat)}
##' @param data A data frame in which to find variables supplied in
##' \code{formula}.  If not given, the variables should be in the working
##' environment.
##' @param weights Optional variable giving case weights.
##' @param bhazard Optional variable giving expected hazards for relative
##' survival models.
##' @param subset Vector of integers or logicals specifying the subset of the
##' observations to be used in the fit.
##' @param na.action a missing-data filter function, applied after any 'subset'
##' argument has been used. Default is \code{options()$na.action}.
##' @param dist Typically, one of the strings in the first column of the
##' following table, identifying a built-in distribution.  This table also
##' identifies the location parameters, and whether covariates on these
##' parameters represent a proportional hazards (PH) or accelerated failure
##' time (AFT) model.  In an accelerated failure time model, the covariate
##' speeds up or slows down the passage of time.  So if the coefficient
##' (presented on the log scale) is log(2), then doubling the covariate value
##' would give half the expected survival time.
##' 
##' \tabular{llll}{ \code{"gengamma"} \tab Generalized gamma (stable) \tab mu
##' \tab AFT \cr \code{"gengamma.orig"} \tab Generalized gamma (original) \tab
##' scale \tab AFT \cr \code{"genf"} \tab Generalized F (stable) \tab mu \tab
##' AFT \cr \code{"genf.orig"} \tab Generalized F (original) \tab mu \tab AFT
##' \cr \code{"weibull"} \tab Weibull \tab scale \tab AFT \cr \code{"gamma"}
##' \tab Gamma \tab rate \tab AFT \cr \code{"exp"} \tab Exponential \tab rate
##' \tab PH \cr \code{"llogis"} \tab Log-logistic \tab scale \tab AFT \cr
##' \code{"lnorm"} \tab Log-normal \tab meanlog \tab AFT \cr \code{"gompertz"}
##' \tab Gompertz \tab rate \tab PH \cr }
##' 
##' \code{"exponential"} and \code{"lognormal"} can be used as aliases for
##' \code{"exp"} and \code{"lnorm"}, for compatibility with
##' \code{\link{survreg}}.
##' 
##' Alternatively, \code{dist} can be a list specifying a custom distribution.
##' See section ``Custom distributions'' below for how to construct this list.
##' 
##' Very flexible spline-based distributions can also be fitted with
##' \code{\link{flexsurvspline}}.
##' 
##' The parameterisations of the built-in distributions used here are the same
##' as in their built-in distribution functions: \code{\link{dgengamma}},
##' \code{\link{dgengamma.orig}}, \code{\link{dgenf}},
##' \code{\link{dgenf.orig}}, \code{\link{dweibull}}, \code{\link{dgamma}},
##' \code{\link{dexp}}, \code{\link{dlnorm}}, \code{\link{dgompertz}},
##' respectively.  The functions in base R are used where available, otherwise,
##' they are provided in this package.
##' 
##' For the Weibull, exponential and log-normal distributions,
##' \code{\link{flexsurvreg}} simply works by calling \code{\link{survreg}} to
##' obtain the maximum likelihood estimates, then calling \code{\link{optim}}
##' to double-check convergence and obtain the covariance matrix for
##' \code{\link{flexsurvreg}}'s preferred parameterisation.
##' 
##' The Weibull parameterisation is different from that in
##' \code{\link[survival]{survreg}}, instead it is consistent with
##' \code{\link{dweibull}}.  The \code{"scale"} reported by
##' \code{\link[survival]{survreg}} is equivalent to \code{1/shape} as defined
##' by \code{\link{dweibull}} and hence \code{\link{flexsurvreg}}.  The first
##' coefficient \code{(Intercept)} reported by \code{\link[survival]{survreg}}
##' is equivalent to \code{log(scale)} in \code{\link{dweibull}} and
##' \code{\link{flexsurvreg}}.
##' 
##' Similarly in the exponential distribution, the rate, rather than the mean,
##' is modelled on covariates.
##' 
##' The object \code{flexsurv.dists} lists the names of the built-in
##' distributions, their parameters, location parameter, functions used to
##' transform the parameter ranges to and from the real line, and the functions
##' used to generate initial values of each parameter for estimation.
##' @param inits An optional numeric vector giving initial values for each
##' unknown parameter.  These are numbered in the order: baseline parameters
##' (in the order they appear in the distribution function, e.g. shape before
##' scale in the Weibull), covariate effects on the location parameter,
##' covariate effects on the remaining parameters.  This is the same order as
##' the printed estimates in the fitted model.
##' 
##' If not specified, default initial values are chosen from a simple summary
##' of the survival or censoring times, for example the mean is often used to
##' initialize scale parameters.  See the object \code{flexsurv.dists} for the
##' exact methods used.  If the likelihood surface may be uneven, it is advised
##' to run the optimisation starting from various different initial values to
##' ensure convergence to the true global maximum.
##' @param fixedpars Vector of indices of parameters whose values will be fixed
##' at their initial values during the optimisation.  The indices are ordered
##' as in \code{inits}.  For example, in a stable generalized Gamma model with
##' two covariates, to fix the third of three generalized gamma parameters (the
##' shape \code{Q}, see the help for \code{\link{GenGamma}}) and the second
##' covariate, specify \code{fixedpars = c(3, 5)}
##' @param dfns An alternative way to define a custom survival distribution
##' (see section ``Custom distributions'' below).  A list whose components may
##' include \code{"d"}, \code{"p"}, \code{"h"}, or \code{"H"} containing the
##' probability density, cumulative distribution, hazard, or cumulative hazard
##' functions of the distribution.  For example,
##' 
##' \code{list(d=dllogis, p=pllogis)}.
##' 
##' If \code{dfns} is used, a custom \code{dlist} must still be provided, but
##' \code{dllogis} and \code{pllogis} need not be visible from the global
##' environment.  This is useful if \code{flexsurvreg} is called within other
##' functions or environments where the distribution functions are also defined
##' dynamically.
##' @param aux A named list of other arguments to pass to custom distribution
##' functions.  This is used, for example, by \code{\link{flexsurvspline}} to
##' supply the knot locations and modelling scale (e.g. hazard or odds).  This
##' cannot be used to fix parameters of a distribution --- use \code{fixedpars}
##' for that.
##' @param cl Width of symmetric confidence intervals for maximum likelihood
##' estimates, by default 0.95.
##' @param integ.opts List of named arguments to pass to
##' \code{\link{integrate}}, if a custom density or hazard is provided without
##' its cumulative version.  For example,
##' 
##' \code{integ.opts = list(rel.tol=1e-12)}
##' @param sr.control For the models which use \code{\link{survreg}} to find
##' the maximum likelihood estimates (Weibull, exponential, log-normal), this
##' list is passed as the \code{control} argument to \code{\link{survreg}}.
##' @param ... Optional arguments to the general-purpose optimisation routine
##' \code{\link{optim}}.  For example, the BFGS optimisation algorithm is the
##' default in \code{\link{flexsurvreg}}, but this can be changed, for example
##' to \code{method="Nelder-Mead"} which can be more robust to poor initial
##' values.  If the optimisation fails to converge, consider normalising the
##' problem using, for example, \code{control=list(fnscale = 2500)}, for
##' example, replacing 2500 by a number of the order of magnitude of the
##' likelihood. If 'false' convergence is reported with a non-positive-definite
##' Hessian, then consider tightening the tolerance criteria for convergence.
##' If the optimisation takes a long time, intermediate steps can be printed
##' using the \code{trace} argument of the control list. See
##' \code{\link{optim}} for details.
##' @return A list of class \code{"flexsurvreg"} containing information about
##' the fitted model.  Components of interest to users may include:
##' \item{call}{A copy of the function call, for use in post-processing.}
##' \item{dlist}{List defining the survival distribution used.}
##' \item{res}{Matrix of maximum likelihood estimates and confidence limits,
##' with parameters on their natural scales.} \item{res.t}{Matrix of maximum
##' likelihood estimates and confidence limits, with parameters all transformed
##' to the real line.  The \code{\link{coef}}, \code{\link{vcov}} and
##' \code{\link{confint}} methods for \code{flexsurvreg} objects work on this
##' scale.} \item{coefficients}{The transformed maximum likelihood estimates,
##' as in \code{res.t}. Calling \code{coef()} on a \code{\link{flexsurvreg}}
##' object simply returns this component.} \item{loglik}{Log-likelihood. This will differ from Stata, where the sum of
##' the log uncensored survival times is added to the log-likelihood in
##' survival models, to remove dependency on the time scale.}
##' \item{logliki}{Vector of individual contributions to the log-likelihood}
##' \item{AIC}{Akaike's information criterion (-2*log likelihood + 2*number of
##' estimated parameters)} \item{cov}{Covariance matrix of the parameters, on
##' the real-line scale (e.g. log scale), which can be extracted with
##' \code{\link{vcov}}.} \item{data}{Data used in the model fit.  To extract
##' this in the standard R formats, use use
##' \code{\link{model.frame.flexsurvreg}} or
##' \code{\link{model.matrix.flexsurvreg}}.}
##' @section Custom distributions: \code{\link{flexsurvreg}} is intended to be
##' easy to extend to handle new distributions.  To define a new distribution
##' for use in \code{\link{flexsurvreg}}, construct a list with the following
##' elements:
##' 
##' \describe{ \item{list("name")}{A string naming the distribution.  If this
##' is called \code{"dist"}, for example, then there must be visible in the
##' working environment, at least, either
##' 
##' a) a function called \code{ddist} which defines the probability density,
##' 
##' or
##' 
##' b) a function called \code{hdist} which defines the hazard.
##' 
##' Ideally, in case a) there should also be a function called \code{pdist}
##' which defines the probability distribution or cumulative density, and in
##' case b) there should be a function called \code{Hdist} defining the
##' cumulative hazard.  If these additional functions are not provided,
##' \pkg{flexsurv} attempts to automatically create them by numerically
##' integrating the density or hazard function.  However, model fitting will be
##' much slower, or may not even work at all, if the analytic versions of these
##' functions are not available.
##' 
##' The functions must accept vector arguments (representing different times,
##' or alternative values for each parameter) and return the results as a
##' vector.  The function \code{\link{Vectorize}} may be helpful for doing
##' this: see the example below.

##' These functions may be in an add-on package (see below for an example) or
##' may be user-written.  If they are user-written they must be defined in the
##' global environment, or supplied explicitly through the \code{dfns} argument
##' to \code{flexsurvreg}.  The latter may be useful if the functions are
##' created dynamically (as in the source of \code{flexsurvspline}) and thus
##' not visible through R's scoping rules.
##' 
##' Arguments other than parameters must be named in the conventional way --
##' for example \code{x} for the first argument of the density function or
##' hazard, as in \code{\link{dnorm}(x, ...)} and \code{q} for the first
##' argument of the probability function.  Density functions should also have
##' an argument \code{log}, after the parameters, which when \code{TRUE},
##' computes the log density, using a numerically stable additive formula if
##' possible.
##' 
##' Additional functions with names beginning with \code{"DLd"} and
##' \code{"DLS"} may be defined to calculate the derivatives of the log density
##' and log survival probability, with respect to the parameters of the
##' distribution.  The parameters are expressed on the real line, for example
##' after log transformation if they are defined as positive.  The first
##' argument must be named \code{t}, representing the time, and the remaining
##' arguments must be named as the parameters of the density function. The
##' function must return a matrix with rows corresponding to times, and columns
##' corresponding to the parameters of the distribution.  The derivatives are
##' used, if available, to speed up the model fitting with \code{\link{optim}}.
##' }\item{:}{A string naming the distribution.  If this is called
##' \code{"dist"}, for example, then there must be visible in the working
##' environment, at least, either
##' 
##' a) a function called \code{ddist} which defines the probability density,
##' 
##' or
##' 
##' b) a function called \code{hdist} which defines the hazard.
##' 
##' Ideally, in case a) there should also be a function called \code{pdist}
##' which defines the probability distribution or cumulative density, and in
##' case b) there should be a function called \code{Hdist} defining the
##' cumulative hazard.  If these additional functions are not provided,
##' \pkg{flexsurv} attempts to automatically create them by numerically
##' integrating the density or hazard function.  However, model fitting will be
##' much slower, or may not even work at all, if the analytic versions of these
##' functions are not available.
##' 
##' The functions must accept vector arguments (representing different times,
##' or alternative values for each parameter) and return the results as a
##' vector.  The function \code{\link{Vectorize}} may be helpful for doing
##' this: see the example below.
##' 
##' These functions may be in an add-on package (see below for an example) or
##' may be user-written.  If they are user-written they must be defined in the
##' global environment, or supplied explicitly through the \code{dfns} argument
##' to \code{flexsurvreg}.  The latter may be useful if the functions are
##' created dynamically (as in the source of \code{flexsurvspline}) and thus
##' not visible through R's scoping rules.
##' 
##' Arguments other than parameters must be named in the conventional way --
##' for example \code{x} for the first argument of the density function or
##' hazard, as in \code{\link{dnorm}(x, ...)} and \code{q} for the first
##' argument of the probability function.  Density functions should also have
##' an argument \code{log}, after the parameters, which when \code{TRUE},
##' computes the log density, using a numerically stable additive formula if
##' possible.
##' 
##' Additional functions with names beginning with \code{"DLd"} and
##' \code{"DLS"} may be defined to calculate the derivatives of the log density
##' and log survival probability, with respect to the parameters of the
##' distribution.  The parameters are expressed on the real line, for example
##' after log transformation if they are defined as positive.  The first
##' argument must be named \code{t}, representing the time, and the remaining
##' arguments must be named as the parameters of the density function. The
##' function must return a matrix with rows corresponding to times, and columns
##' corresponding to the parameters of the distribution.  The derivatives are
##' used, if available, to speed up the model fitting with \code{\link{optim}}.
##' } \item{list("pars")}{Vector of strings naming the parameters of the
##' distribution. These must be the same names as the arguments of the density
##' and probability functions.  }\item{:}{Vector of strings naming the
##' parameters of the distribution. These must be the same names as the
##' arguments of the density and probability functions.  }
##' \item{list("location")}{Name of the main parameter governing the mean of
##' the distribution.  This is the default parameter on which covariates are
##' placed in the \code{formula} supplied to \code{flexsurvreg}. }\item{:}{Name
##' of the main parameter governing the mean of the distribution.  This is the
##' default parameter on which covariates are placed in the \code{formula}
##' supplied to \code{flexsurvreg}. } \item{list("transforms")}{List of R
##' functions which transform the range of values taken by each parameter onto
##' the real line.  For example, \code{c(log, log)} for a distribution with two
##' positive parameters. }\item{:}{List of R functions which transform the
##' range of values taken by each parameter onto the real line.  For example,
##' \code{c(log, log)} for a distribution with two positive parameters. }
##' \item{list("inv.transforms")}{List of R functions defining the
##' corresponding inverse transformations.  Note these must be lists, even for
##' single parameter distributions they should be supplied as, e.g.
##' \code{c(exp)} or \code{list(exp)}. }\item{:}{List of R functions defining
##' the corresponding inverse transformations.  Note these must be lists, even
##' for single parameter distributions they should be supplied as, e.g.
##' \code{c(exp)} or \code{list(exp)}. } \item{list("inits")}{A function of the
##' observed survival times \code{t} (including right-censoring times, and
##' using the halfway point for interval-censored times) which returns a vector
##' of reasonable initial values for maximum likelihood estimation of each
##' parameter.  For example, \code{function(t){ c(1, mean(t)) }} will always
##' initialize the first of two parameters at 1, and the second (a scale
##' parameter, for instance) at the mean of \code{t}.  }\item{:}{A function of
##' the observed survival times \code{t} (including right-censoring times, and
##' using the halfway point for interval-censored times) which returns a vector
##' of reasonable initial values for maximum likelihood estimation of each
##' parameter.  For example, \code{function(t){ c(1, mean(t)) }} will always
##' initialize the first of two parameters at 1, and the second (a scale
##' parameter, for instance) at the mean of \code{t}.  } }
##' 
##' For example, suppose we want to use an extreme value survival distribution.
##' This is available in the CRAN package \pkg{eha}, which provides
##' conventionally-defined density and probability functions called
##' \code{\link[eha]{dEV}} and \code{\link[eha]{pEV}}.  See the Examples below
##' for the custom list in this case, and the subsequent command to fit the
##' model.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{flexsurvspline}} for flexible survival modelling using
##' the spline model of Royston and Parmar.
##' 
##' \code{\link{plot.flexsurvreg}} and \code{\link{lines.flexsurvreg}} to plot
##' fitted survival, hazards and cumulative hazards from models fitted by
##' \code{\link{flexsurvreg}} and \code{\link{flexsurvspline}}.
##' @references Jackson, C. (2016). flexsurv: A Platform for Parametric
##' Survival Modeling in R. Journal of Statistical Software, 70(8), 1-33.
##' doi:10.18637/jss.v070.i08
##' 
##' Cox, C. (2008) The generalized \eqn{F} distribution: An umbrella for
##' parametric survival analysis.  Statistics in Medicine 27:4301-4312.
##' 
##' Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007) Parametric survival
##' analysis and taxonomy of hazard functions for the generalized gamma
##' distribution.  Statistics in Medicine 26:4252-4374
##' 
##' Jackson, C. H. and Sharples, L. D. and Thompson, S. G. (2010) Survival
##' models in health economic evaluations: balancing fit and parsimony to
##' improve prediction. International Journal of Biostatistics 6(1):Article 34.
##' @keywords models survival
##' @examples
##' 
##' data(ovarian)
##' ## Compare generalized gamma fit with Weibull
##' fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
##' fitg
##' fitw <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="weibull")
##' fitw
##' plot(fitg)
##' lines(fitw, col="blue", lwd.ci=1, lty.ci=1)
##' ## Identical AIC, probably not enough data in this simple example for a
##' ## very flexible model to be worthwhile.
##' 
##' ## Custom distribution
##' ## make "dEV" and "pEV" from eha package (if installed)
##' ##   available to the working environment
##' if (require("eha")) {
##' custom.ev <- list(name="EV",
##'                       pars=c("shape","scale"),
##'                       location="scale",
##'                       transforms=c(log, log),
##'                       inv.transforms=c(exp, exp),
##'                       inits=function(t){ c(1, median(t)) })
##' fitev <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian,
##'                     dist=custom.ev)
##' fitev
##' lines(fitev, col="purple", col.ci="purple")
##' }
##' 
##' 
##' ## Custom distribution: supply the hazard function only
##' hexp2 <- function(x, rate=1){ rate } # exponential distribution
##' hexp2 <- Vectorize(hexp2)
##' custom.exp2 <- list(name="exp2", pars=c("rate"), location="rate",
##'                     transforms=c(log), inv.transforms=c(exp),
##'                     inits=function(t)1/mean(t))
##' flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.exp2)
##' flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
##' ## should give same answer
##' 
##' @export
flexsurvreg <- function(formula, anc=NULL, data, weights, bhazard, subset, na.action, dist,
                        inits, fixedpars=NULL, dfns=NULL, aux=NULL, cl=0.95,
                        integ.opts=NULL, sr.control=survreg.control(), ...)
{
    call <- match.call()
    if (missing(dist)) stop("Distribution \"dist\" not specified")
    if (is.character(dist)) {
        # Added: case insensitve matching of distributions
        # Step 1: Use match.arg on lowercase argument, dists.
        # Step 2: Use match to get index in distribution list from
        # value obtained in step 1, and grab corresponding element.
        dist <- match.arg(tolower(dist), tolower(names(flexsurv.dists)))
        dist <- names(flexsurv.dists)[match(dist,tolower(names(flexsurv.dists)))]
        dlist <- flexsurv.dists[[dist]]
    }
    else if (is.list(dist)) {
        dlist <- check.dlist(dist)
    }
    else stop("\"dist\" should be a string for a built-in distribution, or a list for a custom distribution")
    dfns <- form.dp(dlist, dfns, integ.opts)
    parnames <- dlist$pars
    ancnames <- setdiff(parnames, dlist$location)

    check.formula(formula, dlist)
    if (is.null(anc)){
        anc <- vector(mode="list", length=length(ancnames))
        names(anc) <- ancnames
        for (i in ancnames){
            anc[[i]] <- ancpar.formula(formula, i)
        }
    }
    else {
        if (!is.list(anc) || !all(sapply(anc, function(x)inherits(x, "formula"))))
            stop("\"anc\" must be a list of formulae")
    }
    forms <- c(location=get.locform(formula, ancnames), anc)
    names(forms)[[1]] <- dlist$location

    ## a) calling model.frame() directly doesn't work.  it only looks in
    ## "data" or the environment of "formula" for extra variables like
    ## "weights". needs to also look in environment of flexsurvreg.
    ## m <- model.frame(formula=, data=data, weights=weights, subset=subset, na.action=na.action)
    ## b) putting block below in a function doesn't work when calling
    ## flexsurvreg within a function
    ## m <- make.model.frame(call, formula, data, weights, subset, na.action, ancnames)

    ## Make model frame
    indx <- match(c("formula", "data", "weights", "bhazard", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")

    f2 <- concat.formulae(formula,forms)
    temp[["formula"]] <- f2
    if (missing(data)) temp[["data"]] <- environment(formula)
    m <- eval(temp, parent.frame())
    m <- droplevels(m) # remove unused factor levels after subset applied
    attr(m,"covnames") <- attr(f2, "covnames") # for "newdata" in summary
    attr(m,"covnames.orig") <- intersect(colnames(m), attr(f2, "covnames.orig")) # for finding factors in plot method
    Y <- check.flexsurv.response(model.extract(m, "response"))
    mml <- mx <- vector(mode="list", length=length(dlist$pars))
    names(mml) <- names(mx) <- c(dlist$location, setdiff(dlist$pars, dlist$location))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], m)
        mx[[i]] <- length(unlist(mx)) + seq_len(ncol(mml[[i]][,-1,drop=FALSE]))
    }
    X <- compress.model.matrices(mml)

    weights <- model.extract(m, "weights")
    if (is.null(weights)) weights <- m$"(weights)" <- rep(1, nrow(X))
    bhazard <- model.extract(m, "bhazard")
    if (is.null(bhazard)) bhazard <- rep(0, nrow(X))
    dat <- list(Y=Y, m=m, mml=mml)
    ncovs <- length(attr(m, "covnames.orig"))

    ncoveffs <- ncol(X)
    nbpars <- length(parnames) # number of baseline parameters
    npars <- nbpars + ncoveffs

    if (missing(inits) && is.null(dlist$inits))
        stop("\"inits\" not supplied, and no function to estimate them found in the custom distribution list")
    if (missing(inits) || any(is.na(inits))){
        yy <- ifelse(Y[,"status"]==3 & is.finite(Y[,"time2"]), (Y[,"time1"] + Y[,"time2"])/2, Y[,"time1"])
        wt <- yy*weights*length(yy)/sum(weights)
        dlist$inits <- expand.inits.args(dlist$inits)
        inits.aux <- c(aux, list(forms=forms, data=if(missing(data)) NULL else data, weights=temp$weights,
                                 control=sr.control,
                                 counting=(attr(model.extract(m, "response"), "type")=="counting")
                                 ))
        auto.inits <- dlist$inits(t=wt,mf=m,mml=mml,aux=inits.aux)
        if (!missing(inits) && any(is.na(inits))) inits[is.na(inits)] <- auto.inits[is.na(inits)]
        else inits <- auto.inits
    }
    if (!is.numeric(inits)) stop ("initial values must be a numeric vector")
    nin <- length(inits)
    if (nin < npars && ncoveffs > 0)
        inits <- c(inits, rep(0,length=npars-nin))
    else if (nin > npars){
        inits <- inits[1:npars]
        warning("Initial values are a vector length > ", npars, ": using only the first ", npars)
    }
    else if (nin < nbpars){
        stop("Initial values are a vector length ", nin, ", but distribution has ",nbpars, " parameters")
    }

    for (i in 1:nbpars)
        inits[i] <- dlist$transforms[[i]](inits[i])
    outofrange <- which(is.nan(inits) | is.infinite(inits))
    if (any(outofrange)){
        plural <- if(length(outofrange) > 1) "s" else ""
        stop("Initial value", plural, " for parameter", plural, " ",
             paste(outofrange,collapse=","), " out of range")
    }
    names(inits) <- c(parnames, colnames(X))

    if (!is.null(fixedpars) && !is.logical(fixedpars) &&
        (!is.numeric(fixedpars) || any(!(fixedpars %in% 1:npars)))){
        dots <- if(npars>2) "...," else ""
        stop("fixedpars must be TRUE/FALSE or a vector of indices in 1,",dots,npars)
    }

    if ((is.logical(fixedpars) && fixedpars==TRUE) ||
        (is.numeric(fixedpars) && identical(fixedpars, 1:npars))) {
        minusloglik <- minusloglik.flexsurv(inits, Y=Y, X=X,
                                            weights=weights, bhazard=bhazard,
                                            dlist=dlist, inits=inits,
                                            dfns=dfns, aux=aux, mx=mx)
        res.t <- matrix(inits, ncol=1)
        inits.nat <- inits
        for (i in 1:nbpars)
            inits.nat[i] <- dlist$inv.transforms[[i]](inits[i])
        res <- matrix(inits.nat, ncol=1)
        dimnames(res) <- dimnames(res.t) <- list(names(inits), "est")
        ret <- list(res=res, res.t=res.t, npars=0,
                    loglik=-as.vector(minusloglik), logliki=attr(minusloglik,"indiv"))
    }
    else {
        optpars <- inits[setdiff(1:npars, fixedpars)]
        optim.args <- list(...)
        if (is.null(optim.args$method)){
            optim.args$method <- "BFGS"
        }
        gr <- if (dfns$deriv) Dminusloglik.flexsurv else NULL
        optim.args <- c(optim.args,
                        list(par=optpars,
                             fn=logLikFactory(Y=Y, X=X,
                                              weights=weights,
                                              bhazard=bhazard,
                                              inits=inits, dlist=dlist,
                                              dfns=dfns,
                                              aux=aux, mx=mx,
                                              fixedpars=fixedpars),
                             gr=gr,
                             Y=Y, X=X, weights=weights,
                             bhazard=bhazard, dlist=dlist,
                             inits=inits, dfns=dfns, aux=aux,
                             mx=mx, fixedpars=fixedpars,
                             hessian=TRUE))
        opt <- do.call("optim", optim.args)
        est <- opt$par
        if (all(!is.na(opt$hessian)) && all(!is.nan(opt$hessian)) && all(is.finite(opt$hessian)) &&
            all(eigen(opt$hessian)$values > 0))
        {
            cov <- solve(opt$hessian); se <- sqrt(diag(cov))
            if (!is.numeric(cl) || length(cl)>1 || !(cl>0) || !(cl<1))
                stop("cl must be a number in (0,1)")
            lcl <- est - qnorm(1 - (1-cl)/2)*se
            ucl <- est + qnorm(1 - (1-cl)/2)*se
        }
        else {
            warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite. ")
            cov <- lcl <- ucl <- se <- NA
        }
        res <- cbind(est=inits, lcl=NA, ucl=NA, se=NA)
        res[setdiff(1:npars, fixedpars),] <- cbind(est, lcl, ucl, se)
        colnames(res) <- c("est", paste(c("L","U"), round(cl*100), "%", sep=""), "se")
        res.t <- res # results on transformed (log) scale
        for (i in 1:nbpars){ # results on natural scale
            res[i,1:3] <- dlist$inv.transforms[[i]](res[i,1:3])
            if (identical(body(dlist$transforms[[i]]), body(log)))
                res[i,"se"] <- exp(res.t[i,"est"])*res.t[i,"se"]
            else if (identical(body(dlist$transforms[[i]]), body(logh)))
                res[i,"se"] <- dexph(res.t[i,"est"])*res.t[i,"se"]
            else if (!identical(dlist$transforms[[i]], identity))
                res[i,"se"] <- NA
            ## theoretically could also do logit SE(g(x) = exp(x)/(1 + exp(x))) = g'(x) SE(x);  g'(x) = exp(x)/(1 + exp(x))^2
            ## or any interval scale (dglogit) as in msm
        }
        minusloglik <- minusloglik.flexsurv(res.t[,"est"], Y=Y, X=X, weights=weights, bhazard=bhazard,
                                            dlist=dlist, inits=inits, dfns=dfns, aux=aux, mx=mx)
        ret <- list(res=res, res.t=res.t, cov=cov, coefficients=res.t[,"est"],
                    npars=length(est), fixedpars=fixedpars, optpars=setdiff(1:npars, fixedpars),
                    loglik=-opt$value, logliki=attr(minusloglik,"indiv"),
                    cl=cl, opt=opt)
    }
    ret <- c(list(call=call, dlist=dlist, aux=aux,
                  ncovs=ncovs, ncoveffs=ncoveffs,
                  mx=mx, basepars=1:nbpars,
                  covpars=if (ncoveffs>0) (nbpars+1):npars else NULL,
                  AIC=-2*ret$loglik + 2*ret$npars,
                  data = dat, datameans = colMeans(X),
                  N=nrow(dat$Y), events=sum(dat$Y[,"status"]==1), trisk=sum(dat$Y[,"time"]),
                  concat.formula=f2, all.formulae=forms, dfns=dfns),
             ret)
    if (isTRUE(getOption("flexsurv.test.analytic.derivatives"))
        && (dfns$deriv) ) {
        if (is.logical(fixedpars) && fixedpars==TRUE) { optpars <- inits; fixedpars=FALSE }
        ret$deriv.test <- deriv.test(optpars, Y, X, weights, bhazard, dlist, inits, dfns, aux, mx, fixedpars)
    }
    class(ret) <- "flexsurvreg"
    ret
}

##' @export
print.flexsurvreg <- function(x, ...)
{
    covmeans <- colMeans(model.matrix(x))
    covs <- names(covmeans)
    covinds <- match(covs, rownames(x$res))
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    if (x$npars > 0) {
        res <- x$res
        cat ("Estimates: \n")
        if (any(covinds)) {
            ecoefs <- matrix(NA, nrow=nrow(x$res), ncol=3)
            colnames(ecoefs) <- c("exp(est)", colnames(res)[2:3])
            means <- rep(NA,nrow(x$res))
            ecoefs[covinds,] <- exp(x$res[covinds,1:3,drop=FALSE])
            means[covinds] <- covmeans
            res <- cbind(means, res, ecoefs)
            colnames(res)[1] <- "data mean"
        }
        args <- list(...)
        if (is.null(args$digits)) args$digits <- 3
        f <- do.call("format", c(list(x=res), args))
        print(f, print.gap=2, quote=FALSE, na.print="")
    }
    cat("\nN = ", x$N, ",  Events: ", x$events,
        ",  Censored: ", x$N - x$events,
        "\nTotal time at risk: ", x$trisk,
        "\nLog-likelihood = ", x$loglik, ", df = ", x$npars,
        "\nAIC = ", x$AIC, "\n\n", sep="")
}

form.model.matrix <- function(object, newdata){
    mfo <- model.frame(object)

    ## If required covariate missing, give a slightly more informative error message than, e.g.
    ## "Error in eval(expr, envir, enclos) (from flexsurvreg.R#649) : object 'sex' not found"
    covnames <- attr(mfo, "covnames")
    missing.covs <- unique(covnames[!covnames %in% names(newdata)])
    if (length(missing.covs) > 0){
        missing.covs <- sprintf("\"%s\"", missing.covs)
        plural <- if (length(missing.covs)>1) "s" else ""
        stop(sprintf("Value%s of covariate%s ",plural,plural), paste(missing.covs, collapse=", "), " not supplied in \"newdata\"")
    }

    ## as in predict.lm
    tt <- attr(mfo, "terms")
    Terms <- delete.response(tt)
    mf <- model.frame(Terms, newdata, xlev = .getXlevels(tt, mfo))
    if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, mf)

    forms <- object$all.formulae
    mml <- vector(mode="list", length=length(object$dlist$pars))
    names(mml) <- names(forms)
    forms[[1]] <- delete.response(terms(forms[[1]]))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], mf)
    }
    X <- compress.model.matrices(mml)

    attr(X, "newdata") <- mf # newdata with any extra variables stripped.  Used to name components of summary list
    X
}


##' @export
vcov.flexsurvreg <- function (object, ...)
{
    object$cov
}

.onLoad <- function(libname, pkgname) {
    assign("flexsurv.dfns", new.env(), envir=parent.env(environment()))
}



##' Extract original data from \code{flexsurvreg} objects.
##' 
##' Extract the data from a model fitted with \code{flexsurvreg}.
##' 
##' 
##' @aliases model.frame.flexsurvreg model.matrix.flexsurvreg
##' @param formula A fitted model object, as returned by
##' \code{\link{flexsurvreg}}.
##' @param object A fitted model object, as returned by
##' \code{\link{flexsurvreg}}.
##' @param par String naming the parameter whose linear model matrix is
##' desired.
##' 
##' The default value of \code{par=NULL} returns a matrix consisting of the
##' model matrices for all models in the object \code{cbind}ed together, with
##' the intercepts excluded.  This is not really a ``model matrix'' in the
##' usual sense, however, the columns directly correspond to the covariate
##' coefficients in the matrix of estimates from the fitted model.
##' @param ... Further arguments (not used).
##' @return \code{model.frame} returns a data frame with all the original
##' variables used for the model fit.
##' 
##' \code{model.matrix} returns a design matrix for a part of the model that
##' includes covariates.  The required part is indicated by the \code{"par"}
##' argument (see above).
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{flexsurvreg}}, \code{\link{model.frame}},
##' \code{\link{model.matrix}}.
##' @keywords models
##' @export
model.frame.flexsurvreg <- function(formula, ...)
{
    x <- formula
    x$data$m
}

##' @export
##' @rdname model.frame.flexsurvreg
model.matrix.flexsurvreg <- function(object, par=NULL, ...)
{
    x <- object
    if (is.null(par)) compress.model.matrices(x$data$mml) else x$data$mml[[par]]
}

##' @export
logLik.flexsurvreg <- function(object, ...){
    val <- object$loglik
    attr(val, "df") <- object$npars
    attr(val, "nobs") <- nrow(model.frame(object))
    class(val) <- "logLik"
    val
}

##' Extract model coefficients from fitted flexible survival models
##' 
##' Extract model coefficients from fitted flexible survival models.  This
##' presents all parameter estimates, transformed to the real line if necessary.
##' For example, shape or scale parameters, which are constrained to be
##' positive, are returned on the log scale.
##' 
##' This matches the behaviour of \code{coef.default} for standard R model
##' families such as \code{\link[stats]{glm}}, where intercepts in regression
##' models are presented on the same scale as the covariate effects.  Note that
##' any parameter in a distribution fitted by \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvreg}} may be an intercept in a regression model.
#' 
##' @param object Output from \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvspline}}, representing a fitted survival model object.
##' @param ... Further arguments passed to or from other methods.  Currently
##' unused.
##' @return This returns the \code{mod$res.t[,"est"]} component of the fitted
##' model object \code{mod}.  See \code{\link{flexsurvreg}},
##' \code{\link{flexsurvspline}} for full documentation of all components.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \code{\link{flexsurvreg}}, \code{\link{flexsurvspline}}.
##' @keywords models
##' @export
coef.flexsurvreg <- function(object, ...){
    object$coefficients
}
