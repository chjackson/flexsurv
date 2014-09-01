## exp(asinh(y))
exph <- function(y)(y + sqrt(y^2 + 1))
dexph <- function(y)(1 + y / sqrt(y^2 + 1))

## sinh(log(y))
logh <- function(x)((x^2 - 1)/(2*x))

flexsurv.dists <- list(
                       genf = list(
                       name="genf",
                       pars=c("mu","sigma","Q","P"),
                       location="mu",
                       transforms=c(identity, log, identity, log),
                       inv.transforms=c(identity, exp, identity, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 0, 1)
                       }
                       ),
                       genf.orig = list(
                       name="genf.orig",
                       pars=c("mu","sigma","s1","s2"),
                       location="mu",
                       transforms=c(identity, log, log, log),
                       inv.transforms=c(identity, exp, exp, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 1, 1)
                       }
                       ),
                       gengamma = list(
                       name="gengamma",
                       pars=c("mu","sigma","Q"),
                       location="mu",
                       transforms=c(identity, log, identity),
                       inv.transforms=c(identity, exp, identity),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt), 0)
                       }
                       ),
                       gengamma.orig = list(
                       name="gengamma.orig",
                       pars=c("shape","scale","k"),
                       location="scale",
                       transforms=c(log, log, log),
                       inv.transforms=c(exp, exp, exp),
                       inits=function(t){c(1, mean(t), 1)}
                       ),
                       exp = list(
                       name="exp",
                       pars=c("rate"),
                       location="rate",
                       transforms=c(log),
                       inv.transforms=c(exp),
                       inits=function(t){1 / mean(t)}
                       ),
                       weibull = list(
                       name = "weibull.quiet",
                       pars=c("shape","scale"),
                       location="scale",
                       transforms=c(log, log),
                       inv.transforms=c(exp, exp),
                       inits = function(t){
## TODO regress log CH on log time as in splineinits. gets cov effs as well
                           lt <- log(t[t>0])
                           c(1, exp(median(lt)) / log(2))
                       }
                       ),
                       lnorm = list(
                       name="lnorm",
                       pars=c("meanlog","sdlog"),
                       location="meanlog",
                       transforms=c(identity, log),
                       inv.transforms=c(identity, exp),
                       inits=function(t){
                           lt <- log(t[t>0])
                           c(mean(lt), sd(lt))
                       }
                       ),
                       gamma = list(
                       name="gamma",
                       pars=c("shape","rate"),
                       location="rate",
                       transforms=c(log, log),
                       inv.transforms=c(exp, exp),
                       inits=function(t){
                           m=mean(t); v=var(t);
                           c(m^2/v, m/v)
                       }
                       ),
                       gompertz = list(
                       name="gompertz",
                       pars=c("shape","rate"),
                       location="rate",
                       transforms=c(identity, log),
                       inv.transforms=c(identity, exp),
                       inits=function(t){c(0.001,1 / mean(t))}
                       )
                       )

minusloglik.flexsurv <- function(optpars, Y, X=0, weights, bhazard, dlist, inits,
                                 dfns, aux, mx, fixedpars=NULL) {
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    nbpars <- length(dlist$pars)
    pars <- as.list(pars)
    if (npars > nbpars) {
        beta <- unlist(pars[(nbpars+1):npars])
        for (i in dlist$pars)
            pars[[i]] <- pars[[i]] + X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
    }
    pcall <- list()
    for (i in 1:nbpars){
        pcall[[names(pars)[i]]] <- dlist$inv.transforms[[i]](pars[[i]])
    }
    for (i in seq_along(aux)){
        if (!(names(aux)[i] %in% dlist$pars))
            pcall[[names(aux)[i]]] <- aux[[i]]
    }
    pmaxcall <- dcall <- tcall <- pcall
    pcall$q <- Y[,"time1"]
    pmaxcall$q <- Y[,"time2"] # Inf if right-censored, giving pmax=1
    dcall$x <- Y[,"time1"]
    tcall$q <- Y[,"start"]  
    dcall$log <- TRUE
    ## Generic survival model likelihood
    dead <- Y[,"status"]==1
    logdens <- (do.call(dfns$d, dcall))
    pmax <- (do.call(dfns$p, pmaxcall))
    pmin <- (do.call(dfns$p, pcall))
    pobs <- 1 - do.call(dfns$p, tcall) # prob of being observed = 1 unless left-truncated
    
    ## Hazard offset for relative survival models
    if (any(bhazard>0)){
        loghaz <- logdens - log(pmax - pmin)
        offset <- sum(log(1 + bhazard / exp(loghaz)*weights)[dead])
    } else offset <- 0
    
    ret <- - ( sum((logdens*weights)[dead]) +
              sum((log(pmax - pmin)*weights)[!dead]) -
              sum(log(pobs)*weights) + offset)
    ret
}

check.dlist <- function(dlist){
    ## put tests in testthat
    if (is.null(dlist$name)) stop("\"name\" element of custom distribution list not given")
    if (!is.character(dlist$name)) stop("\"name\" element of custom distribution list should be a string")
    if (is.null(dlist$pars)) stop("parameter names \"pars\" not given in custom distribution list")
    if (!is.character(dlist$pars)) stop("parameter names \"pars\" should be a character vector")
    npars <- length(dlist$pars)
    if (is.null(dlist$location)) stop("location parameter not given in custom distribution list")
    if (!(dlist$location %in% dlist$pars)) stop("location parameter \"",dlist$location,"\" not in list of parameters")
    if (is.null(dlist$transforms)) stop("transforms not given in custom distribution list")
    if (is.null(dlist$inv.transforms)) stop("inverse transforms not given in custom distribution list")
    if (!is.list(dlist$transforms)) stop("\"transforms\" must be a list of functions")
    if (!is.list(dlist$inv.transforms)) stop("\"inv.transforms\" must be a list of functions")
    if (length(dlist$transforms) != npars) stop("transforms vector of length ",length(dlist$transforms),", parameter names of length ",npars)
    if (length(dlist$inv.transforms) != npars) stop("inverse transforms vector of length ",length(dlist$transforms),", parameter names of length ",npars) #
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

check.flexsurv.response <- function(Y){
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
### convert Y from Surv object to numeric matrix
### though "time" only used for initial values, printed time at risk, empirical hazard
    if (attr(Y, "type") == "counting")
        Y <- cbind(Y, time=Y[,"stop"] - Y[,"start"], time1=Y[,"stop"], time2=Inf)
    else if (attr(Y, "type") == "interval")
        Y <- cbind(Y, start=0, stop=Y[,"time1"], time=Y[,"time1"])  
    else if (attr(Y, "type") == "right")
        Y <- cbind(Y, start=0, stop=Y[,"time"], time1=Y[,"time"], time2=Inf)
    else stop("Survival object type \"", attr(Y, "type"), "\"", " not supported")
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

flexsurvreg <- function(formula, anc=NULL, data, weights, bhazard, subset, na.action, dist,
                        inits, fixedpars=NULL, dfns=NULL, aux=NULL, cl=0.95,
                        integ.opts=NULL, ...)
{
    call <- match.call()
    if (missing(dist)) stop("Distribution \"dist\" not specified")
    if (is.character(dist)) {
        match.arg(dist, names(flexsurv.dists))
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
    attr(m,"covnames") <- attr(f2, "covnames") # for "newdata" in summary
    Y <- check.flexsurv.response(model.extract(m, "response"))
    mml <- mx <- vector(mode="list", length=length(dlist$pars))
    names(mml) <- names(mx) <- c(dlist$location, setdiff(dlist$pars, dlist$location))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], m)
        mx[[i]] <- length(unlist(mx)) + seq_len(ncol(mml[[i]][,-1,drop=FALSE]))
    }
    X <- compress.model.matrices(mml)
    
    weights <- model.extract(m, "weights")
    if (is.null(weights)) weights <- rep(1, nrow(X))
    bhazard <- model.extract(m, "bhazard")
    if (is.null(bhazard)) bhazard <- rep(0, nrow(X))
    dat <- list(Y=Y, m=m, mml=mml)
    ncovs <- ncol(m) - 1
    ncoveffs <- ncol(X)
    nbpars <- length(parnames) # number of baseline parameters
    npars <- nbpars + ncoveffs
    if (!missing(inits) && (!is.numeric(inits) || (length(inits) != npars)))
        stop("inits must be a numeric vector of length ",npars)
    if (missing(inits) && is.null(dlist$inits))
        stop("\"inits\" not supplied, and no function to estimate them found in the custom distribution list")
    if (missing(inits) || any(is.na(inits))){
        wt <- Y[,"time"]*weights*length(Y[,"time"])/sum(weights)
        dlist$inits <- expand.inits.args(dlist$inits)
        base.inits <- dlist$inits(t=wt,mf=m,mml=mml,aux=aux)
        if (!is.numeric(base.inits)) stop ("\"inits\" function must return a numeric vector")
        if (!(length(base.inits) %in% c(nbpars, npars))){
            if (npars == nbpars)
                stop("\"inits\" function must return a numeric vector of length ", nbpars, " = number of parameters")
            else
                stop("\"inits\" function must return a numeric vector of length ", nbpars, " or ", npars, ", number of parameters including or excluding covariate effects")
        }
        default.inits <- base.inits
        if (length(base.inits) < npars)
            default.inits <- c(default.inits, rep(0,ncoveffs))
    }
    if (missing(inits)) inits <- default.inits
    else if (any(is.na(inits))) inits[is.na(inits)] <- default.inits[is.na(inits)]
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
        minusloglik <- minusloglik.flexsurv(inits, Y=Y, X=X, weights=weights, bhazard=bhazard,
                                            dlist=dlist, inits=inits, dfns=dfns, aux=aux, mx=mx)
        inits.t <- numeric(length(inits))
        for (i in 1:nbpars)
            inits.t[i] <- dlist$inv.transforms[[i]](inits[i])
        res <- matrix(inits.t, ncol=1)
        dimnames(res) <- list(names(inits.t), "est")
        ret <- list(res=res, npars=0, loglik=-minusloglik)
    }
    else {
        optpars <- inits[setdiff(1:npars, fixedpars)]
        optim.args <- list(...)
        if (is.null(optim.args$method)){
            ## L-BFGS-B breaks if lik goes to 0 for finite pars, but BFGS carries on (with warning) if pars zero
            ## if (dlist$name=="weibull"){            
            ##     optim.args$method <- "L-BFGS-B"
            ##     optim.args$lower <- -1e+16
            ## }
            ## else
                optim.args$method <- "BFGS"
        }
        gr <- if (dfns$deriv) Dminusloglik.flexsurv else NULL
        optim.args <- c(optim.args, list(par=optpars, fn=minusloglik.flexsurv, gr=gr,
                                         Y=Y, X=X, weights=weights, bhazard=bhazard, dlist=dlist,
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
        ret <- list(res=res, res.t=res.t, cov=cov, coefficients=res.t[,"est"],
                    npars=length(est), fixedpars=fixedpars, optpars=setdiff(1:npars, fixedpars),
                    mx=mx, ncovs=ncovs, ncoveffs=ncoveffs, basepars=1:nbpars, 
                    covpars=if (ncoveffs>0) (nbpars+1):npars else NULL,
                    loglik=-opt$value, cl=cl, opt=opt)
    }
    ret <- c(list(call=call, dlist=dlist, aux=aux,
                  AIC=-2*ret$loglik + 2*ret$npars,
                  data = dat, datameans = colMeans(X),
                  N=nrow(dat$Y), events=sum(dat$Y[,"status"]), trisk=sum(dat$Y[,"time"]),
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
    
print.flexsurvreg <- function(x, ...)
{
    covmeans <- colMeans(model.matrix(x))
    covs <- names(covmeans)
    covinds <- match(covs, rownames(x$res))
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
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
        print(format(res, digits=3), print.gap=2, quote=FALSE, na.print="")
    }
    cat("\nN = ", x$N, ",  Events: ", x$events,
        ",  Censored: ", x$N - x$events,
        "\nTotal time at risk: ", x$trisk,
        "\nLog-likelihood = ", x$loglik, ", df = ", x$npars,
        "\nAIC = ", x$AIC, "\n\n", sep="")
}

form.model.matrix <- function(object, newdata){
    mfo <- model.frame(object)
    covnames <- attr(mfo, "covnames")
    missing.covs <- unique(covnames[!covnames %in% names(newdata)])
    if (length(missing.covs)>0)
        stop("Values of covariates ", paste(missing.covs, collapse=", "), " not supplied in \"newdata\"")
    ## don't insist on user defining factors in model as factors in newdata, do this for them
    facs <- sapply(mfo, is.factor)
    facnames <- gsub(".+\\((.+)\\)","\\1",names(facs))
    for (i in which(facs))
        newdata[,facnames[i]] <- factor(newdata[,facnames[i]])
    temp <- call("model.frame")
    f2 <- delete.response(terms(object$concat.formula))
    temp[["formula"]] <- f2
    temp[["data"]] <- newdata
    mf <- eval(temp, parent.frame())
    facs <- names(mf)[sapply(mf, is.factor)]
    for (i in facs) mf[,i] <- factor(mf[,i], levels=levels(mfo[,i]))
    forms <- object$all.formulae
    mml <- vector(mode="list", length=length(object$dlist$pars))
    names(mml) <- object$dlist$pars
    forms[[1]] <- delete.response(terms(forms[[1]]))
    for (i in names(forms)){
        mml[[i]] <- model.matrix(forms[[i]], mf)
    }
    X <- compress.model.matrices(mml)
    X
}

summary.flexsurvreg <- function(object, newdata=NULL, X=NULL, type="survival", fn=NULL, 
                                t=NULL, start=0, ci=TRUE, B=1000, cl=0.95,
                                ...)
{
    x <- object
    dat <- x$data
    Xraw <- model.frame(x)[,-1,drop=FALSE]
    isfac <- sapply(Xraw,is.factor)
    ncovs <- x$ncovs
    type <- match.arg(type, c("survival","cumhaz","hazard"))
    if (is.null(newdata)){
        if (is.vector(X)) X <- matrix(X, nrow=1)
        if (ncovs > 0 && is.null(X)) {
            ## if any continuous covariates, calculate fitted survival for "average" covariate value
            if (!all(isfac))
                X <- matrix(colMeans(model.matrix(x)) ,nrow=1)
            ## else calculate for all different factor groupings
            else {
                X <- unique(model.matrix(x))
                ## build names like "COVA=value1,COVB=value2"
                nam <- as.matrix(unique(Xraw))
                for (i in 1:ncol(nam)) nam[,i] <- paste(colnames(nam)[i], nam[,i], sep="=")
                rownames(X) <- apply(nam, 1, paste, collapse=",")
            }
        }
        else if (is.null(X)) X <- as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1))
        else if (!is.matrix(X) || (is.matrix(X) && ncol(X) != x$ncoveffs)) {
            plural <- if (x$ncoveffs > 1) "s" else ""
            stop("expected X to be a matrix with ", x$ncoveffs, " column", plural, " or a vector with ", x$ncoveffs, " elements")
        }
    } else
        X <- form.model.matrix(object, as.data.frame(newdata))
    if (is.null(t))
        t <- sort(unique(dat$Y[,"stop"]))
    if (length(start)==1)
        start <- rep(start, length(t))
    else if (length(start) != length(t))
        stop("length of \"start\" is ",length(start)," should be 1, or length of \"t\" which is ",length(t))

    if (is.null(fn)) {
        fn <- summary.fns(x, type)
    }
    fn <- expand.summfn.args(fn)
    fncall <- list(t,start)
    beta <- if (ncovs==0) 0 else x$res[x$covpars,"est"]
    if (ncol(X) != length(beta)){
        isare <- if(length(beta)==1) "is" else "are"
        plural <- if(ncol(X)==1) "" else "s"
        pluralc <- if(length(beta)==1) "" else "s"
        stop("Supplied X has ", ncol(X), " column",plural," but there ",isare," ",
             length(beta), " covariate effect", pluralc)
    }
    dlist <- x$dlist
    ret <- vector(nrow(X), mode="list")
    if(!is.null(newdata))
        covnames <- apply(as.data.frame(newdata), 1, function(x)paste0(names(newdata), "=", x, collapse=", "))
    else covnames <- rownames(X)
    names(ret) <- covnames
    for (i in 1:nrow(X)) {
        basepars.mat <- add.covs(x, x$res.t[dlist$pars,"est"], beta, X[i,,drop=FALSE], transform=FALSE)
        basepars <- as.list(as.data.frame(basepars.mat))
        fncall[dlist$pars] <- basepars

#        for (j in seq(along=dlist$pars)) {
#            fncall2[[dlist$pars[j]]] <- x$res[dlist$pars[j],"est"]
#            mu <- dlist$transforms[[j]](fncall2[[dlist$pars[j]]]) ## TODO add.covs function?
#            covinds <- x$mx[[dlist$pars[j]]]
#            mu <- mu + X[i,covinds,drop=FALSE] %*% beta[covinds]
#            fncall2[[dlist$pars[j]]] <- dlist$inv.transforms[[j]](mu)
#        }
#        expect_equal(unlist(fncall2), unlist(fncall))
        
        for (j in seq_along(x$aux)){
            fncall[[names(x$aux)[j]]] <- x$aux[[j]]
        }
# browser()
        y <- do.call(fn, fncall)
        if (ci){
            res.ci <- cisumm.flexsurvreg(x, t, start, X[i,,drop=FALSE], fn=fn, B=B, cl=cl)
            ly <- res.ci[,1]
            uy <-  res.ci[,2]
        }
        ret[[i]] <- data.frame(time=t, est=y, row.names=NULL)
        if (ci) { ret[[i]]$lcl <- ly; ret[[i]]$ucl <- uy}
    }
    if (ncovs>0) attr(ret,"X") <- X
    class(ret) <- "summary.flexsurvreg"
    ret
}

summary.fns <- function(x, type){
    switch(type,   # TODO warn for clashing arguments in dfns
           "survival" = function(t,start,...) {
               ret <- (1 - x$dfns$p(t,...))/(1 - x$dfns$p(start,...))
               ret[t<start] <- 1 # prob[t<start] was previously 0
               ret
           },
           "hazard" = function(t,start,...) {
               ret <- x$dfns$h(t,...) * (1 - x$dfns$p(start,...))
               ret[t<start] <- 0
               ret
           },
           "cumhaz" = function(t,start,...) {
               ret <- x$dfns$H(t,...) - x$dfns$H(start,...)
               ret[t<start] <- 0
               ret                         
           })
}

print.summary.flexsurvreg <- function(x, ...){
    for (i in seq_along(x)){
        cat(names(x)[i], "\n")
        print(x[[i]])
        if (i<length(x)) cat("\n")
    }
}

add.covs <- function(x, pars, beta, X, transform=FALSE){  ## TODO option to transform on input 
    nres <- nrow(X)
    if (!is.matrix(pars)) pars <- matrix(pars, nrow=nres, ncol=length(pars), byrow=TRUE)
    if (!is.matrix(beta)) beta <- matrix(beta, nrow=1)
    for (j in seq(along=x$dlist$pars)){
        covinds <- x$mx[[x$dlist$pars[j]]]
        if (length(covinds) > 0){
            pars[,j] <- pars[,j] + beta[,covinds] %*% t(X[,covinds,drop=FALSE])
        }
        if (!transform)            
            pars[,j] <- x$dlist$inv.transforms[[j]](pars[,j])
    }
    colnames(pars) <- x$dlist$pars
    pars
}

## Draw B samples from multivariate normal distribution of baseline
## parameter estimators, for given covariate values

normboot.flexsurvreg <- function(x, B, newdata=NULL, X=NULL, transform=FALSE){
    if (is.null(X)) {
        if (is.null(newdata)) stop("neither \"newdata\" nor \"X\" supplied")
        X <- form.model.matrix(x, as.data.frame(newdata))
    }
    sim <- matrix(nrow=B, ncol=nrow(x$res))
    colnames(sim) <- rownames(x$res)
    sim[,x$optpars] <- rmvnorm(B, x$opt$par, x$cov)
    sim[,x$fixedpars] <- rep(x$res.t[x$fixedpars,"est"],each=B)
    beta <- sim[, x$covpars, drop=FALSE]
    if (nrow(X)==1){
        res <- sim[,x$dlist$pars,drop=FALSE]
        res <- add.covs(x=x, pars=res, beta=beta, X=X, transform=transform)
    }  else {
        res <- vector(nrow(X), mode="list")
        for (i in 1:nrow(X)) {
            res[[i]] <- sim[,x$dlist$pars,drop=FALSE]
            res[[i]] <- add.covs(x=x, pars=res[[i]], beta=beta, X=X[i,,drop=FALSE], transform=transform)
        }
    }
    res
}

### Compute CIs for survival, cumulative hazard or hazard at supplied
### times t and covariates X, using random sample of size B from the
### assumed MVN distribution of MLEs.

normbootfn.flexsurvreg <- function(x, t, start, newdata=NULL, X=NULL, fn, B){
    sim <- normboot.flexsurvreg(x, B, X=X)
    if (!is.list(sim)) sim <- list(sim)
    ret <- array(NA_real_, dim=c(nrow(X), B, length(t)))
    for (k in 1:nrow(X)){
        for (i in seq(length=B)) {
            fncall <- list(t,start)
            for (j in seq(along=x$dlist$pars))
                fncall[[x$dlist$pars[j]]] <- sim[[k]][i,j]
            for (j in seq_along(x$aux))
                fncall[[names(x$aux)[j]]] <- x$aux[[j]]
            ret[k,i,] <- do.call(fn, fncall)
        }
    }
    if (nrow(X)==1) ret[1,,,drop=FALSE] else ret
}

cisumm.flexsurvreg <- function(x, t, start, X, fn, B=1000, cl=0.95) {
    if (any(is.na(x$res[,2])) || (B==0))
        ret <- array(NA, dim=c(length(t), 2))
    else {
        ret <- normbootfn.flexsurvreg(x=x, t=t, start=start, X=X, fn=fn, B=B)
        ret <- apply(ret, c(1,3), function(x)quantile(x, c((1-cl)/2, 1 - (1-cl)/2), na.rm=TRUE))
        ret <- t(ret[,1,])
    }
    ret
}

plot.flexsurvreg <- function(x, newdata=NULL, X=NULL, type="survival", fn=NULL, t=NULL, start=NULL,
                             est=TRUE, ci=NULL, B=1000, cl=0.95,
                             col.obs="black", lty.obs=1, lwd.obs=1,
                             col="red",lty=1,lwd=2,
                             col.ci=NULL,lty.ci=2,lwd.ci=1,ylim=NULL,
                             add=FALSE,...)
{
    ## don't calculate or plot CIs by default if all covs are categorical -> multiple curves
    Xraw <- model.frame(x)[,-1,drop=FALSE]
    if (is.null(ci))
        ci <- ((x$ncovs == 0) || (!(sapply(Xraw,is.factor))))
    if (!ci) B <- 0
    summ <- summary.flexsurvreg(x, newdata=newdata, X=X, type=type, fn=fn, t=t, ci=ci, B=B, cl=cl)
    t <- summ[[1]]$time
    X <- if (is.null(attr(summ,"X"))) as.matrix(0, nrow=1, ncol=max(x$ncoveffs,1)) else attr(summ,"X")
    if (is.null(col.ci)) col.ci <- col
    if (is.null(lwd.ci)) lwd.ci <- lwd
    dat <- x$data
    ncovs <- x$ncovs
    isfac <- sapply(Xraw,is.factor)
    if (!is.null(fn)) type <- ""
    if (!add) {
        mm <- as.data.frame(model.matrix(x))
        form <- "Surv(dat$Y[,\"start\"],dat$Y[,\"stop\"],dat$Y[,\"status\"]) ~ "
        form <- paste(form, if (ncovs > 0 && all(isfac)) paste("mm[,",1:x$ncoveffs,"]", collapse=" + ") else 1)
        form <- as.formula(form)
        ## If any continuous covariates, it is hard to define subgroups
        ## so just plot the population survival
        if (type=="survival") {
            plot(survfit(form, data=mm), col=col.obs, lty=lty.obs, lwd=lwd.obs, ...)
        }
        else if (type=="cumhaz") {
            plot(survfit(form, data=mm), fun="cumhaz", col=col.obs, lty=lty.obs, lwd=lwd.obs, ...)
        }
        else if (type=="hazard") {
            muhaz.args <- list(...)[names(list(...)) %in% names(formals(muhaz))]
            plot.args <- list(...)[!names(list(...)) %in% names(formals(muhaz))]
            if (!all(dat$Y[,"start"]==0)) warning("Left-truncated data not supported by muhaz: ignoring truncation point when plotting observed hazard")
            if (any(dat$Y[,"status"] > 1)) stop("Interval-censored data not supported by muhaz")
            if (!all(isfac)){
                haz <- do.call("muhaz", c(list(times=dat$Y[,"stop"], delta=dat$Y[,"status"]), muhaz.args))
                do.call("plot", c(list(haz), list(col=col.obs, lty=lty.obs, lwd=lwd.obs), plot.args))
            }
            else {
                ## plot hazard for all groups defined by unique combinations of covariates
                group <- if(ncovs>0) do.call("interaction", mm) else factor(rep(1,nrow(dat$Y)))
                haz <- list()
                for (i in 1:nrow(X)) {
                    subset <- (group == unique(group)[i])
                    haz[[i]] <- do.call("muhaz", c(list(times=dat$Y[,"time"], delta=dat$Y[,"status"], subset=subset), muhaz.args))
                }
                if (missing(ylim))
                    ylim <- range(sapply(haz, function(x)range(x$haz.est)))
                do.call("plot", c(list(haz[[1]]), list(col=col.obs, lty=lty.obs, lwd=lwd.obs, ylim=ylim), plot.args))
                if (nrow(X)>1) {
                    for (i in 1:nrow(X)) {
                        lines(haz[[i]], col=col.obs, lty=lty.obs, lwd=lwd.obs)
                    }
                }
            }
        }
    }
    col <- rep(col, length=nrow(X)); lty=rep(lty, length=nrow(X)); lwd=rep(lwd, length=nrow(X))
    col.ci <- rep(col.ci, length=nrow(X)); lty.ci=rep(lty.ci, length=nrow(X)); lwd.ci=rep(lwd.ci, length=nrow(X))
    for (i in 1:nrow(X)) {
        if (est) lines(summ[[i]]$time, summ[[i]]$est, col=col[i], lty=lty[i], lwd=lwd[i])
        if (ci) {
            lines(summ[[i]]$time, summ[[i]]$lcl, col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i])
            lines(summ[[i]]$time, summ[[i]]$ucl, col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i])
        }
    }
}

lines.flexsurvreg <- function(x, newdata=NULL, X=NULL, type="survival", t=NULL,
                              est=TRUE, ci=NULL, B=1000, cl=0.95,
                              col="red",lty=1,lwd=2,
                              col.ci=NULL,lty.ci=2,lwd.ci=1, ...)
{
    plot.flexsurvreg(x, newdata=newdata, X=X, type=type, t=t, est=est, ci=ci, B=B, cl=cl,
                     col=col, lty=lty, lwd=lwd, col.ci=col.ci,lty.ci=lty.ci,lwd.ci=lwd.ci, add=TRUE, ...)
}

vcov.flexsurvreg <- function (object, ...)
{
    object$cov
}

.onLoad <- function(libname, pkgname) {
    assign("flexsurv.dfns", new.env(), envir=parent.env(environment()))
}

model.frame.flexsurvreg <- function(formula, ...)
{
    x <- formula
    x$data$m
}

model.matrix.flexsurvreg <- function(object, par=NULL, ...)
{
    x <- object
    if (is.null(par)) compress.model.matrices(x$data$mml) else x$data$mml[[par]]
}

logLik.flexsurvreg <- function(object, ...){
    val <- object$loglik
    attr(val, "df") <- object$npars
    attr(val, "nobs") <- nrow(model.frame(object))
    class(val) <- "logLik"
    val
}

