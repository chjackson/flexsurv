##' Flexible parametric models for right-truncated, uncensored data defined by times of initial and final events.
##'
##' This function estimates the distribution of the time between an initial and final event, in situations where individuals are only observed if they have experienced both events before a certain time, thus they are right-truncated at this time.   The time of the initial event provides information about the time from initial to final event, given the truncated observation scheme, and initial events are assumed to occur with an exponential growth rate.
##'
##' Covariates are not currently supported.
##'
##' Note that \code{\link{flexsurvreg}}, with an \code{rtrunc} argument, can fit models for a similar kind of data, but those models ignore the information provided by the time of the initial event.
##'
##' A nonparametric estimator of survival under right-truncation is also provided in \code{\link{survrtrunc}}.  See Seaman et al. (2020) for a full comparison of the alternative models. 
##' 
##' 
##' @param t Vector of time differences between an initial and final event for a set of individuals.
##'
##' @param tinit Absolute time of the initial event for each individual.
##'
##' @param rtrunc Individual-specific right truncation points on the same scale as \code{t}, so that each individual's survival time \code{t} would not have been observed if it was greater than the corresponding element of \code{rtrunc}.  Only used in \code{method="joint"}.  In \code{method="final"}, the right-truncation is implicit.
##'
##' @param tmax Maximum possible time between initial and final events that could have been observed.  This is only used in \code{method="joint"}, and is ignored in \code{method="final"}. 
##'
##' @param data Data frame containing \code{t}, \code{rtrunc} and \code{tinit}.
##'
##' @param method If \code{"joint"} then the "joint-conditional" method is used.  If \code{"final"} then the "conditional-on-final" method is used.   The "conditional-on-initial" method can be implemented by using \code{\link{flexsurvreg}} with a \code{rtrunc} argument.  These methods are all described in Seaman et al. (2020).
##'
##' @param theta Initial value (or fixed value) for the exponential growth rate \code{theta}. Defaults to 1.
##'
##' @param inits Initial values for the parameters of the parametric survival distributon. If not supplied, a heuristic is used. as is done in \code{\link{flexsurvreg}}.
##'
##' @param fixedpars Integer indices of the parameters of the survival distribution that should be fixed to their values supplied in \code{inits}.   Same length as \code{inits}.
##'
##' @param fixed.theta Should \code{theta} be fixed at its initial value or estimated.  This only applies to \code{method="joint"}.  For \code{method="final"}, \code{theta} must be fixed.
##'
##' @param optim_control List to supply as the \code{control} argument to \code{\link{optim}} to control the likelihood maximisation.
##'
##' @inheritParams flexsurvreg
##'
##' @seealso \code{\link{flexsurvreg}}, \code{\link{survrtrunc}}.
##' 
##' @examples
##' set.seed(1) 
##' ## simulate time to initial event
##' X <- rexp(1000, 0.2)
##' ## simulate time between initial and final event
##' T <- rgamma(1000, 2, 10) 
##'
##' ## right-truncate to keep only those with final event time
##' ## before tmax
##' tmax <- 40
##' obs <- X + T < tmax 
##' rtrunc <- tmax - X
##' dat <- data.frame(X, T, rtrunc)[obs,]
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gamma", theta=0.2)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gamma", theta=0.2, fixed.theta=FALSE)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gamma", theta=0.2, inits=c(1, 8))
##' 
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gamma", theta=0.2, method="final")
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gamma", fixed.theta=TRUE)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="weibull", fixed.theta=TRUE)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="lnorm", fixed.theta=TRUE)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gengamma", fixed.theta=TRUE)
##'
##' flexsurvrtrunc(t=T, rtrunc=rtrunc, tinit=X, tmax=40, data=dat,
##'                 dist="gompertz", fixed.theta=TRUE)
##'
##' @references Seaman, S., Presanis, A. and Jackson, C. (2020) Estimating a Time-to-Event
##' Distribution from Right-Truncated Data in an Epidemic: a Review of Methods
##'
##' @export
flexsurvrtrunc <- function(t, tinit, rtrunc, tmax, data=NULL, method="joint", dist, 
                           theta=NULL, fixed.theta=TRUE,
                           inits=NULL, fixedpars=NULL, dfns=NULL,
                           integ.opts=NULL, cl=0.95,
                           optim_control = list()){

    tdelay <- eval(substitute(t), data, parent.frame())
    rtrunc <- eval(substitute(rtrunc), data, parent.frame())
    tinit <- eval(substitute(tinit), data, parent.frame())
    dlist <- parse.dist(dist)
    dfns <- form.dp(dlist, dfns, integ.opts)

    if (is.null(inits)){
        if (is.null(theta)) theta <- 1

        ## FIXME for aux arguments.  pop out a function for this 
        ## what is this for.  to fit a survreg weibull and use their heiurisitic
        ## shdnt be a problem, just build the formula up 
        build_inits_aux <- function() {
            list(forms = list(Surv(tdelay) ~ 1),
                 control = list(),
                 counting = FALSE
                 )
        } 
        aux <- build_inits_aux()
        dlist$inits <- expand.inits.args(dlist$inits)
        initsd <- dlist$inits(t=tdelay, aux=aux)
        names(initsd) <- dlist$pars
        inits <- c(initsd, theta=theta)
    } else { 
        if (!is.numeric(inits)) stop("`inits` should be numeric")
        if (length(inits) != length(dlist$pars)) stop(sprintf("%s values supplied in `inits`, but %s parameters in distribution \"%s\"", length(inits), length(dlist$pars), dlist$name))
        inits <- c(inits, theta)
    }
    npars <- length(inits)
    nbpars <- npars - 1
    if (method=="final") {
        if (!fixed.theta) warning("theta cannot be estimated in method=\"final\"")
        fixed.theta <- TRUE
    }

    check.fixedpars(fixedpars, nbpars)
    if (isTRUE(fixedpars)) fixedpars <- 1:nbpars
    
    if (fixed.theta) fixedpars <- unique(c(fixedpars, npars))
    names(inits) <- c(dlist$pars, "theta")
    for (i in 1:nbpars){
        inits[i] <- dlist$transforms[[i]](inits[i])
    }
    inits[nbpars+1] <- inits[nbpars+1]
    optpars <- setdiff(1:npars, fixedpars)
    optvals <- inits[optpars]

    fixed <- length(fixedpars)==npars
    if (!fixed) { 
        opt <- optim(optvals,
                     frtrunc_loglik_factory(inits=inits,tinit=tinit, tdelay=tdelay, tmax=tmax,
                                            method=method, dlist=dlist, dfns=dfns, fixedpars=fixedpars),
                     control=optim_control, hessian=TRUE) 
        est <- opt$par
        covar <- solve(opt$hessian)
        se <- sqrt(diag(covar))
        loglik <- - opt$value
    } else {
        opt <- NULL
        est <- inits
        se <- covar <- NA
        loglik <- - frtrunc_loglik_factory(inits=inits, tinit=tinit, tdelay=tdelay, tmax=tmax,
                                 method=method,dlist=dlist, dfns=dfns, fixedpars=fixedpars)(inits)
    }
    if (!is.numeric(cl) || length(cl)>1 || !(cl>0) || !(cl<1))
                stop("cl must be a number in (0,1)")
    lcl <- est - qnorm(1 - (1-cl)/2)*se
    ucl <- est + qnorm(1 - (1-cl)/2)*se
    est.opt <- cbind(estlog=est, selog=se, lcllog=lcl, ucllog=ucl)

    res <- matrix(nrow=npars, ncol=ncol(est.opt),
                  dimnames=list(names(inits), colnames(est.opt)))
    res[optpars,] <- est.opt
    res[fixedpars,"estlog"] <- inits[fixedpars]
    res <- as.data.frame(res)
    res$ucl <- res$lcl <- res$est <- numeric(npars)
    for (i in 1:npars){
        trf <- if (i==npars) exp else dlist$inv.transforms[[i]]
        res$est[i] <- trf(res$estlog[i])
        res$lcl[i] <- trf(res$lcllog[i])
        res$ucl[i] <- trf(res$ucllog[i])
    }    
    res$pars <- c(dlist$pars, "theta")
    est <- res[,c("pars","est","lcl","ucl","estlog","selog")]
    npars <- length(optpars)
    res <- list(call=match.call(),
                res=est,
                loglik=loglik,
                data=data.frame(t=tdelay, tinit=tinit, rtrunc=rtrunc),
                dfns=dfns,
                dlist=dlist,
                ncovs=0,
                cov=covar,
                opt=opt,
                optpars=optpars,
                fixedpars=fixedpars,
                npars=npars,
                AIC = -2*loglik + 2*npars,
                res.t=est)
    class(res) <- "flexsurvrtrunc"
    res
}

##' @export
print.flexsurvrtrunc <- function(x, ...)
{
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    if (x$npars > 0) {
        res <- x$res
        cat ("Estimates: \n")
        args <- list(...)
        if (is.null(args$digits)) args$digits <- 3
        f <- do.call("format", c(list(x=res), args))
        print(f, print.gap=2, quote=FALSE, na.print="")
    }
    cat("\nLog-likelihood = ", x$loglik, ", df = ", x$npars,
        "\nAIC = ", x$AIC, "\n\n", sep="")
}



frtrunc_loglik_factory <- function(inits, tinit, tdelay, tmax, method, dlist, dfns, fixedpars){
    pars <- inits
    insert.locations <- setdiff(seq_along(pars), fixedpars)

    frtrunc_intfn_joint <- function(a, dpars, theta) {
        pargs <- as.list(dpars)
        pargs$q <- tmax - a
        exp(theta*a) * do.call(dfns$p, pargs)
        
    }

    frtrunc_intfn_final <- function(a, dpars, theta, log=FALSE){
        dargs <- as.list(dpars)
        dargs$x <- a
        dargs$log <- log
        ret <- do.call(dfns$d, dargs)
        if (log)
            ret <- ret - theta*a
        else
            ret <- ret * exp(  - theta*a)
        ret
    }

    function(optpars, ...) {
        pars[insert.locations] <- optpars
        nbpars <- length(pars) - 1  ## NO COVARIATES SO FAR . this needed if covariates
        for (i in 1:nbpars){
            pars[i] <- dlist$inv.transforms[[i]](pars[i])
        }
        dpars <- pars[1:nbpars]
        theta <- pars[nbpars+1]
        nv <- length(tinit)
        dargs <- as.list(dpars)
        dargs$x <- tdelay 
        dargs$log <- TRUE 
        if (method=="joint") { 
            log.numer <- theta*tinit + do.call(dfns$d, dargs)
            integ <- integrate(frtrunc_intfn_joint, lower=0, upper=tmax, dpars = dpars, theta=theta)$value
            logl <- sum(log.numer) - log(integ) * nv
        } else if (method=="final") {
            if (dlist$name=="gamma") {
                logl <- loglik_final_gamma(tdelay, tinit, dpars, theta)
            } else {
                log.numer <- frtrunc_intfn_final(tdelay, dpars, theta, log=TRUE)
                integ <- numeric(nv)
                y <- tinit + tdelay
                for (i in 1:nv){
                    integ[i] <- integrate(frtrunc_intfn_final, lower=0, upper=y[i],
                                          dpars = dpars, theta=theta)$value
                }
                logl <- sum(log.numer) - sum(log(integ))
            }
        }
        - logl
    }
}

loglik_final_gamma <- function(tdelay, tinit, dpars, theta){
    sum(dgamma(tdelay, dpars["shape"], dpars["rate"]+theta, log=TRUE) -
        pgamma(tinit + tdelay, dpars["shape"], dpars["rate"]+theta, log.p=TRUE))
}
