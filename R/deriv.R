## Deriv of loglik wrt transformed parameters p
## loglik(p|x) = sum(log(f(p|xobs))) + sum(log(S(p|xcens))) - sum(log(S(p|xtrunc)))
## dloglik/dp  = sum (df/dp / f(p)) | xobs) + sum(dS/dp / S(p) | xcens) - sum(dS/dp / S(p) | xtrunc)
##             = sum(dlogf/dp | xobs) + sum(dlogS/dp | xcens) - sum(dlogS/dp | xtrunc)

Dminusloglik.flexsurv <- function(optpars, Y, X=0, weights, bhazard, dlist, inits, dfns, aux, mx, fixedpars=NULL) {
    pars <- inits
    npars <- length(pars)
    pars[setdiff(1:npars, fixedpars)] <- optpars
    nbpars <- length(dlist$pars)
    pars <- as.list(pars)
    ncovs <- length(pars) - length(dlist$pars)
    if (ncovs > 0)
        beta <- unlist(pars[(nbpars+1):npars])
    for (i in dlist$pars) {
        if (length(mx[[i]]) > 0)
            pars[[i]] <- pars[[i]] + X[,mx[[i]],drop=FALSE] %*% beta[mx[[i]]]
        else
            pars[[i]] <- rep(pars[[i]], length(Y[,"stop"]))
    }
    dead <- Y[,"status"]==1
    ddcall <- list(t=Y[dead,"stop"])
    dsccall <- list(t=Y[!dead,"stop"])
    dstcall <- list(t=Y[,"start"])
    for (i in 1:nbpars)
        ddcall[[names(pars)[i]]] <-
            dsccall[[names(pars)[i]]] <-
                dstcall[[names(pars)[i]]] <-
                    dlist$inv.transforms[[i]](pars[[i]])
    for (i in seq_along(aux)){
        ddcall[[names(aux)[i]]] <- dsccall[[names(aux)[i]]] <-
            dstcall[[names(aux)[i]]] <- aux[[i]]
    }
    for (i in dlist$pars) {
        ddcall[[i]] <- ddcall[[i]][dead]
        dsccall[[i]] <- dsccall[[i]][!dead]
    }
    dd <- dderiv(dfns$DLd, ddcall, X[dead,,drop=FALSE], mx, dlist)
    dscens <- dderiv(dfns$DLS, dsccall, X[!dead,,drop=FALSE], mx, dlist)
    if (sum(dead) > 0) dd <- dd * weights[dead]
    if (sum(!dead) > 0) dscens <- dscens * weights[!dead]
    dstrunc <- dderiv(dfns$DLS, dstcall, X, mx, dlist) * weights
    res <- - ( colSums(dd) + colSums(dscens) - colSums(dstrunc) )
    ## currently wastefully calculates derivs for fixed pars then discards them
    res[setdiff(1:npars, fixedpars)]
}

dderiv <- function(ddfn, ddcall, X, mx, dlist){
    if (length(ddcall$t) == 0) matrix(0) else { 
        res.base <- do.call(ddfn, ddcall)
        res.beta <- Dcov(res.base, X, mx, dlist)
        cbind(res.base, res.beta)
    }
}

## Derivatives of log density with respect to covariate effects are
## just found by an easy chain rule given the baseline derivatives and
## the covariate value: the parameter on the real-line scale is always
## a linear function of covariates.

Dcov <- function(res, X, mx, dlist){
    ncovs <- length(unlist(mx))
    cres <- matrix(nrow=nrow(res), ncol=ncovs)
    inds <- c(0,cumsum(sapply(mx,length)))
    ## parameters of res ordered as distribution definition, but mx ordered with location first
    for (i in seq_along(mx)) {
        for (j in seq_along(mx[[i]])){
            cres[,inds[i]+j] <- X[,mx[[i]][j]]*res[,match(names(mx)[i], dlist$pars)]
        }
    }
    cres
}

## Derivatives of log density and log survival probability with
## respect to baseline parameters for various distributions.  Baseline
## parameters are on the real-line scale, commonly the log scale.  No
## easy derivatives available for other distributions.

## Exponential

DLdexp <- function(t, rate){
    res <- matrix(nrow=length(t), ncol=1)
    ts <- 1 - t*rate
    res[,1] <- ts
    res
}

DLSexp <- function(t, rate){
    res <- matrix(nrow=length(t), ncol=1)
    res[,1] <- -t*rate
    res
}

## Weibull

DLdweibull <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    res[,1] <- 1 + shape*(log(t/scale) - log(t/scale)*tss) # wrt log shape
    res[,2] <- -1 - (shape-1) + shape*tss
    res
}

DLSweibull <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    res[,1] <- ifelse(t==0, 0, -shape*log(t/scale)*tss)
    res[,2] <- tss*shape
    res
}

DLdweibull.quiet <- DLdweibull
DLSweibull.quiet <- DLSweibull

## Gompertz

DLdgompertz <- function(t, shape, rate){
    res <- matrix(nrow=length(t), ncol=2)
    rs <- rate/shape*exp(shape*t)
    res[shape==0,1] <- 0
    res[shape==0,2] <- 1 - rate[shape==0] * t[shape==0]
    sn0 <- (shape!=0)
    t <- t[sn0]; rs <- rs[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
    res[shape!=0,1] <- t + rs*(1/shape - t) - rate/shape^2
    res[shape!=0,2] <- 1 - rs + rate/shape
    res
}

DLSgompertz <- function(t, shape, rate){
    res <- matrix(nrow=length(t), ncol=2)
    rs <- rate/shape*exp(shape*t)
    res[shape==0,1] <- 0
    res[shape==0,2] <- - rate[shape==0] * t[shape==0]
    sn0 <- (shape!=0)
    t <- t[sn0]; rs <- rs[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
    res[shape!=0,1] <- rs*(1/shape - t) - rate/shape^2
    res[shape!=0,2] <-  - rs + rate/shape
    res
}

DLdsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log"){
    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q

# TODO special value handling
    
    b <- basis(knots, tsfn(t,timescale))
    db <- dbasis(knots, tsfn(t,timescale))
    eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
    ds <- rowSums(db * gamma)
    for (i in 1:ncol(gamma)){
        if (scale=="hazard")
            ret[ind,i] <- db[,i] / ds + b[,i] * (1 - exp(eta))
        else if (scale=="odds"){
            eeta <- 1 - 2*exp(eta)/(1 + exp(eta))
            ret[ind,i] <- db[,i] / ds + b[,i] * eeta
        }
    }
    ret
}

DLSsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log"){

    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q

# TODO special value handling   
    b <- basis(knots, tsfn(t,timescale))
    eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
    for (i in 1:ncol(gamma)){
        if (scale=="hazard")
            ret[ind,i] <- ifelse(t==0, 0, - b[,i] * exp(eta))
        else if (scale=="odds"){
            eeta <- exp(eta)/(1 + exp(eta))
            ret[ind,i] <- ifelse(t==0, 0, - b[,i] * eeta)
        }
    }
    ret
}

deriv.test <- function(optpars, Y, X, weights, bhazard, dlist, inits, dfns, aux, mx, fixedpars){
    an.d <- Dminusloglik.flexsurv(optpars, Y, X, weights, bhazard, dlist, inits, dfns, aux, mx, fixedpars)
    if (requireNamespace("numDeriv", quietly = TRUE))
        num.d <- numDeriv::grad(minusloglik.flexsurv, optpars, Y=Y, X=X, weights=weights, bhazard=bhazard,
                                dlist=dlist, inits=inits, dfns=dfns, aux=aux, mx=mx, fixedpars=fixedpars)
    else stop("\"numDeriv\" package not available")
    res <- cbind(analytic=an.d, numeric=as.vector(num.d))
    rownames(res) <- names(optpars)
    list(res=res, error=mean(abs(an.d - num.d)))
}
