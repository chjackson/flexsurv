## Deriv of loglik wrt transformed parameters p
## loglik(p|x) = sum(log(f(p|xobs))) + sum(log(S(p|xcens))) - sum(log(S(p|xtrunc)))
## dloglik/dp  = sum (df/dp / f(p)) | xobs) + sum(dS/dp / S(p) | xcens) - sum(dS/dp / S(p) | xtrunc)
##             = sum(dlogf/dp | xobs) + sum(dlogS/dp | xcens) - sum(dlogS/dp | xtrunc)

Dminusloglik.flexsurv <- function(optpars, Y, X=0, weights, bhazard, rtrunc, dlist, inits, dfns, aux, mx, fixedpars=NULL) {
    pars <- inits
    npars <- length(pars)
    pars[setdiff(seq_len(npars), fixedpars)] <- optpars
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
    nobs <- nrow(Y)
    dloglik <- matrix(nrow=nobs, ncol=npars)
    dloglik[dead,] <- dderiv(dfns$DLd, ddcall, X[dead,,drop=FALSE], mx, dlist)
    dloglik[!dead,] <- dderiv(dfns$DLS, dsccall, X[!dead,,drop=FALSE], mx, dlist)
    dstrunc <- dderiv(dfns$DLS, dstcall, X, mx, dlist)
    dloglik <- dloglik - dstrunc

    ## Derivatives with baseline hazard.
    ## adjust dd, dscens and dstrunc.  cens and trunc terms don't depend on pars 
    ## Add deriv of log(1 + (1/h) * bh/w)  just for uncensored event 
    ## 1 / (1 + (1/h)*bh/w) * -h^{-2}bh/w * dh 
    ## and h = f * S^-1, so dh = df*s^-1   +   f * -s^-2 * ds)

    if (any(bhazard > 0)) { 
        dcall <- ddcall
        dcall$x <- ddcall$t; dcall$t <- NULL
        dens <- do.call(dfns$d, dcall)
        pcall <- dcall
        pcall$q <- pcall$x; pcall$x <- NULL
        surv <- 1 - do.call(dfns$p, pcall)
        haz <- dens / surv
        offseti <- 1 / (1 + bhazard[dead] / haz)
        ## Note d/dx log(x) = 1/x ddx
        dscense <- dderiv(dfns$DLS, ddcall, X[dead,,drop=FALSE], mx, dlist)
        doff <- - offseti*bhazard[dead] * (dloglik[dead,]  -  dscense)/ haz
        dloglik[dead,] <- dloglik[dead,] + doff
    }
    res <- - colSums(dloglik*weights)
    ## currently wastefully calculates derivs for fixed pars then discards them
    res[setdiff(1:npars, fixedpars)]
}

dderiv <- function(ddfn, ddcall, X, mx, dlist){
    if (length(ddcall$t) == 0) array(dim=c(0,length(dlist$pars))) else { 
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
    ncoveffs <- length(unlist(mx))
    cres <- matrix(nrow=nrow(res), ncol=ncoveffs)
    inds <- c(0,cumsum(sapply(mx,length)))
    ## parameters of res ordered as distribution definition, but mx ordered with location first
    ## mx is a list, with one component per survival distribution parameter
    ## ith component contains indexes of cols of X that are included in regression for ith parameter
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

## Naming: names refer to natural scale pars, but derivatives are taken with respect to transformed scale pars
## Transformed scale is log scale for positive pars.   Unrestricted pars (e.g. Gompertz shape) not transformed

## Exponential

DLdexp <- function(t, rate){
    res <- matrix(nrow=length(t), ncol=1)
    colnames(res) <- c("rate")
    ts <- 1 - t*rate
    res[,"rate"] <- ts
    res
}

DLSexp <- function(t, rate){
    res <- matrix(nrow=length(t), ncol=1)
    colnames(res) <- c("rate")
    res[,"rate"] <- -t*rate
    res
}

## Weibull accelerated failure time

DLdweibull <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","scale")
    tss <- (t/scale)^shape
    lts <- log(t/scale)
    res[,"shape"] <- 1 + shape*lts*(1 - tss)
    res[,"scale"] <- -1 - (shape-1) + shape*tss
    res
}

DLSweibull <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","scale")
    tss <- (t/scale)^shape
    res[,"shape"] <- ifelse(t==0, 0, -shape*log(t/scale)*tss)
    res[,"scale"] <- tss*shape
    res
}

DLdweibull.quiet <- DLdweibull
DLSweibull.quiet <- DLSweibull

## Weibull proportional hazards

DLdweibullPH <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","scale")
    res[,"shape"] <- 1 + shape*log(t)*(1 - scale*t^shape)
    res[,"scale"] <- 1 - scale*t^shape
    res
}

DLSweibullPH <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","scale")
    res[,"shape"] <- ifelse(t==0, 0, -scale*shape*log(t)*t^shape)
    res[,"scale"] <- -scale*t^shape 
    res
}

## Gompertz (note this is derivative wrt shape and log rate)

DLdgompertz <- function(t, shape, rate){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","rate")
    res[shape==0,"shape"] <- 0
    res[shape==0,"rate"] <- 1 - rate[shape==0] * t[shape==0]
    sn0 <- (shape!=0)
    t <- t[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
    est <- exp(shape*t)
    res[sn0,"shape"] <- t + -rate/shape*(1/shape*(1 - est) + t*est)
    res[sn0,"rate"] <- 1 + rate/shape*(1 - est)
    res
}

DLSgompertz <- function(t, shape, rate){
    res <- matrix(nrow=length(t), ncol=2)
    colnames(res) <- c("shape","rate")
    res[shape==0,"shape"] <- 0
    res[shape==0,"rate"] <- - rate[shape==0] * t[shape==0]
    sn0 <- (shape!=0)
    t <- t[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
    est <- exp(shape*t)
    res[sn0,"shape"] <- -rate/shape*(1/shape*(1 - est) + t*est)
    res[sn0,"rate"] <- rate/shape*(1 - est)
    res
}

DLdsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", spline="rp"){
    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE, spline=spline)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q
    b <- basis(knots, tsfn(t,timescale), spline=spline)
    db <- dbasis(knots, tsfn(t,timescale), spline=spline)
    eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
    ds <- rowSums(db * gamma)
    npars <- ncol(gamma)
    parnames <- paste0("gamma", seq_len(npars)-1)
    colnames(ret) <- parnames
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

DLSsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", spline="rp"){
    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE, spline=spline)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q
    b <- basis(knots, tsfn(t,timescale), spline=spline)
    if (any(ind)){
      eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
      for (i in 1:ncol(gamma)){
        if (scale=="hazard")
          ret[ind,i] <- ifelse(t==0, 0, - b[,i] * exp(eta))
        else if (scale=="odds"){
          eeta <- exp(eta)/(1 + exp(eta))
          ret[ind,i] <- ifelse(t==0, 0, - b[,i] * eeta)
        }
      }
    }
    ret
}

deriv.test <- function(optpars, Y, X, weights, bhazard, rtrunc, dlist, inits, dfns, aux, mx, fixedpars){
  an.d <- Dminusloglik.flexsurv(optpars=optpars, Y=Y, X=X, weights=weights, bhazard=bhazard,
                                rtrunc=rtrunc, dlist=dlist, inits=inits,
                                dfns=dfns, aux=aux, mx=mx, fixedpars=fixedpars)
    if (requireNamespace("numDeriv", quietly = TRUE))
      num.d <- numDeriv::grad(minusloglik.flexsurv, optpars, Y=Y, X=X, weights=weights, bhazard=bhazard,
                              rtrunc=rtrunc, 
                              dlist=dlist, inits=inits, dfns=dfns, aux=aux, mx=mx, fixedpars=fixedpars)
    else stop("\"numDeriv\" package not available")
    res <- cbind(analytic=an.d, numeric=as.vector(num.d))
    rownames(res) <- names(optpars)
    list(res=res, error=mean(abs(an.d - num.d)))
}
