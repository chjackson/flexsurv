## Second derivatives of the log density and log survival
## with respect to parameters on transformed scales
## The transformed scale is log scale for positive parameters
## Gompertz shape is unrestricted

D2Ldexp <- function(t, rate){
  res <- array(dim=c(length(t), 1, 1))
  parnames <- "rate"
  dimnames(res) <- list(NULL, parnames, parnames)
  ts <- - t*rate
  res[,"rate","rate"] <- ts
  res
}

D2LSexp <- function(t, rate){
  res <- array(dim=c(length(t), 1, 1))
  parnames <- "rate"
  dimnames(res) <- list(NULL, parnames, parnames)
  res[,"rate","rate"] <- -t*rate
  res
}

D2Ldweibull <- function(t, shape, scale){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","scale")
  dimnames(res) <- list(NULL, parnames, parnames)
  tss <- (t/scale)^shape
  lts <- log(t/scale)
  res[,"shape","shape"] <- ifelse(t==0, 0, shape*lts*(1 - tss - shape*lts*tss))
  res[,"scale","scale"] <- -shape^2*tss
  res[,"shape","scale"] <- res[,"scale","shape"] <-
    ifelse(t==0, 0, shape*(tss - 1 + shape*lts*tss))
  res
}

D2LSweibull <- function(t, shape, scale){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","scale")
  dimnames(res) <- list(NULL, parnames, parnames)
  tss <- (t/scale)^shape
  lts <- log(t/scale)
  res[,"shape","shape"] <- ifelse(t==0, 0, -shape*lts*tss*(1 + shape*lts))
  res[,"scale","scale"] <- -shape^2*tss
  res[,"shape","scale"] <-  res[,"scale","shape"] <-
    ifelse(t==0, 0, shape*tss*(1 + shape*lts))
  res
}

D2Ldweibull.quiet <- D2Ldweibull
D2LSweibull.quiet <- D2LSweibull

D2LdweibullPH <- function(t, shape, scale){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","scale")
  dimnames(res) <- list(NULL, parnames, parnames)
  tss <- (t/scale)^shape
  lts <- log(t/scale)
  res[,"shape","shape"] <- ifelse(t==0, 0, shape*log(t)*(1 - scale*t^shape*(1 + shape*log(t))))
  res[,"scale","scale"] <- -scale*t^shape
  res[,"shape","scale"] <-  res[,"scale","shape"] <-
    ifelse(t==0, 0, -scale*shape*log(t)*t^shape)
  res
}

D2LSweibullPH <- function(t, shape, scale){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","scale")
  dimnames(res) <- list(NULL, parnames, parnames)
  res[,"shape","shape"] <- ifelse(t==0, 0, -shape*scale*log(t)*(t^shape + shape*log(t)*t^shape))
  res[,"scale","scale"] <- -scale*t^shape
  res[,"shape","scale"] <-  res[,"scale","shape"] <-
    ifelse(t==0, 0, - scale * shape* log(t) * t^shape)
  res
}

D2Ldgompertz <- function(t, shape, rate){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","rate")
  dimnames(res) <- list(NULL, parnames, parnames)
  res[shape==0,"shape","shape"] <- 0
  res[shape==0,"rate","rate"] <- - rate[shape==0]*t[shape==0]
  res[shape==0,"shape","rate"] <- res[shape==0,"rate","shape"] <-  0
  est <- exp(shape*t)
  sn0 <- (shape!=0)
  t <- t[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
  res[sn0,"shape","shape"] <-  2*rate/shape^3*(1 - est) + 2*rate/shape^2*t*est - rate/shape*t^2*est
  res[sn0,"rate","rate"] <- rate/shape*(1 - est)
  res[sn0,"shape","rate"] <- res[sn0,"rate","shape"] <-
    -rate/shape*(1/shape*(1 - est) + t*est)
  res
}

D2LSgompertz <- function(t, shape, rate){
  res <- array(dim=c(length(t), 2, 2))
  parnames <- c("shape","rate")
  dimnames(res) <- list(NULL, parnames, parnames)
  res[shape==0,"shape","shape"] <- 0
  res[shape==0,"rate","rate"] <- -rate[shape==0]*t[shape==0]
  res[shape==0,"shape","rate"] <- res[shape==0,"rate","shape"] <-  0
  est <- exp(shape*t)
  sn0 <- (shape!=0)
  t <- t[sn0]; rate <- rate[sn0]; shape <- shape[sn0]
  res[sn0,"shape","shape"] <- 2*rate/shape^3*(1 - est) + 2*rate/shape^2*t*est - rate/shape*t^2*est
  res[sn0,"rate","rate"] <- rate/shape*(1 - est)
  res[sn0,"shape","rate"] <- res[sn0,"rate","shape"] <- 
    -rate/shape*(1/shape*(1 - est) + t*est)
  res
}

D2Ldsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", spline="rp"){
    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE, spline=spline)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q
    b <- basis(knots, tsfn(t,timescale), spline=spline)
    db <- dbasis(knots, tsfn(t,timescale), spline=spline)
    eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
    ds <- rowSums(db * gamma)
    if (scale=="odds") {
      eeta <- 2*exp(eta)/(1 + exp(eta))^2
      db <- dbasis(knots, tsfn(t,timescale), spline=spline)
    }
    npars <- ncol(gamma)
    parnames <- paste0("gamma",seq_len(npars)-1)
    res <- array(dim=c(length(t), npars, npars),
                 dimnames = list(NULL, parnames, parnames))
    for (i in 1:npars){
      for (j in 1:npars){
        if (scale=="hazard") 
          res[ind,i,j] <- -db[,i]*db[,j] / ds^2 - b[,i]*b[,j] * exp(eta)
        else if (scale=="odds") 
          res[ind,i,j] <- -db[,i]*db[,j] / ds^2 - b[,i]*b[,j] * eeta
      }
    }
    res
}

D2LSsurvspline <- function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", spline="rp"){
    d <- dbase.survspline(q=t, gamma=gamma, knots=knots, scale=scale, deriv=TRUE, spline=spline)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); t <- q
    b <- basis(knots, tsfn(t,timescale), spline=spline)
    npars <- ncol(gamma)
    parnames <- paste0("gamma",seq_len(npars)-1)
    res <- array(dim=c(length(t), npars, npars),
                 dimnames = list(NULL, parnames, parnames))
    if (length(t) > 0){
      eta <- rowSums(b * gamma) + as.numeric(X %*% beta)
      if (scale=="odds") {
        eeta <- exp(eta)/(1 + exp(eta))^2
      }
      for (i in 1:npars){
        for (j in 1:npars){
          if (scale=="hazard") 
            res[ind,i,j] <- ifelse(t==0, 0, - b[,i]*b[,j]*exp(eta))
          else if (scale=="odds") 
            res[ind,i,j] <- ifelse(t==0, 0, - b[,i]*b[,j] * eeta)
        }
      }
    }
    res
}

d2deriv <- function(d2dfn, ddcall, X, mx, dlist){
  res.base <- do.call(d2dfn, ddcall)
  res.beta <- D2cov(res.base, X, mx, dlist) # does this have to be in a different function? 
  nobs <- dim(res.base)[1]
  nbasepars <- dim(res.base)[2]
  ncoveffs <- dim(res.beta[[1]])[2]
  npars <- nbasepars + ncoveffs
  res <- array(dim = c(nobs, npars, npars))
  if (length(ddcall$t) > 0){
    dimnames(res)[[2]] <- dimnames(res)[[3]] <- c(colnames(res.base), colnames(X))
    baseinds <- seq_len(nbasepars)
    res[,baseinds,baseinds] <- res.base
    betainds <- nbasepars + seq_len(ncoveffs)
    res[,betainds,betainds] <- res.beta$ccres
    res[,baseinds,betainds] <- res.beta$bcres
    res[,betainds,baseinds] <- aperm(res.beta$bcres,c(1,3,2))
  }
  res
}

D2cov <- function(res, X, mx, dlist){
    ncoveffs <- length(unlist(mx))
    nbpars <- length(mx)
    bcres  <- array(dim = c(nrow(res), nbpars, ncoveffs))
    basepars <- dimnames(res)[[2]] # ordered as in the distribution (used for baseline pars)
    basepars_locfirst <- names(mx) # ordered with location parameter first (used for covariate effects)
    dimnames(bcres)[[2]] <- basepars_locfirst
    inds <- c(0,cumsum(sapply(mx[basepars_locfirst],length)))
    ## second derivs wrt two covariate effects
    ccres <- array(dim = c(nrow(res), ncoveffs, ncoveffs))
    for (i in seq_len(nbpars)) {
      mi <- mx[[basepars_locfirst[i]]]
      for (j in seq_along(mi)){
          for (k in seq_len(nbpars)) {
            mk <- mx[[basepars_locfirst[k]]]
            for (l in seq_along(mk)){
              ccres[,inds[i]+j,inds[k]+l] <-
                X[,mi[j]]*X[,mk[l]]*res[,basepars_locfirst[i],basepars_locfirst[k]]
            }
          }
        }
    }
    ## second derivs wrt baseline par x covariate effect
    for (i in seq_len(nbpars)) {
      for (j in seq_len(nbpars)) {
        mj <- mx[[basepars_locfirst[j]]]
        for (k in seq_along(mj)){
          bcres[,i,inds[j]+k] <- X[,mj[k]] * res[,basepars[i],basepars_locfirst[j]]
        }
      }
    }
    list(ccres=ccres, bcres=bcres)
}


D2minusloglik.flexsurv <- function(optpars, Y, X=0, weights, bhazard, rtrunc, dlist, inits, dfns, aux, mx, fixedpars=NULL) {
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

    dd <- d2deriv(dfns$D2Ld, ddcall, X[dead,,drop=FALSE], mx, dlist)
    dscens <- d2deriv(dfns$D2LS, dsccall, X[!dead,,drop=FALSE], mx, dlist)
    if (sum(dead) > 0) dd <- dd * weights[dead]
    if (sum(!dead) > 0) dscens <- dscens * weights[!dead]
    dstrunc <- d2deriv(dfns$D2LS, dstcall, X, mx, dlist) * weights
    res <- - ( apply(dd,2:3,sum) + apply(dscens,2:3,sum) - apply(dstrunc,2:3,sum) )

    if (any(bhazard > 0)) { 
        dcall <- ddcall
        dcall$x <- ddcall$t; dcall$t <- NULL
        dens <- do.call(dfns$d, dcall)
        pcall <- dcall
        pcall$q <- pcall$x; pcall$x <- NULL
        surv <- 1 - do.call(dfns$p, pcall)
        haz <- dens / surv
        bw <- bhazard[dead]
        offseti <- 1 / (1 + bw/haz)
        d1Ld <- dderiv(dfns$DLd, ddcall, X[dead,,drop=FALSE], mx, dlist)
        d2Ld <- d2deriv(dfns$D2Ld, ddcall, X[dead,,drop=FALSE], mx, dlist)
        d1LS <- dderiv(dfns$DLS, ddcall, X[dead,,drop=FALSE], mx, dlist)
        d2LS <- d2deriv(dfns$D2LS, ddcall, X[dead,,drop=FALSE], mx, dlist)
        dhazinv <- surv/dens*(d1LS - d1Ld)
        
        if (any(dead)){
            ## given two matrices with dims (n,p), 
            ## return array of dim(n,p,p) with outer product of each pair of matrix rows
            vouter <- function(y,z){
                rows <- seq_len(nrow(y))
                res <- mapply(outer, split(y, rows), split(z, rows), SIMPLIFY=FALSE)
                res <- simplify2array(res, except = NULL)
                aperm(res, c(3,1,2))
            }
            
            doff <- - offseti*bw*(offseti*vouter(dhazinv, dhazinv)*bw +
                                      (d2Ld - d2LS)*surv/dens +
                                      vouter(d1Ld - d1LS, dhazinv))
            
            res <- res - apply(doff*weights[dead],2:3,sum)
        }
    }
    ## currently wastefully calculates derivs for fixed pars then discards them
    optpars <- setdiff(1:npars, fixedpars)
    res[optpars,optpars]
}

hess.test <- function(optpars, Y, X, weights, bhazard, rtrunc,
                      dlist, inits, dfns, aux, mx, fixedpars){
  an.d <- D2minusloglik.flexsurv(optpars=optpars, Y=Y, X=X,
                                 weights=weights, bhazard=bhazard,
                                 rtrunc=rtrunc, dlist=dlist,
                                 inits=inits, dfns=dfns, aux=aux,
                                 mx=mx, fixedpars=fixedpars)
  if (requireNamespace("numDeriv", quietly = TRUE))
    num.d <- numDeriv::hessian(minusloglik.flexsurv, optpars, Y=Y,
                               X=X, weights=weights, bhazard=bhazard,
                               rtrunc=rtrunc, dlist=dlist,
                               inits=inits, dfns=dfns, aux=aux, mx=mx,
                               fixedpars=fixedpars)
  else stop("\"numDeriv\" package not available")
  res <- list(analytic=an.d, numeric=num.d)
  list(res=res, error=mean(abs(an.d - num.d)))
}
