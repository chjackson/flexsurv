##' Royston/Parmar spline survival distribution
##' 
##' Probability density, distribution, quantile, random generation, hazard
##' cumulative hazard, mean and restricted mean functions for the Royston/Parmar
##' spline model.
##' 
##' @aliases dsurvspline psurvspline qsurvspline rsurvspline
##' hsurvspline Hsurvspline mean_survspline rmst_survspline
##' @param x,q,t Vector of times.
##' @param p Vector of probabilities.
##' @param n Number of random numbers to simulate.
##' @param gamma Parameters describing the baseline spline function, as
##' described in \code{\link{flexsurvspline}}.  This may be supplied as a
##' vector with number of elements equal to the length of \code{knots}, in
##' which case the parameters are common to all times.  Alternatively a matrix
##' may be supplied, with rows corresponding to different times, and columns
##' corresponding to \code{knots}.
##' @param start Optional left-truncation time or times.  The returned
##' restricted mean survival will be conditioned on survival up to
##' this time.
##' @param beta Vector of covariate effects (deprecated).
##' @param X Matrix of covariate values (deprecated).
##' @param knots Locations of knots on the axis of log time, supplied in
##' increasing order.  Unlike in \code{\link{flexsurvspline}}, these include
##' the two boundary knots.  If there are no additional knots, the boundary
##' locations are not used.  If there are one or more additional knots, the
##' boundary knots should be at or beyond the minimum and maximum values of the
##' log times.  In \code{\link{flexsurvspline}} these are exactly at the
##' minimum and maximum values.
##' 
##' This may in principle be supplied as a matrix, in the same way as for
##' \code{gamma}, but in most applications the knots will be fixed.
##' @param scale \code{"hazard"}, \code{"odds"}, or \code{"normal"}, as
##' described in \code{\link{flexsurvspline}}.  With the default of no knots in
##' addition to the boundaries, this model reduces to the Weibull, log-logistic
##' and log-normal respectively.  The scale must be common to all times.
##' @param timescale \code{"log"} or \code{"identity"} as described in
##' \code{\link{flexsurvspline}}.
##' @param offset An extra constant to add to the linear predictor
##' \eqn{\eta}{eta}.
##' @param log,log.p Return log density or probability.
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##' @return \code{dsurvspline} gives the density, \code{psurvspline} gives the
##' distribution function, \code{hsurvspline} gives the hazard and
##' \code{Hsurvspline} gives the cumulative hazard, as described in
##' \code{\link{flexsurvspline}}.
##' 
##' \code{qsurvspline} gives the quantile function, which is computed by crude
##' numerical inversion (using \code{\link{qgeneric}}).
##' 
##' \code{rsurvspline} generates random survival times by using
##' \code{qsurvspline} on a sample of uniform random numbers.  Due to the
##' numerical root-finding involved in \code{qsurvspline}, it is slow compared
##' to typical random number generation functions.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{flexsurvspline}}.
##' @references Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' @keywords distribution
##' @examples
##' 
##' ## reduces to the weibull
##' regscale <- 0.786; cf <- 1.82
##' a <- 1/regscale; b <- exp(cf)
##' dweibull(1, shape=a, scale=b)
##' dsurvspline(1, gamma=c(log(1 / b^a), a)) # should be the same
##' 
##' ## reduces to the log-normal
##' meanlog <- 1.52; sdlog <- 1.11
##' dlnorm(1, meanlog, sdlog) 
##' dsurvspline(1, gamma = c(-meanlog/sdlog, 1/sdlog), scale="normal")
##' # should be the same
##' @name Survspline
NULL

## Things to do that are common to d/p/q functions
## could be generalized to any function with vector of arguments 
## TODO put in h,H functions
## TODO more special value handling

dbase.survspline <- function(q, gamma, knots, scale, deriv=FALSE){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    if(!is.matrix(knots)) knots <- matrix(knots, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(q), lg)
    q <- rep(q, length=nret)

    gamma <- apply(gamma, 2, rep, length=nret)
    knots <- apply(knots, 2, rep, length=nret)
    
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=nret)
    if(!is.matrix(knots)) knots <- matrix(knots, nrow=nret)
    if (ncol(gamma) != ncol(knots)) {
        stop("length of gamma should equal number of knots")
    }
    scale <- match.arg(scale, c("hazard","odds","normal"))
    if (deriv){
        ret <- matrix(0, nrow=nret, ncol=ncol(gamma))
        ret[is.na(q),] <- NA
    } else {
        ret <- numeric(nret)
        ret[is.na(q)] <- NA
    }
    ind <- !is.na(q) & q > 0
    q <- q[ind]
    gamma <- gamma[ind,,drop=FALSE]
    knots <- knots[ind,,drop=FALSE]
    list(ret=ret, gamma=gamma, q=q, scale=scale, ind=ind, knots=knots)
}

dlink <- function(scale){
    switch(scale,
           hazard=function(x){exp(x - exp(x))},
           odds=function(x){exp(x) / (1 + exp(x))^2},
           normal=function(x){dnorm(x)}
           )
}

ldlink <- function(scale){
    switch(scale,
           hazard=function(x){x - exp(x)},
           odds=function(x){x - 2*log(1 + exp(x))},
           normal=function(x){dnorm(x, log=TRUE)}
           )
}

## probability density function.

##' @export
dsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0, log=FALSE){
    d <- dbase.survspline(q=x, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        eta <- rowSums(basis(knots, tsfn(q, timescale)) * gamma) + as.numeric(X %*% beta) + offset # log cumulative hazard/odds
        eeta <- exp(ldlink(scale)(eta))
        ret[ind][eeta==0] <- 0
        ret[ind][is.nan(eeta)] <- NaN
        ind2 <- !(eeta==0 | is.nan(eeta))
        q <- q[ind2]
        gamma <- gamma[ind2,,drop=FALSE]
        knots <- knots[ind2,,drop=FALSE] 
        eeta <- eeta[ind2]
        ind[ind] <- ind[ind] & ind2
        dsum <- rowSums(dbasis(knots, tsfn(q, timescale)) * gamma)  # ds/dx
        ret[ind] <- dtsfn(q,timescale) * dsum * eeta
        ## derivative of log cum haz cannot be negative by definition, but
        ## optimisation doesn't constrain gamma to respect this, so set
        ## likelihood to zero then (assuming at least one death)
        ret[ind][ret[ind]<=0] <- 0
    }
    if (log) {ret <- log(ret)}
    as.numeric(ret)
}

tsfn <- function(x, timescale){
    switch(timescale,
           log = log(x),
           identity = x)
}

dtsfn <- function(x, timescale){
    switch(timescale,
           log = 1/x,
           identity = 1)
}

Slink <- function(scale){
    switch(scale,
           hazard=function(x){exp(-exp(x))},
           odds=function(x){1 / (1 + exp(x))},
           normal=function(x){pnorm(-x)}
           )
}

## cumulative distribution function

##' @export
psurvspline <- function(q, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0, lower.tail=TRUE, log.p=FALSE){
    d <- dbase.survspline(q=q, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        eta <- rowSums(basis(knots, tsfn(q,timescale)) * gamma) + as.numeric(X %*% beta) + offset
        surv <- Slink(scale)(eta)
        ret[ind] <- as.numeric(1 - surv)
        ret[ind][q==0] <- 0
        ret[ind][q==Inf] <- 1
    }
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

##' @export
qsurvspline <- function(p, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0, lower.tail=TRUE, log.p=FALSE){
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    qgeneric(psurvspline, p=p, matargs=c("gamma","knots"), gamma=gamma, beta=beta, X=X, knots=knots, scale=scale, timescale=timescale, offset=offset)
}

##' @export
rsurvspline <- function(n, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0){
    if (length(n) > 1) n <- length(n)
    ret <- qsurvspline(p=runif(n), gamma=gamma, beta=beta, X=X, knots=knots, scale=scale, timescale=timescale, offset=offset)
    ret
}
  

Hlink <- function(scale){
    switch(scale,
           hazard=function(x){exp(x)}, # log cum haz, or log cum odds is a spline function of log time
           odds=function(x){log1p(exp(x))},
           normal=function(x){-pnorm(-x, log.p=TRUE)}
           )
}

## cumulative hazard function

##' @export
Hsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0){
    match.arg(scale, c("hazard","odds","normal"))
    d <- dbase.survspline(q=x, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        eta <- rowSums(basis(knots, tsfn(q,timescale)) * gamma) + as.numeric(X %*% beta) + offset
        ret[ind] <- as.numeric(Hlink(scale)(eta))
    }
    ret
}

hlink <- function(scale){
    switch(scale,
           hazard=function(x){exp(x)},
           odds=function(x){plogis(x)},
           normal=function(x){dnorm(-x)/pnorm(-x)}
           )
}

## hazard function

##' @export
hsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0){
## value for x=0?  currently zero, should it be limit as x reduces to 0? 
    match.arg(scale, c("hazard","odds","normal"))
    d <- dbase.survspline(q=x, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        eta <- rowSums(basis(knots, tsfn(q,timescale)) * gamma) + as.numeric(X %*% beta) + offset
        eeta <- hlink(scale)(eta)
        ret[ind] <- dtsfn(q, timescale) * rowSums(dbasis(knots, tsfn(q, timescale)) * gamma) * eeta
    }
    as.numeric(ret)
}

##' @export
rmst_survspline = function(t, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0, start=0){
  rmst_generic(psurvspline, t, start=start, gamma=gamma, knots=knots, beta=beta, X=X, scale=scale, timescale=timescale, offset=offset)
}

##' @export
mean_survspline = function(gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", timescale="log", offset=0){
  rmst_generic(psurvspline, Inf, start=0, gamma=gamma, knots=knots, beta=beta, X=X, scale=scale, timescale=timescale, offset=offset)
}

##' Natural cubic spline basis
##' 
##' Compute a basis for a natural cubic spline, using the parameterisation
##' described by Royston and Parmar (2002).  Used for flexible parametric
##' survival models.
##' 
##' The exact formula for the basis is given in \code{\link{flexsurvspline}}.
##' 
##' @aliases basis dbasis fss dfss
##' @param knots Vector of knot locations in increasing order, including the
##' boundary knots at the beginning and end.
##' @param x Vector of ordinates to compute the basis for.
##' @return A matrix with one row for each ordinate and one column for each
##' knot.
##' 
##' \code{basis} returns the basis, and \code{dbasis} returns its derivative
##' with respect to \code{x}.
##' 
##' \code{fss} and \code{dfss} are the same, but with the order of the
##' arguments swapped around for consistency with similar functions in other R
##' packages.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{flexsurvspline}}.
##' @references Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' @keywords models
##' @export
basis <- function(knots, x) {
    if (is.matrix(knots)) {
        basis_matrix(knots, x)
    } else {
        basis_vector(knots, x)
    }
}

##' @export
dbasis <- function(knots, x) {
    if (is.matrix(knots)) {
        dbasis_matrix(knots, x)
    } else {
        dbasis_vector(knots, x)
    }
}

##' @export
fss <- function(x, knots) { basis(knots, x) }

##' @export
dfss <- function(x, knots) { dbasis(knots, x) }


flexsurv.splineinits <- function(t=NULL, mf, mml, aux)
{
    Y <- check.flexsurv.response(model.extract(mf, "response"))

    ## Impute interval-censored data, where Cox doesn't work.
    intcens <- Y[,"status"]==3 & (Y[,"start"] < Y[,"stop"])
    Y[intcens,"status"] <- 1
    Y[intcens,"stop"] <- Y[intcens,"start"] + (Y[intcens,"stop"] - Y[intcens,"start"])/2  
    Y[,"status"] <- ifelse(Y[,"status"]==1, 1, 0)
    X <- mml[[1]][,-1,drop=FALSE]
   
    ## Use only uncensored observations, unless < 5 of those, in which
    ## case use all (coxph doesn't work if all censored). Impute as
    ## alternative?
    dead <- Y[,"status"]==1
    inc <- if (sum(dead) >= 5) dead else rep(TRUE, nrow(Y))
    inc <- inc & Y[,"start"] < Y[,"stop"]
   
    ## Estimate empirical cumulative hazard for each observation 
    formdat <- as.data.frame(cbind(Y, X, weights=model.extract(mf, "weights")))
    names(formdat)[1:ncol(Y)] <- colnames(Y)
    form <- c("Surv(start, stop, status) ~")
    if (ncol(X)>0){
        names(formdat)[ncol(Y) + 1:ncol(X)] <- paste0("X", 1:ncol(X))
        form <- paste(form, paste(paste0("X", 1:ncol(X)), collapse=" + "))
    }
    else form <- paste(form, "1")
    formdat <- formdat[inc,]
    cox <- coxph(as.formula(form), weights=weights, data=formdat)
    
    if (ncol(X)>0){
        newdata <- as.data.frame(matrix(rep(0, ncol(X)), nrow=1))
        names(newdata) <- paste0("X", 1:ncol(X))
        surv <- survfit(cox, newdata=newdata)
    } else surv <- survfit(cox)
    surv <- surv$surv[match(Y[inc,"stop"], surv$time)]
    surv <- surv^exp(cox$linear.predictors)
    if (aux$scale=="hazard")
        logH <- log(-log(surv))
    else if (aux$scale=="odds")
        logH <- log((1 - surv)/surv)
    else if (aux$scale=="normal")
        logH <- qnorm(1 - surv)
    x <- tsfn(Y[inc,"time"], aux$timescale)
    b <- if (!is.null(aux$knots)) basis(aux$knots, x) else cbind(1, x)
    
    ## Regress empirical logH on covariates and spline basis to obtain
    ## initial values.
    ## Constrain initial spline function to be increasing, using
    ## quadratic programming
    Xq <- cbind(b, X[inc,,drop=FALSE])
    kx <- x
    kr <- diff(range(x))
    kx[1] <- x[1] - 0.01*kr
    kx[length(kx)] <- x[length(kx)] + 0.01*kr
    db <- if (!is.null(aux$knots)) dbasis(aux$knots, kx) else cbind(0, kx)
    dXq <- cbind(db, matrix(0, nrow=nrow(db), ncol=ncol(X)))
    nints <- length(mml)-1
    for (i in seq_len(nints)){
        Xi <- mml[[i+1]][inc,-1,drop=FALSE]
        Xq <- cbind(Xq, Xi*b[,i+1])
        dXq <- cbind(dXq, Xi*db[,i+1])
    }
    y <- logH[is.finite(logH)]
    Xq <- Xq[is.finite(logH),,drop=FALSE]
    dXq <- dXq[is.finite(logH),,drop=FALSE]
    eps <- rep(1e-09, length(y)) # so spline is strictly increasing

    Dmat <- t(Xq) %*% Xq
    posdef <- all(eigen(Dmat)$values > 0) # TODO work out why non-positive definite matrices can occur here
    if (!posdef && is.null(aux$skip.pd.check))
        inits <- flexsurv.splineinits.cox(t=t, mf=mf, mml=mml, aux=aux)
    else if (posdef) inits <- solve.QP(Dmat=Dmat, dvec=t(t(y) %*% Xq), Amat=t(dXq), bvec=eps)$solution
    else inits <- rep(0, ncol(Xq)) # if Cox fails to get pos def
    inits
}

## Initialise coefficients on spline intercepts (standard hazard
## ratios in PH models) using standard Cox regression, in cases where
## the default inits routine above results in parameters that give a
## zero likelihood, e.g. GBSG2 example

flexsurv.splineinits.cox <- function(t=NULL, mf, mml, aux)
{
    aux$skip.pd.check <- TRUE # avoid infinite recursion
    inits <- flexsurv.splineinits(t=t, mf=mf, mml=mml, aux=aux)
    Y <- check.flexsurv.response(model.extract(mf, "response"))
    ## Impute interval-censored data, where Cox doesn't work.
    intcens <- Y[,"status"]==3 & (Y[,"start"] < Y[,"stop"])
    Y[intcens,"status"] <- 1
    Y[intcens,"stop"] <- Y[intcens,"start"] + (Y[intcens,"stop"] - Y[intcens,"start"])/2  
    Y[,"status"] <- ifelse(Y[,"status"]==1, 1, 0)
    inc <- Y[,"start"] < Y[,"stop"]   
    X <- mml[[1]][,-1,drop=FALSE]
    formdat <- as.data.frame(cbind(Y, X, weights=model.extract(mf, "weights")))
    names(formdat)[1:ncol(Y)] <- colnames(Y)
    form <- c("Surv(start, stop, status) ~")
    if (ncol(X)>0){
        names(formdat)[ncol(Y) + 1:ncol(X)] <- paste0("X", 1:ncol(X))
        form <- paste(form, paste(paste0("X", 1:ncol(X)), collapse=" + "))
        formdat <- formdat[inc,]
        cox <- coxph(as.formula(form), weights=weights, data=formdat)
        covinds <- length(aux$knots) + 1:ncol(X)
        inits[covinds] <- coef(cox)
    }
    inits
}



##' Flexible survival regression using the Royston/Parmar spline model.
##' 
##' Flexible parametric modelling of time-to-event data using the spline model
##' of Royston and Parmar (2002).
##' 
##' This function works as a wrapper around \code{\link{flexsurvreg}} by
##' dynamically constructing a custom distribution using
##' \code{\link{dsurvspline}}, \code{\link{psurvspline}} and
##' \code{\link{unroll.function}}.
##' 
##' In the spline-based survival model of Royston and Parmar (2002), a
##' transformation \eqn{g(S(t,z))} of the survival function is modelled as a
##' natural cubic spline function of log time \eqn{x = \log(t)}{x = log(t)}
##' plus linear effects of covariates \eqn{z}.
##' 
##' \deqn{g(S(t,z)) = s(x, \bm{\gamma}) + \bm{\beta}^T \mathbf{z}}{g(S(t,z)) =
##' s(x, gamma) + beta^T z}
##' 
##' The proportional hazards model (\code{scale="hazard"}) defines
##' \eqn{g(S(t,\mathbf{z})) = \log(-\log(S(t,\mathbf{z}))) =
##' \log(H(t,\mathbf{z}))}{g(S(t,z)) = log(-log(S(t,z))) = log(H(t,z))}, the
##' log cumulative hazard.
##' 
##' The proportional odds model (\code{scale="odds"}) defines
##' \eqn{g(S(t,\mathbf{z})) }{g(S(t,z)) = log(1/S(t,z) - 1)}\eqn{ =
##' \log(S(t,\mathbf{z})^{-1} - 1)}{g(S(t,z)) = log(1/S(t,z) - 1)}, the log
##' cumulative odds.
##' 
##' The probit model (\code{scale="normal"}) defines \eqn{g(S(t,\mathbf{z})) =
##' }{g(S(t,z)) = -InvPhi(S(t,z))}\eqn{ -\Phi^{-1}(S(t,\mathbf{z}))}{g(S(t,z))
##' = -InvPhi(S(t,z))}, where \eqn{\Phi^{-1}()}{InvPhi()} is the inverse normal
##' distribution function \code{\link{qnorm}}.
##' 
##' With no knots, the spline reduces to a linear function, and these models
##' are equivalent to Weibull, log-logistic and lognormal models respectively.
##' 
##' The spline coefficients \eqn{\gamma_j: j=1, 2 \ldots }{gamma_j: j=1, 2
##' \ldots}, which are called the "ancillary parameters" above, may also be
##' modelled as linear functions of covariates \eqn{\mathbf{z}}, as
##' 
##' \deqn{\gamma_j(\mathbf{z}) = \gamma_{j0} + \gamma_{j1}z_1 + \gamma_{j2}z_2
##' + ... }{gamma_j(z) = \gamma_{j0} + \gamma_{j1}z_1 + \gamma_{j2}z_2 + ... }
##' 
##' giving a model where the effects of covariates are arbitrarily flexible
##' functions of time: a non-proportional hazards or odds model.
##' 
##' Natural cubic splines are cubic splines constrained to be linear beyond
##' boundary knots \eqn{k_{min},k_{max}}{kmin,kmax}.  The spline function is
##' defined as
##' 
##' \deqn{s(x,\bm{\gamma}) = \gamma_0 + \gamma_1 x + \gamma_2 v_1(x) + \ldots +
##' }{s(x,gamma) = gamma0 + gamma1 x + gamma2 v1(x) + ... + gamma_{m+1}
##' vm(x)}\deqn{ \gamma_{m+1} v_m(x)}{s(x,gamma) = gamma0 + gamma1 x + gamma2
##' v1(x) + ... + gamma_{m+1} vm(x)}
##' 
##' where \eqn{v_j(x)}{vj(x)} is the \eqn{j}th basis function
##' 
##' \deqn{v_j(x) = (x - k_j)^3_+ - \lambda_j(x - k_{min})^3_+ - (1 - }{vj(x) =
##' (x - kj)^3_+ - \lambda_j(x - kmin)^3_+ - (1 -\lambda_j) (x -
##' kmax)^3_+}\deqn{ \lambda_j) (x - k_{max})^3_+}{vj(x) = (x - kj)^3_+ -
##' \lambda_j(x - kmin)^3_+ - (1 -\lambda_j) (x - kmax)^3_+}
##' 
##' \deqn{\lambda_j = \frac{k_{max} - k_j}{k_{max} - k_{min}}}{\lambda_j =
##' (kmax - kj) / (kmax - kmin)}
##' 
##' and \eqn{(x - a)_+ = max(0, x - a)}.
##' 
##' @param formula A formula expression in conventional R linear modelling
##' syntax. The response must be a survival object as returned by the
##' \code{\link{Surv}} function, and any covariates are given on the right-hand
##' side.  For example,
##' 
##' \code{Surv(time, dead) ~ age + sex}
##' 
##' specifies a model where the log cumulative hazard (by default, see
##' \code{scale}) is a linear function of the covariates \code{age} and
##' \code{sex}.
##' 
##' If there are no covariates, specify \code{1} on the right hand side, for
##' example \code{Surv(time, dead) ~ 1}.
##' 
##' Time-varying covariate effects can be specified using the method described
##' in \code{\link{flexsurvreg}} for placing covariates on ancillary
##' parameters.  The ancillary parameters here are named \code{gamma1},
##' \ldots{}, \code{gammar} where \code{r} is the number of knots \code{k} plus
##' one (the ``degrees of freedom'' as defined by Royston and Parmar).  So for
##' the default Weibull model, there is just one ancillary parameter
##' \code{gamma1}.
##' 
##' Therefore a model with one internal spline knot, where the equivalents of
##' the Weibull shape and scale parameters, but not the higher-order term
##' \code{gamma2}, vary with age and sex, can be specified as:
##' 
##' \code{Surv(time, dead) ~ age + sex + gamma1(age) + gamma1(sex)}
##' 
##' or alternatively (and more safely, see \code{flexsurvreg})
##' 
##' \code{Surv(time, dead) ~ age + sex, anc=list(gamma1=~age + sex)}
##' 
##' \code{Surv} objects of \code{type="right"},\code{"counting"},
##' \code{"interval1"} or \code{"interval2"} are supported, corresponding to
##' right-censored, left-truncated or interval-censored observations.
##' @param data A data frame in which to find variables supplied in
##' \code{formula}.  If not given, the variables should be in the working
##' environment.
##' @param weights Optional variable giving case weights.
##' @param bhazard Optional variable giving expected hazards for relative
##' survival models.
##' @param subset Vector of integers or logicals specifying the subset of the
##' observations to be used in the fit.
##' @param k Number of knots in the spline. The default \code{k=0} gives a
##' Weibull, log-logistic or lognormal model, if \code{"scale"} is
##' \code{"hazard"}, \code{"odds"} or \code{"normal"} respectively.  \code{k}
##' is equivalent to \code{df-1} in the notation of \code{stpm} for Stata.  The
##' knots are then chosen as equally-spaced quantiles of the log uncensored
##' survival times, for example, at the median with one knot, or at the 33\%
##' and 67\% quantiles of log time (or time, see \code{"timescale"}) with two
##' knots.  To override this default knot placement, specify \code{knots}
##' instead.
##' @param knots Locations of knots on the axis of log time (or time, see
##' \code{"timescale"}).  If not specified, knot locations are chosen as
##' described in \code{k} above.  Either \code{k} or \code{knots} must be
##' specified. If both are specified, \code{knots} overrides \code{k}.
##' @param bknots Locations of boundary knots, on the axis of log time (or
##' time, see \code{"timescale"}).  If not supplied, these are are chosen as
##' the minimum and maximum log death time.
##' @param scale If \code{"hazard"}, the log cumulative hazard is modelled as a
##' spline function.
##' 
##' If \code{"odds"}, the log cumulative odds is modelled as a spline function.
##' 
##' If \code{"normal"}, \eqn{-\Phi^{-1}(S(t))}{-InvPhi(S(t))} is modelled as a
##' spline function, where \eqn{\Phi^{-1}()}{InvPhi()} is the inverse normal
##' distribution function \code{\link{qnorm}}.
##' @param timescale If \code{"log"} (the default) the log cumulative hazard
##' (or alternative) is modelled as a spline function of log time.  If
##' \code{"identity"}, it is modelled as a spline function of time.
##' @param ...  Any other arguments to be passed to or through
##' \code{\link{flexsurvreg}}, for example, \code{anc}, \code{inits},
##' \code{fixedpars}, \code{weights}, \code{subset}, \code{na.action}, and any
##' options to control optimisation.  See \code{\link{flexsurvreg}}.
##' @return A list of class \code{"flexsurvreg"} with the same elements as
##' described in \code{\link{flexsurvreg}}, and including extra components
##' describing the spline model.  See in particular:
##' 
##' \item{k}{Number of knots.} \item{knots}{Location of knots on the log time
##' axis.} \item{scale}{The \code{scale} of the model, hazard, odds or normal.}
##' \item{res}{Matrix of maximum likelihood estimates and confidence limits.
##' Spline coefficients are labelled \code{"gamma..."}, and covariate effects
##' are labelled with the names of the covariates.
##' 
##' Coefficients \code{gamma1,gamma2,...} here are the equivalent of
##' \code{s0,s1,...} in Stata \code{streg}, and \code{gamma0} is the equivalent
##' of the \code{xb} constant term.  To reproduce results, use the
##' \code{noorthog} option in Stata, since no orthogonalisation is performed on
##' the spline basis here.
##' 
##' In the Weibull model, for example, \code{gamma0,gamma1} are
##' \code{-shape*log(scale), shape} respectively in \code{\link{dweibull}} or
##' \code{\link{flexsurvreg}} notation, or (\code{-Intercept/scale},
##' \code{1/scale}) in \code{\link[survival]{survreg}} notation.
##' 
##' In the log-logistic model with shape \code{a} and scale \code{b} (as in
##' \code{\link[eha]{dllogis}} from the \pkg{eha} package), \code{1/b^a} is
##' equivalent to \code{exp(gamma0)}, and \code{a} is equivalent to
##' \code{gamma1}.
##' 
##' In the log-normal model with log-scale mean \code{mu} and standard
##' deviation \code{sigma}, \code{-mu/sigma} is equivalent to \code{gamma0} and
##' \code{1/sigma} is equivalent to \code{gamma1}.  } \item{loglik}{The
##' maximised log-likelihood.  This will differ from Stata, where the sum of
##' the log uncensored survival times is added to the log-likelihood in
##' survival models, to remove dependency on the time scale.}
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{flexsurvreg}} for flexible survival modelling using
##' general parametric distributions.
##' 
##' \code{\link{plot.flexsurvreg}} and \code{\link{lines.flexsurvreg}} to plot
##' fitted survival, hazards and cumulative hazards from models fitted by
##' \code{\link{flexsurvspline}} and \code{\link{flexsurvreg}}.
##' @references Royston, P. and Parmar, M. (2002).  Flexible parametric
##' proportional-hazards and proportional-odds models for censored survival
##' data, with application to prognostic modelling and estimation of treatment
##' effects. Statistics in Medicine 21(1):2175-2197.
##' 
##' Jackson, C. (2016). flexsurv: A Platform for Parametric Survival Modeling
##' in R. Journal of Statistical Software, 70(8), 1-33.
##' doi:10.18637/jss.v070.i08
##' @keywords models survival
##' @examples
##' 
##' ## Best-fitting model to breast cancer data from Royston and Parmar (2002)
##' ## One internal knot (2 df) and cumulative odds scale
##' 
##' spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="odds")
##' 
##' ## Fitted survival
##' 
##' plot(spl, lwd=3, ci=FALSE)
##' 
##' ## Simple Weibull model fits much less well
##' 
##' splw <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=0, scale="hazard")
##' lines(splw, col="blue", ci=FALSE)
##' 
##' ## Alternative way of fitting the Weibull
##' 
##' \dontrun{
##' splw2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
##' }
##' 
##' @export
flexsurvspline <- function(formula, data, weights, bhazard, subset,
                           k=0, knots=NULL, bknots=NULL, scale="hazard", timescale="log", ...){
    ## Get response matrix from the formula.  Only need this to obtain
    ## default knots.  Largely copied from flexsurvreg - ideally
    ## should be in separate function, but can't make scoping work.

    call <- match.call()
    indx <- match(c("formula", "data", "weights", "bhazard", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    ## remove the predictors
    f2 <- as.formula(gsub("(.*~).*", "\\1 1", Reduce(paste, deparse(formula))))
    environment(f2) <- environment(formula)
    temp[["formula"]] <- f2
    if (missing(data)) temp[["data"]] <- environment(formula)
    if (missing(data)) data <- environment(formula) # TESTME to pass to flexsurvreg
    m <- eval(temp, parent.frame())
    Y <- check.flexsurv.response(model.extract(m, "response"))
    
    dtimes <- Y[,"stop"][Y[,"status"]==1]
    if (is.null(knots)) {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        if (is.null(k)) stop("either \"knots\" or \"k\" must be specified")
        if (!is.numeric(k)) stop("k must be numeric")
        if (!is.wholenumber(k) || (k<0)) stop("number of knots \"k\" must be a non-negative integer")
        knots <- quantile(tsfn(dtimes,timescale), seq(0, 1, length=k+2)[-c(1,k+2)])
    }
    else {
        if (!is.numeric(knots)) stop("\"knots\" must be a numeric vector")
        minlogtime <- min(tsfn(Y[,"stop"], timescale))
        if (any(knots <= minlogtime)) {
            badknots <- knots[knots < min(tsfn(Y[,"stop"],timescale))]
            plural <- if (length(badknots) > 1) "s" else ""
            stop(sprintf("knot%s %s less than or equal to minimum %stime", plural, paste(badknots,collapse=", "), (if (timescale=="log") "log " else "")))
        }
        maxlogtime <- max(tsfn(Y[,"stop"], timescale))
        if (any(knots >= maxlogtime)) {
            badknots <- knots[knots > maxlogtime]
            plural <- if (length(badknots) > 1) "s" else ""
            stop(sprintf("knot%s %s greater than or equal to maximum %stime", plural, paste(badknots,collapse=", "), (if (timescale=="log") "log " else "")))
        }
    }
    if (is.null(bknots)) {
        ## boundary knots chosen from min/max observed death times...
        if (length(dtimes)>0) { 
            bt <- dtimes
        } else {
            ## unless all observations censored, where censoring times used instead
            ## "time" used with right censoring 
            ## "time1", "time2" used with interval censoring
            bt <- c(Y[,"time1"], Y[,"time2"], Y[,"time"])
            bt <- bt[is.finite(bt)]
        }
        bknots <- c(min(tsfn(bt,timescale)), max(tsfn(bt,timescale)))
        if (bknots[1] == bknots[2])
            warning("minimum and maximum log death times are the same: knot and boundary knot locations should be supplied explicitly")
    } else 
        if (!is.numeric(bknots) || (length(bknots) !=2) ) stop("boundary knots should be a numeric vector of length 2")
    knots <- c(bknots[1], knots, bknots[2])
    
    nk <- length(knots)
    custom.fss <- list(
        name = "survspline", # unused, d,p functions passed through
        pars = c(paste0("gamma",0:(nk-1))),
        location = c("gamma0"),
        transforms = rep(c(identity), nk), inv.transforms=rep(c(identity), nk),
        inits = flexsurv.splineinits
        )
    aux <- list(knots=knots, scale=scale, timescale=timescale)
    dfn <- unroll.function(dsurvspline, gamma=0:(nk-1))
    pfn <- unroll.function(psurvspline, gamma=0:(nk-1))
    rfn <- unroll.function(rsurvspline, gamma=0:(nk-1))
    hfn <- unroll.function(hsurvspline, gamma=0:(nk-1))
    Hfn <- unroll.function(Hsurvspline, gamma=0:(nk-1))
    qfn <- unroll.function(qsurvspline, gamma=0:(nk-1))
    meanfn <- unroll.function(mean_survspline, gamma=0:(nk-1))
    rmstfn <- unroll.function(rmst_survspline, gamma=0:(nk-1))
    Ddfn <- if (scale=="normal") NULL else unroll.function(DLdsurvspline, gamma=0:(nk-1))
    DSfn <- if (scale=="normal") NULL else unroll.function(DLSsurvspline, gamma=0:(nk-1))
    args <- c(list(formula=formula, data=data, dist=custom.fss,
                   dfns=list(d=dfn,p=pfn,r=rfn,h=hfn,H=Hfn,rmst=rmstfn,mean=meanfn, q=qfn,
                             DLd=Ddfn,DLS=DSfn,deriv=!(scale=="normal")), aux=aux), list(...))

    ## Try an alternative initial value routine if the default one gives zero likelihood
    fpold <- args$fixedpars
    args$fixedpars <- TRUE
    if (is.infinite(do.call("flexsurvreg", args)$loglik)){
        args$dist$inits <- flexsurv.splineinits.cox
    }
    args$fixedpars <- fpold
    args$weights <- temp$weights
    args$bhazard <- temp$bhazard
    args$subset <- temp$subset

    ret <- do.call("flexsurvreg", args) # faff to make ... args work within functions
    ret <- c(ret, list(k=length(knots) - 2, knots=knots, scale=scale))
    ret$call <- call
    class(ret) <- "flexsurvreg"
    ret
}
