### FUNCTIONS FOR MULTI-STATE MODELLING 

# in the future:
#  \pkg{flexsurv} makes the \code{msfit} function generic, defines
#  the default method to be \code{\link[mstate]{msfit}} from \pkg{mstate},
#  and adds this new method for flexsurvreg objects.
#S3method(msfit, default)
#S3method(msfit, flexsurvreg)
#msfit <- function(object, ...) UseMethod("msfit")
#msfit.default <- function(object, ...) mstate::msfit(object, ...)

form.msm.newdata <- function(x, newdata=NULL, tvar="trans", trans){
    tr <- sort(unique(na.omit(as.vector(trans))))
    ntr <- length(tr)
    mfo <- model.frame(x)
    if (!(tvar %in% colnames(mfo))){
        if (missing(tvar))
            stop("\"tvar\" not supplied and variable \"", tvar, "\" not in model")
        else stop("\"variable \"", tvar, "\" not in model")
    }
    trobs <- unique(mfo[,tvar])
    if (!all(trobs %in% tr)) stop("\"tvar\" contains elements not in the transition indicator matrix \"trans\"")
    if(is.null(newdata)){
        newdata <- data.frame(trans=trobs); names(newdata) <- tvar
    } else {
        newdata <- as.data.frame(newdata)
        if (nrow(newdata)==1) newdata <- newdata[rep(1,ntr),,drop=FALSE]
        else if (nrow(newdata) != ntr) stop(sprintf("length of variables in \"newdata\" must be either 1 or number of transitions, %d", ntr))
        newdata[,tvar] <- trobs
    }
    newdata
}



##' Cumulative intensity function for parametric multi-state models
##' 
##' Cumulative transition-specific intensity/hazard functions for
##' fully-parametric multi-state or competing risks models, using a
##' piecewise-constant approximation that will allow prediction using the
##' functions in the \pkg{mstate} package.
##' 
##' 
##' @param object Output from \code{\link{flexsurvreg}} or
##' \code{\link{flexsurvspline}}, representing a fitted survival model object.
##' 
##' The model should have been fitted to data consisting of one row for each
##' observed transition and additional rows corresponding to censored times to
##' competing transitions.  This is the "long" format, or counting process
##' format, as explained in the \pkg{flexsurv} vignette.
##' 
##' The model should contain a categorical covariate indicating the transition.
##' In \code{flexsurv} this variable can have any name, indicated here by the
##' \code{tvar} argument.  In the Cox models demonstrated by \pkg{mstate} it is
##' usually included in model formulae as \code{strata(trans)}, but note that
##' the \code{strata} function does not do anything in \pkg{flexsurv}.  The
##' formula supplied to \code{\link{flexsurvreg}} should be precise about which
##' parameters are assumed to vary with the transition type.
##' 
##' Alternatively, if the parameters (including covariate effects) are assumed
##' to be different between different transitions, then a list of
##' transition-specific models can be formed.  This list has one component for
##' each permitted transition in the multi-state model.  This is more
##' computationally efficient, particularly for larger models and datasets.
##' See the example below, and the vignette.
##' @param t Vector of times.  These do not need to be the same as the observed
##' event times, and since the model is parametric, they can be outside the
##' range of the data.  A grid of more frequent times will provide a better
##' approximation to the cumulative hazard trajectory for prediction with
##' \code{\link[mstate]{probtrans}} or \code{\link[mstate]{mssample}}, at the
##' cost of greater computational expense.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  This must be specified if
##' there are other covariates. The variable names should be the same as those
##' in the fitted model formula.  There must be either one value per covariate
##' (the typical situation) or \eqn{n} values per covariate, a different one
##' for each of the \eqn{n} allowed transitions.
##' @param variance Calculate the variances and covariances of the transition
##' cumulative hazards (\code{TRUE} or \code{FALSE}).  This is based on
##' simulation from the normal asymptotic distribution of the estimates, which
##' is computationally-expensive.
##' @param tvar Name of the categorical variable in the model formula that
##' represents the transition number. The values of this variable should
##' correspond to elements of \code{trans}, conventionally a sequence of
##' integers starting from 1.  Not required if \code{x} is a list of
##' transition-specific models.
##' @param trans Matrix indicating allowed transitions in the multi-state
##' model, in the format understood by \pkg{mstate}: a matrix of integers whose
##' \eqn{r,s} entry is \eqn{i} if the \eqn{i}th transition type (reading across
##' rows) is \eqn{r,s}, and has \code{NA}s on the diagonal and where the
##' \eqn{r,s} transition is disallowed.
##' @param B Number of simulations from the normal asymptotic distribution used
##' to calculate variances.  Decrease for greater speed at the expense of
##' accuracy.
##' @return An object of class \code{"msfit"}, in the same form as the objects
##' used in the \pkg{mstate} package.  The \code{\link[mstate]{msfit}} method
##' from \pkg{mstate} returns the equivalent cumulative intensities for Cox
##' regression models fitted with \code{\link{coxph}}.
##' @author C. H. Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}
##' @seealso \pkg{flexsurv} provides alternative functions designed
##' specifically for predicting from parametric multi-state models without
##' calling \pkg{mstate}.  These include \code{\link{pmatrix.fs}} and
##' \code{\link{pmatrix.simfs}} for the transition probability matrix, and
##' \code{\link{totlos.fs}} and \code{\link{totlos.simfs}} for expected total
##' lengths of stay in states.  These are generally more efficient than going
##' via \pkg{mstate}.
##' @references Liesbeth C. de Wreede, Marta Fiocco, Hein Putter (2011).
##' \pkg{mstate}: An R Package for the Analysis of Competing Risks and
##' Multi-State Models. \emph{Journal of Statistical Software}, 38(7), 1-30.
##' \url{http://www.jstatsoft.org/v38/i07}
##' 
##' Mandel, M. (2013). "Simulation based confidence intervals for functions
##' with complicated derivatives." The American Statistician 67(2):76-81
##' @keywords models
##' @examples
##' 
##' ## 3 state illness-death model for bronchiolitis obliterans
##' ## Compare clock-reset / semi-Markov multi-state models
##' 
##' ## Simple exponential model (reduces to Markov)
##' 
##' bexp <- flexsurvreg(Surv(years, status) ~ trans,
##'                     data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' mexp <- msfit.flexsurvreg(bexp, t=seq(0,12,by=0.1),
##'                           trans=tmat, tvar="trans", variance=FALSE)
##' 
##' ## Cox semi-parametric model within each transition
##' 
##' bcox <- coxph(Surv(years, status) ~ strata(trans), data=bosms3)
##' 
##' if (require("mstate")){
##' 
##' mcox <- mstate::msfit(bcox, trans=tmat)
##' 
##' ## Flexible parametric spline-based model 
##' 
##' bspl <- flexsurvspline(Surv(years, status) ~ trans + gamma1(trans),
##'                        data=bosms3, k=3)
##' mspl <- msfit.flexsurvreg(bspl, t=seq(0,12,by=0.1),
##'                          trans=tmat, tvar="trans", variance=FALSE)
##' 
##' ## Compare fit: exponential model is OK but the spline is better
##' 
##' plot(mcox, lwd=1, xlim=c(0, 12), ylim=c(0,4))
##' cols <- c("black","red","green")
##' for (i in 1:3){
##'     lines(mexp$Haz$time[mexp$Haz$trans==i], mexp$Haz$Haz[mexp$Haz$trans==i],
##'              col=cols[i], lwd=2, lty=2)
##'     lines(mspl$Haz$time[mspl$Haz$trans==i], mspl$Haz$Haz[mspl$Haz$trans==i],
##'              col=cols[i], lwd=3)
##' }
##' legend("topright", lwd=c(1,2,3), lty=c(1,2,1),
##'    c("Cox", "Exponential", "Flexible parametric"), bty="n")
##' 
##' }
##' 
##' ## Fit a list of models, one for each transition
##' ## More computationally efficient, but only valid if parameters
##' ## are different between transitions.
##' 
##' \dontrun{
##' bexp.list <- vector(3, mode="list")
##' for (i in 1:3) { 
##'   bexp.list[[i]] <- flexsurvreg(Surv(years, status) ~ 1, subset=(trans==i),
##'                                 data=bosms3, dist="exp")
##' }
##' 
##' ## The list of models can be passed to this and other functions,
##' ## as if it were a single multi-state model. 
##' 
##' msfit.flexsurvreg(bexp.list, t=seq(0,12,by=0.1), trans=tmat)
##' }
##' 
##' @export
msfit.flexsurvreg <- function(object, t, newdata=NULL, variance=TRUE, tvar="trans",
                              trans, B=1000){
    tr <- sort(unique(na.omit(as.vector(trans))))
    ntr <- length(tr)
    if (is.flexsurvlist(object)) {
        Haz <- vector(ntr, mode="list")
        for (i in seq_len(ntr))
            Haz[[i]] <- summary(object[[i]], type="cumhaz", t=t, newdata=newdata, ci=FALSE)[[1]]
    } else {
        newdata <- form.msm.newdata(object, newdata=newdata, tvar=tvar, trans=trans)
        X <- form.model.matrix(object, newdata)
        Haz <- summary(object, type="cumhaz", t=t, X=X, ci=FALSE)
    }
    Haz <- do.call("rbind",Haz[seq_along(tr)])
    rownames(Haz) <- NULL
    Haz$trans <- rep(seq_along(tr), each=length(t))
    names(Haz)[names(Haz)=="est"] <- "Haz"
    res <- list(Haz=Haz, trans=trans)
    foundse <- if (is.flexsurvlist(object)) all(!is.na(sapply(object, function(x)x$cov[[1]]))) else !is.na(object$cov[1])
    if (variance && foundse){
        boot <- array(dim=c(B, length(t), ntr))
        for (i in seq_along(tr))
            boot[,,i] <-
                if (is.flexsurvlist(object))
                    normbootfn.flexsurvreg(object[[i]], t=t, start=0, newdata=newdata, B=B,
                                           fn=summary.fns(object[[i]],"cumhaz"))
                else
                    normbootfn.flexsurvreg(object, t=t, start=0, X=X[i,,drop=FALSE], B=B,
                                           fn=summary.fns(object,"cumhaz"))
        ntr2 <- 0.5*ntr*(ntr+1)
        nt <- length(t)
        mat <- matrix(nrow=ntr, ncol=ntr)
        trans1 <- rep(t(row(mat))[!t(lower.tri(mat))], each=nt)
        trans2 <- rep(t(col(mat))[!t(lower.tri(mat))], each=nt)
        res$varHaz <- data.frame(time=rep(t, ntr2),  varHaz=numeric(ntr2*nt),
                                 trans1=trans1, trans2=trans2)
        for (i in 1:ntr){
            for (j in i:ntr){
                res$varHaz$varHaz[trans1==i & trans2==j] <-
                    mapply(cov, split(t(boot[,,i]),seq_along(t)), split(t(boot[,,j]),seq_along(t)))
            }            
        }        
    } 
    class(res) <- "msfit"
    res
}

## Matrix with ntrans rows, npars columns giving transition-specific
## baseline parameters at given covariate values



##' Transition-specific parameters in a flexible parametric multi-state model
##' 
##' Matrix of maximum likelihood estimates of transition-specific parameters in
##' a flexible parametric multi-state model, at given covariate values.
##' 
##' 
##' @param x A multi-state model fitted with \code{\link{flexsurvreg}}.  See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.
##' 
##' \code{x} can also be a list of \code{\link{flexsurvreg}} models, with one
##' component for each permitted transition in the multi-state model, as
##' illustrated in \code{\link{msfit.flexsurvreg}}.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @return A matrix with one row for each permitted transition, and one column
##' for each parameter of the parametric distribution that generates each event
##' in the multi-state model.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @keywords models,survival
##' @export
pars.fmsm <- function(x, trans, newdata=NULL, tvar="trans")
{
    if (is.flexsurvlist(x)){
        ntr <- length(x) # number of allowed transitions
        if (ntr != length(na.omit(as.vector(trans)))) stop(sprintf("x is a list of %s flexsurvreg objects, but trans indicates %s transitions", ntr, length(na.omit(as.vector(trans)))))
        basepar <- matrix(nrow=ntr, ncol=length(x[[1]]$dlist$pars), dimnames=list(NULL,x[[1]]$dlist$pars))
        newdata <- as.data.frame(newdata)
        for (i in 1:ntr){
            if (x[[i]]$ncovs==0)
                X <- matrix(0)
            else {
                if(nrow(newdata) == 1L) {
                    X <- form.model.matrix(x[[i]], as.data.frame(newdata))
                } else if(nrow(newdata) == ntr){
                    X <- form.model.matrix(x[[i]], as.data.frame(newdata[i, ,drop = FALSE]))
                } else stop(sprintf("`newdata` has %s rows. It must either have one row, or one row for each of the %s allowed transitions",nrow(newdata),ntr))
              }
            beta <- if (x[[i]]$ncovs==0) 0 else x[[i]]$res.t[x[[i]]$covpars,"est"]
            basepar[i,] <- add.covs(x[[i]], x[[i]]$res.t[x[[i]]$dlist$pars,"est"], beta, X, transform=FALSE)
        }
    } else if (inherits(x, "flexsurvreg")) {
        newdata <- form.msm.newdata(x, newdata=newdata, tvar=tvar, trans=trans)
        X <- form.model.matrix(x, newdata)   
        basepar <- add.covs(x, pars=x$res.t[x$dlist$pars,"est"], beta=x$res.t[x$covpars,"est"], X=X)
    } else
        stop("expected x to be a flexsurvreg object or list of flexsurvreg objects") 
    basepar
}

## TODO allow time dependent covs when computing the hazard.



##' Transition probability matrix from a fully-parametric, time-inhomogeneous
##' Markov multi-state model
##' 
##' The transition probability matrix for time-inhomogeneous Markov multi-state
##' models fitted to time-to-event data with \code{\link{flexsurvreg}}.  This
##' has \eqn{r,s} entry giving the probability that an individual is in state
##' \eqn{s} at time \eqn{t}, given they are in state \eqn{r} at time \eqn{0}.
##' 
##' This is computed by solving the Kolmogorov forward differential equation
##' numerically, using the methods in the \code{\link{deSolve}} package.  The
##' equation is
##' 
##' \deqn{\frac{dP(t)}{dt} = P(t) Q(t)}
##' 
##' where \eqn{P(t)} is the transition probability matrix for time \eqn{t}, and
##' \eqn{Q(t)} is the transition hazard or intensity as a function of \eqn{t}.
##' The initial condition is \eqn{P(0) = I}.
##' 
##' Note that the package \pkg{msm} has a similar method \code{pmatrix.msm}.
##' \code{pmatrix.fs} should give the same results as \code{pmatrix.msm} when
##' both of these conditions hold:
##' 
##' \itemize{ \item the time-to-event distribution is exponential for all
##' transitions, thus the \code{flexsurvreg} model was fitted with
##' \code{dist="exp"} and the model is time-homogeneous.  \item the \pkg{msm}
##' model was fitted with \code{exacttimes=TRUE}, thus all the event times are
##' known, and there are no time-dependent covariates.  }
##' 
##' \pkg{msm} only allows exponential or piecewise-exponential time-to-event
##' distributions, while \pkg{flexsurvreg} allows more flexible models.
##' \pkg{msm} however was designed in particular for panel data, where the
##' process is observed only at arbitrary times, thus the times of transition
##' are unknown, which makes flexible models difficult.
##' 
##' This function is only valid for Markov ("clock-forward") multi-state
##' models, though no warning or error is currently given if the model is not
##' Markov.  See \code{\link{pmatrix.simfs}} for the equivalent for semi-Markov
##' ("clock-reset") models.
##' 
##' @param x A model fitted with \code{\link{flexsurvreg}}.  See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.  Additionally, this must be a Markov / clock-forward model, but can
##' be time-inhomogeneous.  See the package vignette for further explanation.
##' 
##' \code{x} can also be a list of models, with one component for each
##' permitted transition, as illustrated in \code{\link{msfit.flexsurvreg}}.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param t Time or vector of times to predict state occupancy probabilities
##' for.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param ci Return a confidence interval calculated by simulating from the
##' asymptotic normal distribution of the maximum likelihood estimates.  Turned
##' off by default, since this is computationally intensive.  If turned on,
##' users should increase \code{B} until the results reach the desired
##' precision.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @param sing.inf If there is a singularity in the observed hazard, for
##' example a Weibull distribution with \code{shape < 1} has infinite hazard at
##' \code{t=0}, then as a workaround, the hazard is assumed to be a large
##' finite number, \code{sing.inf}, at this time.  The results should not be
##' sensitive to the exact value assumed, but users should make sure by
##' adjusting this parameter in these cases.
##' @param B Number of simulations from the normal asymptotic distribution used
##' to calculate variances.  Decrease for greater speed at the expense of
##' accuracy.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @param ... Arguments passed to \code{\link{ode}} in \pkg{deSolve}.
##' @return The transition probability matrix, if \code{t} is of length 1, or a
##' list of matrices if \code{t} is longer.
##' 
##' If \code{ci=TRUE}, each element has attributes \code{"lower"} and
##' \code{"upper"} giving matrices of the corresponding confidence limits.
##' These are formatted for printing but may be extracted using \code{attr()}.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @seealso \code{\link{pmatrix.simfs}}, \code{\link{totlos.fs}},
##' \code{\link{msfit.flexsurvreg}}.
##' @keywords models,survival
##' @examples
##' 
##' # BOS example in vignette, and in msfit.flexsurvreg
##' bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans,
##'                     data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' # more likely to be dead (state 3) as time moves on, or if start with
##' # BOS (state 2)
##' pmatrix.fs(bexp, t=c(5,10), trans=tmat)
##' @export
pmatrix.fs <- function(x, trans, t=1, newdata=NULL, ci=FALSE,
                       tvar="trans", sing.inf=1e+10, B=1000, cl=0.95, ...){
    ntr <- sum(!is.na(trans))
    n <- nrow(trans)
    dp <- function(t, y, parms, ...){
        if (is.flexsurvlist(x)) x <- x[[1]] 
        P <- matrix(y, nrow=n, ncol=n)
        haz <- numeric(n)
        for (i in 1:ntr){
            hcall <- list(x=t)
            for (j in seq(along=x$dlist$pars))
                hcall[[x$dlist$pars[j]]] <- parms$par[i,j]
            haz[i] <- do.call(x$dfns$h, hcall)
        }
        Q <- haz[trans]
        Q[is.na(Q)] <- 0
        Q[is.infinite(Q) & Q>0] <- sing.inf
        Q <- matrix(Q, nrow=n, ncol=n)
        diag(Q) <- -rowSums(Q)
        list(P %*% Q)
    }
    nt <- length(t)
    if (nt<1) stop("number of times should be at least one")
    basepar <- pars.fmsm(x=x, trans=trans, newdata=newdata, tvar=tvar)
    res <- ode(y=diag(n), times=c(0,t), func=dp, parms=list(par=basepar), ...)[-1,-1]
    res <- lapply(split(res,1:nt), function(x)matrix(x,nrow=n))
    names(res) <- t
    if (ci){
        resci <- bootci.fmsm(x, B, nt*n*n, match.call(), cl)
        resl <- lapply(split(resci[1,],rep(1:nt, each=n*n)), function(x)matrix(x,nrow=n))
        resu <- lapply(split(resci[2,],rep(1:nt, each=n*n)), function(x)matrix(x,nrow=n))
        names(resl) <- names(resu) <- t
        for (i in 1:nt){
            attr(res[[i]], "lower") <- resl[[i]]
            attr(res[[i]], "upper") <- resu[[i]]
            class(res[[i]]) <- "fs.msm.est"
        }
    }
    if(nt==1) res[[1]] else res
}

## Obtains matrix T(t) of expected times spent in state (col) starting
## from state (row) up to time t.
## Solves the second order linear ODE system T''(t) = P(t) Q(t)
## Express as dT/dt = P(t), dP/dt = P(t)Q(t), solve for both P and T at once



##' Total length of stay in particular states for a fully-parametric,
##' time-inhomogeneous Markov multi-state model
##' 
##' The matrix whose \eqn{r,s} entry is the expected amount of time spent in
##' state \eqn{s} for a time-inhomogeneous, continuous-time Markov multi-state
##' process that starts in state \eqn{r}, up to a maximum time \eqn{t}. This is
##' defined as the integral of the corresponding transition probability up to
##' that time.
##' 
##' This is computed by solving a second order extension of the Kolmogorov
##' forward differential equation numerically, using the methods in the
##' \code{\link{deSolve}} package.  The equation is expressed as a linear
##' system
##' 
##' \deqn{\frac{dT(t)}{dt} = P(t)} \deqn{\frac{dP(t)}{dt} = P(t) Q(t)}
##' 
##' and solved for \eqn{T(t)} and \eqn{P(t)} simultaneously, where \eqn{T(t)}
##' is the matrix of total lengths of stay, \eqn{P(t)} is the transition
##' probability matrix for time \eqn{t}, and \eqn{Q(t)} is the transition
##' hazard or intensity as a function of \eqn{t}.  The initial conditions are
##' \eqn{T(0) = 0} and \eqn{P(0) = I}.
##' 
##' Note that the package \pkg{msm} has a similar method \code{totlos.msm}.
##' \code{totlos.fs} should give the same results as \code{totlos.msm} when
##' both of these conditions hold:
##' 
##' \itemize{ \item the time-to-event distribution is exponential for all
##' transitions, thus the \code{flexsurvreg} model was fitted with
##' \code{dist="exp"}, and is time-homogeneous.  \item the \pkg{msm} model was
##' fitted with \code{exacttimes=TRUE}, thus all the event times are known, and
##' there are no time-dependent covariates.  }
##' 
##' \pkg{msm} only allows exponential or piecewise-exponential time-to-event
##' distributions, while \pkg{flexsurvreg} allows more flexible models.
##' \pkg{msm} however was designed in particular for panel data, where the
##' process is observed only at arbitrary times, thus the times of transition
##' are unknown, which makes flexible models difficult.
##' 
##' This function is only valid for Markov ("clock-forward") multi-state
##' models, though no warning or error is currently given if the model is not
##' Markov.  See \code{\link{totlos.simfs}} for the equivalent for semi-Markov
##' ("clock-reset") models.
##' 
##' @param x A model fitted with \code{\link{flexsurvreg}}.  See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.  Additionally, this must be a Markov / clock-forward model, but can
##' be time-inhomogeneous.  See the package vignette for further explanation.
##' 
##' \code{x} can also be a list of models, with one component for each
##' permitted transition, as illustrated in \code{\link{msfit.flexsurvreg}}.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param t Time or vector of times to predict up to.  Must be finite.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param ci Return a confidence interval calculated by simulating from the
##' asymptotic normal distribution of the maximum likelihood estimates.  Turned
##' off by default, since this is computationally intensive.  If turned on,
##' users should increase \code{B} until the results reach the desired
##' precision.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @param sing.inf If there is a singularity in the observed hazard, for
##' example a Weibull distribution with \code{shape < 1} has infinite hazard at
##' \code{t=0}, then as a workaround, the hazard is assumed to be a large
##' finite number, \code{sing.inf}, at this time.  The results should not be
##' sensitive to the exact value assumed, but users should make sure by
##' adjusting this parameter in these cases.
##' @param B Number of simulations from the normal asymptotic distribution used
##' to calculate variances.  Decrease for greater speed at the expense of
##' accuracy.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @param ... Arguments passed to \code{\link{ode}} in \pkg{deSolve}.
##' @return The matrix of lengths of stay \eqn{T(t)}, if \code{t} is of length
##' 1, or a list of matrices if \code{t} is longer.
##' 
##' If \code{ci=TRUE}, each element has attributes \code{"lower"} and
##' \code{"upper"} giving matrices of the corresponding confidence limits.
##' These are formatted for printing but may be extracted using \code{attr()}.
##' 
##' The result also has an attribute \code{P} giving the transition probability
##' matrices, since these are unavoidably computed as a side effect.  These are
##' suppressed for printing, but can be extracted with \code{attr(...,"P")}.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @seealso \code{\link{totlos.simfs}}, \code{\link{pmatrix.fs}},
##' \code{\link{msfit.flexsurvreg}}.
##' @keywords models,survival
##' @examples
##' 
##' # BOS example in vignette, and in msfit.flexsurvreg
##' bexp <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans,
##'                     data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' 
##' # predict 4 years spent without BOS, 3 years with BOS, before death
##' # As t increases, this should converge
##' 
##' totlos.fs(bexp, t=10, trans=tmat)
##' totlos.fs(bexp, t=1000, trans=tmat)
##' totlos.fs(bexp, t=c(5,10), trans=tmat)
##' 
##' # Answers should match results in help(totlos.simfs) up to Monte Carlo
##' # error there / ODE solving precision here, since with an exponential
##' # distribution, the "semi-Markov" model there is the same as the Markov
##' # model here
##' @export
totlos.fs <- function(x, trans, t=1, newdata=NULL, ci=FALSE,
                       tvar="trans", sing.inf=1e+10, B=1000, cl=0.95, ...){
    ntr <- sum(!is.na(trans))
    n <- nrow(trans)
    nsq <- n*n
    dp <- function(t, y, parms, ...){
        if (is.flexsurvlist(x)) x <- x[[1]] 
        P <- matrix(y[nsq + 1:nsq], nrow=n, ncol=n)
        haz <- numeric(n)
        for (i in 1:ntr){
            hcall <- list(x=t)
            for (j in seq(along=x$dlist$pars))
                hcall[[x$dlist$pars[j]]] <- parms$par[i,j]
            haz[i] <- do.call(x$dfns$h, hcall)
        }
        Q <- haz[trans]
        Q[is.na(Q)] <- 0
        Q[is.infinite(Q) & Q>0] <- sing.inf
        Q <- matrix(Q, nrow=n, ncol=n)
        diag(Q) <- -rowSums(Q)
        list(cbind(P, P %*% Q))
    }
    nt <- length(t)
    if (nt<1) stop("number of times should be at least one")
    basepar <- pars.fmsm(x=x, trans=trans, newdata=newdata, tvar=tvar)
    init <- cbind(matrix(0, nrow=n, ncol=n), diag(n))
    res <- ode(y=init, times=c(0,t), func=dp, parms=list(par=basepar), ...)[-1,-1]
    res.t <- lapply(split(res,1:nt), function(x)matrix(x[1:nsq],nrow=n))
    res.p <- lapply(split(res,1:nt), function(x)matrix(x[nsq + 1:nsq],nrow=n))
    names(res.t) <- names(res.p) <- t
    if (ci){
        resci <- bootci.fmsm(x, B, nt*2*n*n, match.call(), cl)
        tind <- rep(rep(1:nt,each=n*n), 2)
        res.tl <- lapply(split(resci[1,],tind), function(x)matrix(x[1:nsq],nrow=n))
        res.tu <- lapply(split(resci[2,],tind), function(x)matrix(x[1:nsq],nrow=n))
        res.pl <- lapply(split(resci[1,],tind), function(x)matrix(x[nsq + 1:nsq],nrow=n))
        res.pu <- lapply(split(resci[2,],tind), function(x)matrix(x[nsq + 1:nsq],nrow=n))
        names(res.tl) <- names(res.tu) <- names(res.pl) <- names(res.pu) <- t
        for (i in 1:nt){
            attr(res.t[[i]], "lower") <- res.tl[[i]]
            attr(res.t[[i]], "upper") <- res.tu[[i]]
            class(res.t[[i]]) <- "fs.msm.est"
            attr(res.p[[i]], "lower") <- res.pl[[i]]
            attr(res.p[[i]], "upper") <- res.pu[[i]]
            class(res.p[[i]]) <- "fs.msm.est"
        }
    }
    if(nt==1) {res.t <- res.t[[1]]; res.p <- res.p[[1]]}
    attr(res.t, "P") <- res.p
    class(res.t) <- "totlos.fs"
    res.t
}

##' @export
print.totlos.fs <- function(x, ...){attr(x, "P") <- NULL; print(unclass(x),...)}

# TODO make pmatrix generic
# pmatrix.flexsurvreg <- pmatrix.fs

format.ci <- function(x, l, u, digits=NULL, ...)
{
    if (is.null(digits)) digits <- 4
    ## note format() aligns nicely on point, unlike formatC
    est <- format(x, digits=digits, ...)
    if (!is.null(l)) {
        low <- format(l, digits=digits, ...)
        upp <- format(u, digits=digits, ...)
        res <- paste(est, " (", low, ",", upp, ")", sep="")
        res[x==0] <- 0
    }
    else res <- est
    dim(res) <- dim(x)
    dimnames(res) <- dimnames(x)
    names(res) <- names(x)
    res
}

print.ci <- function(x, l, u, digits=NULL){
    res <- format.ci(x, l, u, digits)
    print(res, quote=FALSE)
}

print.fs.msm.est <- function(x, digits=NULL, ...)
{
    if (!is.null(attr(x, "lower")))
        print.ci(x, attr(x, "lower"), attr(x, "upper"), digits=digits)
    else print(unclass(x))
}

absorbing <- function(trans){
    which(apply(trans, 1, function(x)all(is.na(x))))
}

transient <- function(trans){
    which(apply(trans, 1, function(x)any(!is.na(x))))
}

is.flexsurvlist <- function(x){
    is.list(x) &&
        (length(x) > 0) && 
            inherits(x[[1]], "flexsurvreg") && 
                all(sapply(x, inherits, "flexsurvreg"))
}

## Handle predictable time-dependent covariates in simulating from
## semi-Markov models.  Assume the covariate changes at same rate as
## time (e.g. age), but the covariate values used in simulation only
## change when the clock resets, at each change of state.

form.basepars.tcovs <- function(x, transi, # index of allowed transition
                                newdata, tcovs,
                                t # time increment
                                ){
    if (is.flexsurvlist(x)){
        x <- x[[transi]]
        dat <- as.list(newdata)
    } else if (inherits(x, "flexsurvreg")) {
        dat <- as.list(newdata[transi,,drop=FALSE])
    }
    for (i in tcovs) { dat[[i]] <- dat[[i]] + t}
    dat <- as.data.frame(dat)
    X <- form.model.matrix(x, dat)
    beta <- if (x$ncovs==0) 0 else x$res.t[x$covpars,"est"]
    basepars.mat <- add.covs(x, x$res.t[x$dlist$pars,"est"], beta, X, transform=FALSE)
    as.list(as.data.frame(basepars.mat))
    ## for distribution with npars parameters, whose values are required at nt different times
    ## returns (as list) data frame with nt rows, npars cols
}

## TODO Unclear how to check for semi Markov vs nonhomogenous Markov
## model.  attr(model.response(model.frame(x)), "type") will be
## "counting" for a nonhomogeneous model, but also if there are
## time-dependent covariates



##' Simulate paths through a fully parametric semi-Markov multi-state model
##' 
##' Simulate changes of state and transition times from a semi-Markov
##' multi-state model fitted using \code{\link{flexsurvreg}}.
##' 
##' \code{sim.fmsm} relies on the presence of a function to sample random
##' numbers from the parametric survival distribution used in the fitted model
##' \code{x}, for example \code{\link{rweibull}} for Weibull models. If
##' \code{x} was fitted using a custom distribution, called \code{dist} say,
##' then there must be a function called (something like) \code{rdist} either
##' in the working environment, or supplied through the \code{dfns} argument to
##' \code{\link{flexsurvreg}}.  This must be in the same format as standard R
##' functions such as \code{\link{rweibull}}, with first argument \code{n}, and
##' remaining arguments giving the parameters of the distribution.  It must be
##' vectorised with respect to the parameter arguments.
##' 
##' This function is only valid for semi-Markov ("clock-reset") models, though
##' no warning or error is currently given if the model is not of this type. An
##' equivalent for time-inhomogeneous Markov ("clock-forward") models has
##' currently not been implemented.
##' 
##' Note the random sampling method for \code{flexsurvspline} models is
##' currently very inefficient, so that looping over the \code{M} individuals
##' will be very slow.
##' 
##' @param x A model fitted with \code{\link{flexsurvreg}}. See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.
##' 
##' Alternatively \code{x} can be a list of fitted \code{\link{flexsurvreg}}
##' model objects.  The \code{i}th element of this list is the model
##' corresponding to the \code{i}th transition in \code{trans}.  This is a more
##' efficient way to fit a multi-state model, but only valid if the parameters
##' are different between different transitions.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param t Time, or vector of times for each of the \code{M} individuals, to
##' simulate trajectories until.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param start Starting state, or vector of starting states for each of the
##' \code{M} individuals.
##' @param M Number of individual trajectories to simulate.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @param tcovs Names of "predictable" time-dependent covariates in
##' \code{newdata}, i.e. those whose values change at the same rate as time.
##' Age is a typical example.  During simulation, their values will be updated
##' after each transition time, by adding the current time to the value
##' supplied in \code{newdata}.  This assumes the covariate is measured in the
##' same unit as time. \code{tcovs} is supplied as a character vector.
##' @param debug Print intermediate outputs: for development use.
##' @return A list of two matrices named \code{st} and \code{t}.  The rows of
##' each matrix represent simulated individuals.  The columns of \code{t}
##' contain the times when the individual changes state, to the corresponding
##' states in \code{st}.
##' 
##' The first columns will always contain the starting states and the starting
##' times. The last column of \code{t} represents either the time when the
##' individual moves to an absorbing state, or right-censoring in a transient
##' state at the time given in the \code{t} argument to \code{sim.fmsm}.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @seealso \code{\link{pmatrix.simfs}},\code{\link{totlos.simfs}}
##' @keywords models,survival
##' @examples
##' 
##' bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' sim.fmsm(bexp, M=10, t=5, trans=tmat)
##' @export
sim.fmsm <- function(x, trans, t, newdata=NULL, start=1, M=10, tvar="trans", tcovs=NULL, debug=FALSE){
    if (length(t)==1) t <- rep(t, M)
    else if (length(t)!=M) stop("length of t should be 1 or M=",M)
    if (length(start)==1) start <- rep(start, M)
    else if (length(start)!=M) stop("length of start should be 1 or M=",M)

    basepars.mat <- pars.fmsm(x=x, trans=trans, newdata=newdata, tvar=tvar)
    xbase <- if (is.flexsurvlist(x)) x[[1]] else x

    nst <- nrow(trans)
    ## TODO only need a max time if model is transient, else if absorbing, can allocate these up front
    res.st <- cur.st <- start
    res.t <- cur.t <- rep(0, M)
    todo <- seq_len(M)
    while (any(todo)){
        if (debug) { cat("cur.t\n"); cat(cur.t); cat("\n") }
        if (debug) { cat("cur.st\n"); cat(cur.st); cat("\n") }
        if (debug) { cat("TODO\n"); cat(todo); cat("\n") }
        cur.st.out <- cur.st[todo]
        cur.t.out <- cur.t[todo]
        done <- numeric()
        for (i in unique(cur.st[todo])){            
            if (i %in% transient(trans)) {
                ## simulate next time and states for people whose current state is i
                transi <- na.omit(trans[i,])
                ni <- sum(cur.st[todo]==i)
                t.trans1 <- matrix(0, nrow=ni, ncol=length(transi))
                ## simulate times to all potential destination states
                for (j in seq_along(transi)) {         
                    if (length(tcovs)>0){
                        basepars <- form.basepars.tcovs(x, transi[j], newdata, tcovs, cur.t.out)
                    } else 
                        basepars <- as.list(as.data.frame(basepars.mat)[transi[j],])
                    fncall <- c(list(n=ni), basepars, xbase$aux)
                    if (is.null(xbase$dfns$r)) stop("No random sampling function found for this model")
                    t.trans1[,j] <- do.call(xbase$dfns$r, fncall)
                }
                if (debug) { print(t(t.trans1)) }
                ## simulated next state is the one with minimum simulated time
                mc <- max.col(-t.trans1)
                if (debug) { cat("mc\n"); cat(mc); cat("\n") }
                if (debug) { cat("transi\n"); cat(transi); cat("\n") }
                if (debug) { cat("trans[i,]\n"); cat(trans[i,]); cat("\n") }
                next.state <- match(transi[mc], trans[i,])
                if (debug) { cat("next.state\n"); cat(next.state); cat("\n") }
                next.time <- t.trans1[cbind(seq_along(next.state), mc)]
                if (debug) { cat("next.time\n"); cat(next.time); cat("\n") }
                inds <- which(cur.st[todo]==i)
                cur.t.out[inds] <- cur.t.out[inds] + next.time
                ## if final simulated state is greater than target time, censor at target time
                cens <- cur.t.out[inds] > t[inds]
                cur.t.out[inds][cens] <- t[inds][cens]               
                cur.st.out[!cens] <- next.state[!cens]
                done <- todo[inds][cens]
            }
        }
        cur.st[todo] <- cur.st.out
        cur.t[todo] <- cur.t.out
        res.st <- cbind(res.st, cur.st)
        res.t <- cbind(res.t, cur.t)
        done <- union(done, which(cur.st %in% absorbing(trans)))
        todo <- setdiff(todo, done)
        if (debug) { cat("\n") }
    }
    list(st=unname(res.st), t=unname(res.t))
}    

### Generic CIs for multi-state output functions.
### Replace the parameters in the fitted model object with a draw from the MVN distribution of the MLEs
### Then calls the output function again on this tweaked model object, with ci=FALSE. Repeat B times.

bootci.fmsm <- function(x, B, ncols, fncall, cl){
    res.rep <- array(NA_real_, dim=c(B, ncols))
    fncall$ci <- FALSE
    if (is.flexsurvlist(x)){
        sim <- vector("list", length(x))
        for (j in seq_along(x)){
            sim[[j]] <- normboot.flexsurvreg(x=x[[j]], B=B, raw=TRUE, transform=TRUE)
        }
        for (i in 1:B){
            x.rep <- x
            for (j in seq_along(x))
                x.rep[[j]]$res.t[,"est"] <- sim[[j]][i,]
            fncall$x <- x.rep
            resi <- eval(fncall)
            ## "P" attribute needed for totlos.fs. need to extend if ever return any other results as attributes.
            res.rep[i,] <- c(unlist(resi), unlist(attr(resi, "P")))
        }
    } else {
        sim <- normboot.flexsurvreg(x=x, B=B, raw=TRUE, transform=TRUE)
        for (i in 1:B){
            x.rep <- x
            x.rep$res.t[,"est"] <- sim[i,]
            fncall$x <- x.rep
            resi <- eval(fncall)
            res.rep[i,] <- c(unlist(resi), unlist(attr(resi, "P")))
        }
    }
    resci <- apply(res.rep, 2, quantile, c((1-cl)/2, 1 - (1-cl)/2))
    resci
}



##' Transition probability matrix from a fully-parametric, semi-Markov
##' multi-state model
##' 
##' The transition probability matrix for semi-Markov multi-state models fitted
##' to time-to-event data with \code{\link{flexsurvreg}}.  This has \eqn{r,s}
##' entry giving the probability that an individual is in state \eqn{s} at time
##' \eqn{t}, given they are in state \eqn{r} at time \eqn{0}.
##' 
##' This is computed by simulating a large number of individuals \code{M} using
##' the maximum likelihood estimates of the fitted model and the function
##' \code{\link{sim.fmsm}}.  Therefore this requires a random sampling function
##' for the parametric survival model to be available: see the "Details"
##' section of \code{\link{sim.fmsm}}.  This will be available for all built-in
##' distributions, though users may need to write this for custom models.
##' 
##' Note the random sampling method for \code{flexsurvspline} models is
##' currently very inefficient, so that looping over the \code{M} individuals
##' will be very slow.
##' 
##' \code{\link{pmatrix.fs}} is a more efficient method based on solving the
##' Kolmogorov forward equation numerically, which requires the multi-state
##' model to be Markov.  No error or warning is given if running
##' \code{\link{pmatrix.simfs}} with a Markov model, but this is still invalid.
##' 
##' @param x A model fitted with \code{\link{flexsurvreg}}.  See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.  Additionally this should be semi-Markov, so that the time variable
##' represents the time since the last transition.  In other words the response
##' should be of the form \code{Surv(time,status)}. See the package vignette
##' for further explanation.
##' 
##' \code{x} can also be a list of models, with one component for each
##' permitted transition, as illustrated in \code{\link{msfit.flexsurvreg}}.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param t Time to predict state occupancy probabilities for.  This must be a
##' single number, unlike \code{\link{pmatrix.fs}}.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param ci Return a confidence interval calculated by simulating from the
##' asymptotic normal distribution of the maximum likelihood estimates.  This
##' is turned off by default, since two levels of simulation are required.  If
##' turned on, users should adjust \code{B} and/or \code{M} until the results
##' reach the desired precision.  The simulation over \code{M} is generally
##' vectorised, therefore increasing \code{B} is usually more expensive than
##' increasing \code{M}.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @param tcovs Predictable time-dependent covariates such as age, see
##' \code{\link{sim.fmsm}}.
##' @param M Number of individuals to simulate in order to approximate the
##' transition probabilities.  Users should adjust this to obtain the required
##' precision.
##' @param B Number of simulations from the normal asymptotic distribution used
##' to calculate variances.  Decrease for greater speed at the expense of
##' accuracy.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @return The transition probability matrix.  If \code{ci=TRUE}, there are
##' attributes \code{"lower"} and \code{"upper"} giving matrices of the
##' corresponding confidence limits.  These are formatted for printing but may
##' be extracted using \code{attr()}.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @seealso
##' \code{\link{pmatrix.fs}},\code{\link{sim.fmsm}},\code{\link{totlos.simfs}},
##' \code{\link{msfit.flexsurvreg}}.
##' @keywords models,survival
##' @examples
##' 
##' # BOS example in vignette, and in msfit.flexsurvreg
##' 
##' bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' 
##' # more likely to be dead (state 3) as time moves on, or if start with
##' # BOS (state 2)
##' 
##' pmatrix.simfs(bexp, t=5, trans=tmat)
##' pmatrix.simfs(bexp, t=10, trans=tmat)
##' 
##' # these results should converge to those in help(pmatrix.fs), as M
##' # increases here and ODE solving precision increases there, since with
##' # an exponential distribution, the semi-Markov model is the same as the
##' # Markov model.
##' @export
pmatrix.simfs <- function(x, trans, t=1, newdata=NULL, ci=FALSE,
                          tvar="trans", tcovs=NULL, M=100000, B=1000, cl=0.95)
{
    n <- nrow(trans)
    res <- matrix(0, nrow=n, ncol=n)
    if (length(t)>1) stop("\"t\" must be a single number")
    for (i in seq_len(n)){
        sim <- sim.fmsm(x=x, trans=trans, t=t, newdata=newdata,
                      start=i, M=M, tvar=tvar, tcovs=tcovs, debug=FALSE)
        last.st <- sim$st[,ncol(sim$st)]
        res[i,] <- prop.table(table(factor(last.st, levels=seq_len(n))))
    }
    if (ci){
        resci <- bootci.fmsm(x, B, n*n, match.call(), cl)
        resl <- matrix(resci[1,], nrow=n)
        resu <- matrix(resci[2,], nrow=n)
        attr(res, "lower") <- resl
        attr(res, "upper") <- resu
        class(res) <- "fs.msm.est"
    }  
    res
}



##' Expected total length of stay in specific states, from a fully-parametric,
##' semi-Markov multi-state model
##' 
##' The expected total time spent in each state for semi-Markov multi-state
##' models fitted to time-to-event data with \code{\link{flexsurvreg}}.  This
##' is defined by the integral of the transition probability matrix, though
##' this is not analytically possible and is computed by simulation.
##' 
##' This is computed by simulating a large number of individuals \code{M} using
##' the maximum likelihood estimates of the fitted model and the function
##' \code{\link{sim.fmsm}}.  Therefore this requires a random sampling function
##' for the parametric survival model to be available: see the "Details"
##' section of \code{\link{sim.fmsm}}.  This will be available for all built-in
##' distributions, though users may need to write this for custom models.
##' 
##' Note the random sampling method for \code{flexsurvspline} models is
##' currently very inefficient, so that looping over \code{M} will be very
##' slow.
##' 
##' The equivalent function for time-inhomogeneous Markov models is
##' \code{\link{totlos.fs}}.  Note neither of these functions give errors or
##' warnings if used with the wrong type of model, but the results will be
##' invalid.
##' 
##' @param x A model fitted with \code{\link{flexsurvreg}}.  See
##' \code{\link{msfit.flexsurvreg}} for the required form of the model and the
##' data.  Additionally this should be semi-Markov, so that the time variable
##' represents the time since the last transition.  In other words the response
##' should be of the form \code{Surv(time,status)}. See the package vignette
##' for further explanation.
##' 
##' \code{x} can also be a list of models, with one component for each
##' permitted transition, as illustrated in \code{\link{msfit.flexsurvreg}}.
##' @param trans Matrix indicating allowed transitions.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param t Maximum time to predict to.
##' @param start Starting state.
##' @param newdata A data frame specifying the values of covariates in the
##' fitted model, other than the transition number.  See
##' \code{\link{msfit.flexsurvreg}}.
##' @param ci Return a confidence interval calculated by simulating from the
##' asymptotic normal distribution of the maximum likelihood estimates.  This
##' is turned off by default, since two levels of simulation are required.  If
##' turned on, users should adjust \code{B} and/or \code{M} until the results
##' reach the desired precision.  The simulation over \code{M} is generally
##' vectorised, therefore increasing \code{B} is usually more expensive than
##' increasing \code{M}.
##' @param tvar Variable in the data representing the transition type. Not
##' required if \code{x} is a list of models.
##' @param tcovs Predictable time-dependent covariates such as age, see
##' \code{\link{sim.fmsm}}.
##' @param group Optional grouping for the states.  For example, if there are
##' four states, and \code{group=c(1,1,2,2)}, then \code{\link{totlos.simfs}}
##' returns the expected total time in states 1 and 2 combined, and states 3
##' and 4 combined.
##' @param M Number of individuals to simulate in order to approximate the
##' transition probabilities.  Users should adjust this to obtain the required
##' precision.
##' @param B Number of simulations from the normal asymptotic distribution used
##' to calculate variances.  Decrease for greater speed at the expense of
##' accuracy.
##' @param cl Width of symmetric confidence intervals, relative to 1.
##' @return The expected total time spent in each state (or group of states
##' given by \code{group}) up to time \code{t}, and corresponding confidence
##' intervals if requested.
##' @author Christopher Jackson \email{chris.jackson@@mrc-bsu.cam.ac.uk}.
##' @seealso
##' \code{\link{pmatrix.simfs}},\code{\link{sim.fmsm}},\code{\link{msfit.flexsurvreg}}.
##' @keywords models,survival
##' @examples
##' 
##' # BOS example in vignette, and in msfit.flexsurvreg
##' bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp")
##' tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
##' 
##' # predict 4 years spent without BOS, 3 years with BOS, before death
##' # As t increases, this should converge
##' totlos.simfs(bexp, t=10, trans=tmat)
##' totlos.simfs(bexp, t=1000, trans=tmat)
##' @export
totlos.simfs <- function(x, trans, t=1, start=1, newdata=NULL, ci=FALSE,
                          tvar="trans", tcovs=NULL, group=NULL, M=100000, B=1000, cl=0.95)
{
    if (length(t)>1) stop("\"t\" must be a single number")
    sim <- sim.fmsm(x=x, trans=trans, t=t, newdata=newdata,
                    start=start, M=M, tvar=tvar, tcovs=tcovs, debug=FALSE)
    dt <- diff(t(cbind(sim$t, t)))
    st <- factor(t(sim$st), levels=1:nrow(trans))
    res <- tapply(dt, st, sum) / M
    res[is.na(res)] <- 0
    if (!is.null(group)) {
        if(length(group) != nrow(trans))
            stop("\"group\" must be a vector of length ",nrow(trans), " = number of states")
        res <- tapply(res, group, sum)
    }
    if (ci){
        resci <- bootci.fmsm(x, B, length(res), match.call(), cl)
        resl <- resci[1,]
        resu <- resci[2,]
        names(resl) <- names(resu) <- t
        attr(res, "lower") <- resl
        attr(res, "upper") <- resu
        class(res) <- "fs.msm.est"
        res <- cbind(est=res, L=resl, U=resu)
    }
    res
}
