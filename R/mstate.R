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

pars.fmsm <- function(x, trans, newdata=NULL, tvar="trans")
{
    if (is.flexsurvlist(x)){
        ntr <- length(x) # number of allowed transitions
        if (ntr != length(na.omit(as.vector(trans)))) stop(sprintf("x is a list of %s flexsurvreg objects, but trans indicates %s transitions", ntr, length(na.omit(as.vector(trans)))))
        basepar <- matrix(nrow=ntr, ncol=length(x[[1]]$dlist$pars), dimnames=list(NULL,x[[1]]$dlist$pars))
        for (i in 1:ntr){
            X <- if (x[[i]]$ncovs==0) matrix(0) else form.model.matrix(x[[i]], as.data.frame(newdata))
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
        dat <- as.list(newdata[transi,])
    }
    for (i in tcovs) { dat[[i]] <- dat[[i]] + 0  } # t}
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

pmatrix.simfs <- function(x, trans, t=1, newdata=NULL, ci=FALSE,
                          tvar="trans", tcovs=NULL, M=100000, B=1000, cl=0.95)
{
    n <- nrow(trans)
    res <- matrix(0, nrow=n, ncol=n)
    if (length(t)>1) stop("\"t\" must be a single number")
    for (i in seq_len(n)){
        sim <- sim.fmsm(x=x, trans=trans, t=t, newdata=newdata,
                      start=i, M=M, tvar="trans", tcovs=tcovs, debug=FALSE)
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

totlos.simfs <- function(x, trans, t=1, start=1, newdata=NULL, ci=FALSE,
                          tvar="trans", tcovs=NULL, group=NULL, M=100000, B=1000, cl=0.95)
{
    if (length(t)>1) stop("\"t\" must be a single number")
    sim <- sim.fmsm(x=x, trans=trans, t=t, newdata=newdata,
                    start=start, M=M, tvar="trans", tcovs=tcovs, debug=FALSE)
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
