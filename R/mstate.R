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
    newdata <- form.msm.newdata(object, newdata=newdata, tvar=tvar, trans=trans)
    X <- form.model.matrix(object, newdata)
    Haz <- summary(object, type="cumhaz", t=t, X=X, ci=FALSE)
    Haz <- do.call("rbind",Haz[seq_along(tr)])
    Haz$trans <- rep(seq_along(tr), each=length(t))
    names(Haz)[names(Haz)=="est"] <- "Haz"
    res <- list(Haz=Haz, trans=trans)
    if (variance & !is.na(object$cov[1])){
        boot <- array(dim=c(B, length(t), ntr))
        for (i in seq_along(tr))
            boot[,,i] <- normbootfn.flexsurvreg(object, t=t, start=0, X=X[i,,drop=FALSE], B=B,
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

pmatrix.fs <- function(x, trans, t=1, newdata=NULL, ci=FALSE,
                       tvar="trans", sing.inf=1e+10, B=1000, cl=0.95, ...){
    newdata <- form.msm.newdata(x, newdata=newdata, tvar=tvar, trans=trans)
    X <- form.model.matrix(x, newdata)   
    n <- nrow(trans)
    dp <- function(t, y, parms, ...){
        P <- matrix(y, nrow=n, ncol=n)
        haz <- numeric(n)
        for (i in 1:n){
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
    basepar <- add.covs(x, pars=x$res.t[x$dlist$pars,"est"], beta=x$res.t[x$covpars,"est"], X=X)
    res <- ode(y=diag(n), times=c(0,t), func=dp, parms=list(par=basepar), ...)[-1,-1]
    res <- lapply(split(res,1:nt), function(x)matrix(x,nrow=n))
    names(res) <- t
    if (ci){
        sim <- normboot.flexsurvreg(x=x, B=B, X=X)
        if (!is.list(sim)) sim <- list(sim)        
        res.rep <- array(NA_real_, dim=c(B, nt, n*n))
        for (i in 1:B){
            pari <- do.call("rbind", lapply(sim, function(x)x[i,,drop=FALSE]))
            res.rep[i,,] <- ode(y=diag(n), times=c(0,t), func=dp, parms=list(pars=pari),...)[-1,-1]
        }
        resci <- apply(res.rep, c(2,3), quantile, c((1-cl)/2, 1 - (1-cl)/2))
        resl <- lapply(split(resci[1,,],1:nt), function(x)matrix(x,nrow=n))
        resu <- lapply(split(resci[2,,],1:nt), function(x)matrix(x,nrow=n))
        names(resl) <- names(resu) <- t
        for (i in 1:nt){
            attr(res[[i]], "lower") <- resl[[i]]
            attr(res[[i]], "upper") <- resu[[i]]
            class(res[[i]]) <- "fs.msm.est"
        }
    }
    if(nt==1) res[[1]] else res
}

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
