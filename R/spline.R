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
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    knots <- apply(knots, 2, function(x)rep(x,length=nret))
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=nret)
    if(!is.matrix(knots)) knots <- matrix(knots, nrow=nret)
    if (ncol(gamma) != ncol(knots)) {
        stop("length of gamma should equal number of knots")
    }
    scale <- match.arg(scale, c("hazard","odds","normal"))
    if (deriv){
        ret <- matrix(nrow=nret, ncol=ncol(gamma))
        ret[is.na(q),] <- NA
        ret[!is.na(q) & q <= 0,] <- 0
    } else {
        ret <- numeric(nret)
        ret[is.na(q)] <- NA
        ret[!is.na(q) & q <= 0] <- 0
    }
    ind <- !is.na(q) & q > 0
    q <- q[ind]; gamma <- gamma[ind,,drop=FALSE]
    list(ret=ret, gamma=gamma, q=q, scale=scale, ind=ind)
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

dsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0, log=FALSE){
    d <- dbase.survspline(q=x, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]]); x <- q
    eta <- rowSums(basis(knots, log(x)) * gamma) + as.numeric(X %*% beta) + offset # log cumulative hazard/odds
    eeta <- exp(ldlink(scale)(eta))
    ret[ind][eeta==0] <- 0
    ret[ind][is.nan(eeta)] <- NaN
    ind2 <- !(eeta==0 || is.nan(eeta))
    x <- x[ind2]; gamma <- gamma[ind2,,drop=FALSE]; eeta <- eeta[ind2]
    ind <- ind & ind2
    dsum <- rowSums(dbasis(knots, log(x)) * gamma)  # ds/dx
    ret[ind] <- 1 / x * dsum * eeta
    ## derivative of log cum haz cannot be negative by definition, but
    ## optimisation doesn't constrain gamma to respect this, so set
    ## likelihood to zero then (assuming at least one death)
    ret[ind][ret[ind]<=0] <- 0
    if (log) {ret <- log(ret)}
    as.numeric(ret)
}

Slink <- function(scale){
    switch(scale,
           hazard=function(x){exp(-exp(x))},
           odds=function(x){1 / (1 + exp(x))},
           normal=function(x){pnorm(-x)}
           )
}

## cumulative distribution function

psurvspline <- function(q, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0, lower.tail=TRUE, log.p=FALSE){
    d <- dbase.survspline(q=q, gamma=gamma, knots=knots, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    eta <- rowSums(basis(knots, log(q)) * gamma) + as.numeric(X %*% beta) + offset
    surv <- Slink(scale)(eta)
    ret[ind] <- as.numeric(1 - surv)
    ret[ind][q==0] <- 0
    ret[ind][q==Inf] <- 1
    if (!lower.tail) ret <- 1 - ret
    if (log.p) ret <- log(ret)
    ret
}

qsurvspline <- function(p, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0, lower.tail=TRUE, log.p=FALSE){
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    qgeneric(psurvspline, p=p, gamma=gamma, beta=beta, X=X, knots=knots, scale=scale, offset=offset)
}

rsurvspline <- function(n, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0){
    if (length(n) > 1) n <- length(n)
    ret <- qsurvspline(p=runif(n), gamma=gamma, beta=beta, X=X, knots=knots, scale=scale, offset=offset)
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

Hsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0){
# TODO base on dp H=-log(1-p) to get error handling etc
    match.arg(scale, c("hazard","odds","normal"))
    eta <- basis(knots, log(x)) %*% gamma + as.numeric(X %*% beta) + offset
    as.numeric(Hlink(scale)(eta))
}

hlink <- function(scale){
    switch(scale,
           hazard=function(x){exp(x)},
           odds=function(x){plogis(x)},
           normal=function(x){dnorm(-x)/pnorm(-x)}
           )
}

## hazard function

hsurvspline <- function(x, gamma, beta=0, X=0, knots=c(-10,10), scale="hazard", offset=0){
# TODO base on dp:  d / (1 - p)
    match.arg(scale, c("hazard","odds","normal"))
    eta <- basis(knots, log(x)) %*% gamma + as.numeric(X %*% beta) + offset
    eeta <- hlink(scale)(eta)
    haz <- 1 / x * dbasis(knots, log(x)) %*% gamma * eeta
    as.numeric(haz)
}

basis <- function(knots, x) {
    nx <- length(x)
    if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
    nk <- ncol(knots)
    b <- matrix(nrow=length(x), ncol=nk)
    if (nk>0){
        b[,1] <- 1
        b[,2] <- x
    }
    if (nk>2) {
        lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
        for (j in 1:(nk-2)) {
            b[,j+2] <- pmax(x - knots[,j+1], 0)^3 - lam[,j+1]*pmax(x - knots[,1], 0)^3 -
                (1 - lam[,j+1])*pmax(x - knots[,nk], 0)^3
        }
    }
    b
}

dbasis <- function(knots, x) {
    nx <- length(x) 
   if (!is.matrix(knots)) knots <- matrix(rep(knots, nx), byrow=TRUE, ncol=length(knots))
    nk <- ncol(knots)
    b <- matrix(nrow=length(x), ncol=nk)
    if (nk>0){
        b[,1] <- 0
        b[,2] <- 1
    }
    if (nk>2) {
        lam <- (knots[,nk] - knots)/(knots[,nk] - knots[,1])
        for (j in 3:nk) {
            b[,j] <- 3*pmax(x - knots[,j-1], 0)^2 - 3*lam[,j-1]*pmax(x - knots[,1], 0)^2 -
                3*(1 - lam[,j-1])*pmax(x - knots[,nk], 0)^2
        }
    }
    b
}

fss <- function(x, knots) { basis(knots, x) }
dfss <- function(x, knots) { dbasis(knots, x) }

## Given a function with matrix arguments (e.g. matrix.fn <-
## function(..., gamma, knots), where "gamma" and "knots" have 2
## columns each, say), this makes an equivalent function with sets of
## vector arguments, defined something like

## function (..., gamma1, gamma2, knots1, knots2)
##{
##    gamma <- do.call("cbind", list(gamma1, gamma2))
##    knots <- do.call("cbind", list(knots1, knots2))
##    do.call(matrix.fn, list(..., gamma = gamma, knots = knots))
##}

## where ... represent the arguments which are unchanged
## Called as, e.g. 
## vector.fn <- unroll.function(matrix.fn, gamma=2, knots=2)

## Used, e.g. to return d or p functions for spline distribution with
## a variable number nk of parameters gamma1, gamma2, ...  in a format
## that can be used in flexsurvreg

unroll.function <- function(mat.fn, ...){
    fargs <- formals(mat.fn)
    vargs <- list(...) # list of names and numbers
    if (length(vargs)==0) return(mat.fn)
    badargs <- paste0("\"",names(vargs)[!(names(vargs) %in% names(fargs))],"\"")
    argerr <- if (length(badargs) > 1) "arguments" else "an argument"
    if (badargs!="\"\"")
        stop("\"",deparse(substitute(mat.fn)),"\" does not have ",argerr, " named ", paste(badargs,collapse=","))
    sargs <- fargs[setdiff(names(fargs), names(vargs))] # arguments to not change
    ## converts e.g. list(gamma=2,knots=2) to c("gamma1","gamma2","knots1","knots2")
    unames <- mapply(function(x,y)paste0(x, y), names(vargs), vargs)
    unamesv <- as.vector(unames)
    ## makes an alist(gamma1=,gamma2=,knots1=,knots2=)
    uargs <- eval(parse(text=paste0("alist(",paste(unamesv, collapse="=, "),"=)")))
    args <- as.pairlist(c(sargs,uargs))
    ## copy the old function definition into the body of the new one
    ## so it remains visible
    dpa <- deparse(args(mat.fn))
    basefn.lines <- paste("base.fn <-",
                          paste(paste(dpa[-length(dpa)]), collapse="\n"),
                          paste(deparse(body(mat.fn)), collapse="\n"))
    ## build the function body 
    cbind.lines <- character(length=length(vargs))
    for (i in seq_along(vargs)){
        ## make statements like:  gamma <- do.call("cbind", list(gamma1, gamma2))
        cbind.lines[i] <- sprintf("%s <- do.call(\"cbind\", list(%s))",
                                  names(vargs)[i], paste(unames[,i],collapse=","))
    }
    ## make the return statement of the new function
    sargsn <- paste(paste(names(sargs),names(sargs),sep="="), collapse=", ")
    vargsn <- paste(paste0(names(vargs),"=",names(vargs)), collapse=", ")
    ret.line <- if (sargsn=="")
        sprintf("do.call(base.fn, list(%s))",vargsn) else
        sprintf("do.call(base.fn, list(%s, %s))",sargsn,vargsn)
    body <- as.call(c(as.name("{"),
                      parse(text = paste(
                            basefn.lines,
                            paste(cbind.lines, collapse="\n"),
                            ret.line,
                            sep="\n")
                            )))
    ## thanks to http://adv-r.had.co.nz/Expressions.html#pairlists for this
    res <- eval(call("function", args, body), envir=parent.frame())
    environment(res) <- environment(mat.fn)
    res
}

flexsurv.splineinits <- function(t=NULL, mf, mml, aux)
{
    Y <- check.flexsurv.response(model.extract(mf, "response"))
    weights <- model.extract(mf, "weights")
    dead <- Y[,"status"]==1
    X <- mml[[1]][dead,-1,drop=FALSE]
    Y <- Y[dead,,drop=FALSE]
    wt <- weights[dead]
    ## using coxph on original formula followed by survfit.coxph fails
    ## due to scoping
    form <- paste("Surv(Y[,\"start\"], Y[,\"stop\"], Y[,\"status\"]) ~ ")
    if (ncol(X)>0)
        form <- paste(form, paste(paste("X[,",1:ncol(X),"]",sep=""), collapse=" + "))
    else form <- paste(form, "1")

    cox <- coxph(as.formula(form), weights=wt)
    surv <- survfit(cox, data=cbind(Y, X))
    surv <- surv$surv[match(Y[,"stop"], surv$time)]
    if (aux$scale=="hazard")
        logH <- log(-log(surv))
    else if (aux$scale=="odds")
        logH <- log((1 - surv)/surv)
    else if (aux$scale=="normal")
        logH <- qnorm(1 - surv)
    b <- basis(aux$knots, log(Y[,"time"]))
    form <- paste("logH ~ ",
                  paste(paste("b[,",2:ncol(b),"]",sep=""), collapse=" + "))
    if (ncol(X)>0)
        form <- paste(form, "+", paste(paste("X[,",1:ncol(X),"]",sep=""), collapse=" + "))
    ## include interactions between covs and spline basis
    nints <- length(mml)-1
    inter <- lapply(mml[-1], function(x)NULL)
    for (i in seq_along(inter)){
        Xi <- mml[[i+1]][dead,-1,drop=FALSE]
        ncovsi <- if (is.null(Xi)) 0 else ncol(Xi)
        nam <- paste0("X",i)
        assign(nam, Xi)
        inter[[i]] <- if (ncovsi==0) NULL else paste0("b[,",i+1,"]:", paste0(nam,"[,",seq_len(ncovsi),"]"))
    }
    inter <- paste(unlist(inter), collapse=" + ")
    if (inter != "") form <- paste(form, inter, sep="+")
    inits <- coef(lm(as.formula(form), weights=wt))
    inits
}

flexsurvspline <- function(formula, data, k=0, knots=NULL, bknots=NULL, scale="hazard", ...){
    ## Get response matrix from the formula.  Only need this to obtain
    ## default knots.  Largely copied from flexsurvreg - ideally
    ## should be in separate function, but can't make scoping work.   
    call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), names(call), nomatch = 0)
    if (indx[1] == 0)
        stop("A \"formula\" argument is required")
    temp <- call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    ## remove the predictors
    f2 <- as.formula(gsub("(.*~).*", "\\1 1", Reduce(paste, deparse(formula))))
    environment(f2) <- environment(formula)
    temp[["formula"]] <- f2
    if (missing(data)) temp[["data"]] <- environment(formula)
    m <- eval(temp, parent.frame())
    Y <- check.flexsurv.response(model.extract(m, "response"))
    
    dtimes <- Y[,"stop"][Y[,"status"]==1]
    if (is.null(knots)) {
        is.wholenumber <-
            function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
        if (is.null(k)) stop("either \"knots\" or \"k\" must be specified")
        if (!is.numeric(k)) stop("k must be numeric")
        if (!is.wholenumber(k) || (k<0)) stop("number of knots \"k\" must be a non-negative integer")
        knots <- quantile(log(dtimes), seq(0, 1, length=k+2)[-c(1,k+2)])
    }
    else {
        if (!is.numeric(knots)) stop("\"knots\" must be a numeric vector")
        minlogtime <- min(log(Y[,"stop"]))
        if (any(knots <= minlogtime)) {
            badknots <- knots[knots < min(log(Y[,"stop"]))]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " less than or equal to minimum log time ", minlogtime)
        }
        maxlogtime <- max(log(Y[,"stop"]))
        if (any(knots >= maxlogtime)) {
            badknots <- knots[knots > maxlogtime]
            plural <- if (length(badknots) > 1) "s" else ""
            stop("knot", plural, " ", paste(badknots,collapse=", "), " greater than or equal to maximum log time ", maxlogtime)
        }
    }
    if (is.null(bknots)) {
        bknots <- c(min(log(dtimes)), max(log(dtimes)))
        if (bknots[1] == bknots[2])
            warning("minimum and maximum log death times are the same: knot and boundary knot locations should be supplied explicitly")
    } else 
        if (!is.numeric(bknots) || (length(bknots) !=2) ) stop("boundary knots should be a numeric vector of length 2")
    knots <- c(bknots[1], knots, bknots[2])
    
    nk <- length(knots)
    custom.fss <- list(
        name = "fn", # unused, d,p functions passed through
        pars = c(paste0("gamma",0:(nk-1))),
        location = c("gamma0"),
        transforms = rep(c(identity), nk), inv.transforms=rep(c(identity), nk),
        inits = flexsurv.splineinits
        )
    aux <- list(knots=knots, scale=scale)
    dfn <- unroll.function(dsurvspline, gamma=0:(nk-1))
    pfn <- unroll.function(psurvspline, gamma=0:(nk-1))
    Ddfn <- if (scale=="normal") NULL else unroll.function(DLdsurvspline, gamma=0:(nk-1))
    DSfn <- if (scale=="normal") NULL else unroll.function(DLSsurvspline, gamma=0:(nk-1))
    args <- c(list(formula=formula, data=data, dist=custom.fss,
                   dfns=list(d=dfn,p=pfn,DLd=Ddfn,DLS=DSfn), aux=aux), list(...))
    ret <- do.call("flexsurvreg", args) # faff to make ... args work within functions
    ret <- c(ret, list(k=length(knots) - 2, knots=knots, scale=scale))
    ret$call <- call
    class(ret) <- "flexsurvreg"
    ret
}

## do S3 class stuff: label the object as flexsurvspline, inheriting from flexsurvreg?
