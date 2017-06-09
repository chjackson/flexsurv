## Get density, probability, hazard and cumulative hazard functions
## and return them as a list of functions

form.dp <- function(dlist, dfns, integ.opts){
    
    ## TODO check for format of dfn (args x, log) 
    ## FIXME bug if object called d is found in global env
    ## check for existence in current frame.  inherits false? 
    
    name <- dlist$name
    hname <- paste0("h",name); Hname <- paste0("H",name)
    dname <- paste0("d",name); pname <- paste0("p",name)
    rmstname <- paste0("rmst_",name)
    meanname <- paste0("mean_",name)
    qname <- paste0("q",name)
    rname <- paste0("r",name)
    if (is.function(dfns$d)) d <- dfns$d
    if (is.function(dfns$p)) p <- dfns$p    
    if (is.function(dfns$h)) h <- dfns$h
    if (is.function(dfns$H)) H <- dfns$H
    if (is.function(dfns$r)) r <- dfns$r
    if (is.function(dfns$q)) q <- dfns$q
    if (is.function(dfns$mean)) meanf <- dfns$mean
    if (is.function(dfns$rmst)) rmst <- dfns$rmst
    if (!exists("h", inherits=FALSE)){
        if (exists(hname)) h <- get(hname)
        else {
            if (!exists("d")){
                if (exists(dname)) d <- get(dname)
                else stop("Neither density function \"",dname,
                          "\" nor hazard function \"", hname, "\" found")
            }
            if (!exists("p")){
                if (exists(pname)) p <- get(pname)
                else {
                    message("Forming cumulative distribution function...")
                    p <- integrate.dh(d, dlist, integ.opts, what="density")
                }
            }
            h <- function(x, ...){
                d(x,...)/(1 - p(x,...))
            }
        }
    }
    if (!exists("H", inherits=FALSE)){
        if (exists(Hname)) H <- get(Hname)
        else {
            if (!exists("p")) { if (exists(pname)) p <- get(pname) }
            if (exists("p")){
                H <- function(x, ...){
                    -log(1 - p(x, ...))
                }
            } else {
                message("Forming integrated hazard function...")
                H <- integrate.dh(h, dlist, integ.opts, what="hazard") 
            }
        }
    }
    if (!exists("p", inherits=FALSE)){
        if (exists(pname)) p <- get(pname)
        else {
            p <- function(q, ...) {
                ret <- 1 - exp(-H(q, ...))
#                ret[q==Inf] <- 1 # should have been handled already in cum.fn
#                ret[q==0] <- 0
                ret
### TODO special values in other functions
            }
        }
    }
    if (!exists("q", inherits=FALSE)){
      if (exists(qname)) q <- get(qname)
      else {
        # giving this another name to avoid scoping issues
        # w/ name p also being an argument to q functions
        pfun <- p
        q <- function(p, ...) qgeneric(pfun, p)
      } 
    }
    if (!exists("d", inherits=FALSE)){
        if (exists(dname)) d <- get(dname)
        else { 
            d <- function(x, log=FALSE, ...) {
                if (log)
                    log(h(x,...)) + log(1 - p(x, ...))
                else h(x,...) * (1 - p(x, ...))
            }
        }
    }
    if (!exists("rmst", inherits=FALSE)){
      if (exists(rmstname)) rmst <- get(rmstname)
      else {
        message("Forming integrated rmst function...")
        rmst <- function(t, start=0, ...) rmst_generic(p, t=t, start=start, ...) 
      }
    }
    if (!exists("meanf", inherits=FALSE)){
      if (exists(rmstname)) meanf <- get(meanname)
      else {
        message("Forming integrated mean function...")
        meanf <- function(start=0, ...) rmst(t=Inf, start=start, ...)
      }
    }
    if (!exists("r", inherits=FALSE)){
        if (exists(rname)) r <- get(rname)
        else r <- NULL
        ## random sampling function is currently only used for multi-state models
    }
    ## Check for existence of derivative functions
    ## conventionally called DLd, DLs
    ## if dfns$deriv set to FALSE on entry, derivatives not available
    if (is.function(dfns$DLd)) DLd <- dfns$DLd
    else if (is.null(dfns$deriv) && exists(paste0("DLd",name)))
        DLd <- get(paste0("DLd",name))
    else DLd <- NULL
    if (is.function(dfns$DLS)) DLS <- dfns$DLS
    else if (is.null(dfns$deriv) && exists(paste0("DLS",name)))
        DLS <- get(paste0("DLS",name))
    else DLS <- NULL
    
    list(p=p, d=d, h=h, H=H, r=r, DLd=DLd, DLS=DLS, rmst=rmst, mean= meanf,
         q=q, deriv = !is.null(DLd) && !is.null(DLS))
}


## Produce cumulative version of hazard function or density function
## by numerical integration

integrate.dh <- function(fn, dlist, integ.opts, what="dens"){

    cum.fn <- function(q, ...){
        args <- list(...)
        pars <- as.list(dlist$pars)
        names(pars) <- dlist$pars
        args.done <- numeric()
        ## if argument is unnamed, assume it is supplied in the default order
        for (i in seq(along=dlist$pars)){
            if(any(names(args)==dlist$pars[i])) {
                pars[[i]] <- args[[dlist$pars[i]]]
                args.done <- c(args.done, match(dlist$pars[i], names(args)))
            } else {
                pars[[i]] <- args[[i]]
                args.done <- c(args.done, i)
            }
        }
        ## any auxiliary arguments not in main distribution parameters
        rest <- args[setdiff(seq_along(args), args.done)] 
        ## replicate all arguments to have the length of the longest one (=n)
        n <- max(sapply(c(list(q),pars), length))
        q <- rep(q, length=n)
        for (i in seq_along(pars)) pars[[i]] <- rep(pars[[i]], length=n)
        ret <- numeric(n)                       
        du <- function(u, ...)fn(u,...)
        ## then return a vector of length n
        for (i in 1:n){
            parsi <- lapply(pars, function(x)x[i])
            int.args <- c(list(f=du, lower=0, upper=q[i]), parsi, rest, integ.opts)
            if (q[i]==0) ret[i] <- 0
            else if (q[i]==Inf) {
                if (what=="density") ret[i] <- 1
                else if (what=="hazard") ret[i] <- Inf
            }
            else {
                int <- try(do.call("integrate", int.args))
#                 if (inherits(int, "try-error")) browser()
                ret[i] <- int$value
            }
        }
        ret
    }

    cum.fn
}
