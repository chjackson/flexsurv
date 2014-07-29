## Get density, probability, hazard and cumulative hazard functions
## and return them as a list of functions

form.dp <- function(dlist, dfns, integ.opts){
    
    ## TODO check for format of dfn (args x, log) 
    
    name <- dlist$name
    hname <- paste0("h",name); Hname <- paste0("H",name)
    dname <- paste0("d",name); pname <- paste0("p",name)
    if (is.function(dfns$d)) d <- dfns$d
    if (is.function(dfns$p)) p <- dfns$p    
    if (is.function(dfns$h)) h <- dfns$h
    if (is.function(dfns$H)) H <- dfns$H
    if (!exists("h")){
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
                    p <- function(q, ...){
                        args <- list(...)
                        pars <- as.list(dlist$pars)
                        for (i in seq(along=dlist$pars))
                            pars[[i]] <- if(any(names(args)==dlist$pars[i])) args[[dlist$pars[i]]] else args[[i]]
                        ## TODO do we pass through any other arguments.  if so, easiest to name them beforehand and treat in above line like extra dlist$pars
                        n <- max(sapply(c(list(q),pars), length))
                        q <- rep(q, length=n)
                        for (i in seq_along(pars)) pars[[i]] <- rep(pars[[i]], length=n)
                        ret <- numeric(n)                       
                        du <- function(u, ...)d(u,...)
                        for (i in 1:n){
                            parsi <- lapply(pars, function(x)x[i])
                            int.args <- c(list(f=du, lower=0, upper=q[i]), parsi, integ.opts)
                            int <- do.call("integrate", int.args)
                            ret[i] <- int$value
                        }
                        ret
###                        int.args <- c(list(f=du, lower=0, upper=q), integ.opts, list(...))
###                        do.call("integrate", int.args)$value
                    }
###                    p <- Vectorize(p)
                }
            }
            h <- function(x, ...){
                d(x,...)/(1 - p(x,...))
            }
        }
    }
    if (!exists("H")){
        if (exists(Hname)) H <- get(Hname)
        else {
            if (!exists("p")) { if (exists(pname)) p <- get(pname) }
            if (exists("p")){
                H <- function(x, ...){
                    -log(1 - p(x, ...))
                }
            } else {
                message("Forming integrated hazard function...")
                H <- function(x, ...)
                {
                    args <- list(...)
                    pars <- as.list(dlist$pars)
                    for (i in seq(along=dlist$pars))
                        pars[[i]] <- if(any(names(args)==dlist$pars[i])) args[[dlist$pars[i]]] else args[[i]]
                    ## TODO do we pass through any other arguments.  if so, easiest to name them beforehand and treat in above line like extra dlist$pars
                    n <- max(sapply(c(list(x),pars), length))
                    x <- rep(x, length=n)
                    for (i in seq_along(pars)) pars[[i]] <- rep(pars[[i]], length=n)
                    ret <- numeric(n)
                    hu <- function(u, ...)h(u,...)
                    for (i in 1:n){
                        parsi <- lapply(pars, function(x)x[i])
                        int.args <- c(list(f=hu, lower=0, upper=x[i]), parsi, integ.opts)
                        int <- do.call("integrate", int.args)
                        ret[i] <- int$value
                    }
                    ret
                }
                ## could copy same code to get p from d
            }
        }
    }
    if (!exists("p")){
        if (exists(pname)) p <- get(pname)
        else {
            p <- function(q, ...) {
                1 - exp(-H(q, ...))
            }
        }
    }
    if (!exists("d")){
        if (exists(dname)) d <- get(dname)
        else { 
            d <- function(x, log=FALSE, ...) {
                if (log)
                    log(h(x,...)) + log(1 - p(x, ...))
                else h(x,...) * (1 - p(x, ...))
            }
        }
    }

    ## Check for existence of derivative functions
    ## conventionally called DLd, DLs
    if (is.function(dfns$DLd)) DLd <- dfns$DLd
    else if (exists(paste0("DLd",name))) DLd <- get(paste0("DLd",name))
    else DLd <- NULL
    if (is.function(dfns$DLS)) DLS <- dfns$DLS
    else if (exists(paste0("DLS",name))) DLS <- get(paste0("DLS",name))
    else DLS <- NULL

    list(p=p, d=d, h=h, H=H, DLd=DLd, DLS=DLS,
         deriv = !is.null(DLd) && !is.null(DLS))
}
