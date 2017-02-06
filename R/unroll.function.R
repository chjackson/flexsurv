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
