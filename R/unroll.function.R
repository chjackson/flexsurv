##' Convert a function with matrix arguments to a function with vector
##' arguments.
##' 
##' Given a function with matrix arguments, construct an equivalent function
##' which takes vector arguments defined by the columns of the matrix.  The new
##' function simply uses \code{cbind} on the vector arguments to make a matrix,
##' and calls the old one.
##' 
##' 
##' @param mat.fn A function with any number of arguments, some of which are
##' matrices.
##' @param \dots A series of other arguments.  Their names define which
##' arguments of \code{mat.fn} are matrices.  Their values define a vector of
##' strings to be appended to the names of the arguments in the new function.
##' For example
##' 
##' \code{fn <- unroll.function(oldfn, gamma=1:3, alpha=0:1)}
##' 
##' will make a new function \code{fn} with arguments
##' \code{gamma1},\code{gamma2},\code{gamma3},\code{alpha0},\code{alpha1}.
##' 
##' Calling
##' 
##' \code{fn(gamma1=a,gamma2=b,gamma3=c,alpha0=d,alpha1=e)}
##' 
##' should give the same answer as
##' 
##' \code{oldfn(gamma=cbind(a,b,c),alpha=cbind(d,e))}
##' @return The new function, with vector arguments.
##' @section Usage in \pkg{flexsurv}:
##' 
##' This is used by \code{\link{flexsurvspline}} to allow spline models, which
##' have an arbitrary number of parameters, to be fitted using
##' \code{\link{flexsurvreg}}.
##' 
##' The ``custom distributions'' facility of \code{\link{flexsurvreg}}
##' expects the user-supplied probability density and distribution
##' functions to have one explicitly named argument for each scalar
##' parameter, and given R vectorisation, each of those arguments
##' could be supplied as a vector of alternative parameter values.
##' 
##' However, spline models have a varying number of scalar parameters,
##' determined by the number of knots in the spline.
##' \code{\link{dsurvspline}} and \code{\link{psurvspline}} have an
##' argument called \code{gamma}.  This can be supplied as a matrix,
##' with number of columns \code{n} determined by the number of knots
##' (plus 2), and rows referring to alternative parameter values.  The
##' following statements are used in the source of
##' \code{flexsurvspline}: \preformatted{ dfn <-
##' unroll.function(dsurvspline, gamma=0:(nk-1)) pfn <-
##' unroll.function(psurvspline, gamma=0:(nk-1)) }
##' 
##' to convert these into functions with arguments \code{gamma0},
##' \code{gamma1},\ldots{},\code{gamman}, corresponding to the columns
##' of \code{gamma}, where \code{n = nk-1}, and with other arguments
##' in the same format.
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##' @seealso \code{\link{flexsurvspline}},\code{\link{flexsurvreg}}
##' @examples
##' 
##' fn <- unroll.function(ncol, x=1:3)
##' fn(1:3, 1:3, 1:3) # equivalent to...
##' ncol(cbind(1:3,1:3,1:3))
##' @export
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
