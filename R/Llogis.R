dllogis <- function(x, shape=1, scale=1, log = FALSE)
{
    d <- dbase("llogis", log=log, x=x, shape=shape, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    logdens <- log(shape) - log(scale) + (shape-1)*(log(x) - log(scale)) - 2*log(1 + (x/scale)^shape)
    ret[ind] <- if (log) logdens else exp(logdens)
    ret
}

pllogis <- function(q, shape=1, scale=1, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("llogis", lower.tail=lower.tail, log=log.p, q=q, shape=shape, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    prob <- 1 - 1/(1 + (q/scale)^shape)
    if (!lower.tail) prob <- 1 - prob
    if (log.p) prob <- log(prob)
    ret[ind] <- prob
    ret
}

qllogis <- function(p, shape=1, scale=1, lower.tail = TRUE, log.p = FALSE) {
    d <- dbase("llogis", lower.tail=lower.tail, log=log.p, p=p, shape=shape, scale=scale)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    ret[ind] <- exp(qlogis(p, log(scale), 1/shape, lower.tail, log.p))
    ret
}

rllogis <- function(n, shape=1, scale=1){
    r <- rbase("llogis", n=n, shape=shape, scale=scale)
    for (i in seq_along(r)) assign(names(r)[i], r[[i]])
    ret[ind] <- qllogis(p=runif(sum(ind)), shape=shape, scale=scale)
    ret
} 

hllogis <- function(x, shape=1, scale=1, log = FALSE) 
{
    h <- dbase("llogis", log=log, x=x, shape=shape, scale=scale)
    for (i in seq_along(h)) assign(names(h)[i], h[[i]])
    if (log) ret[ind] <- log(shape) - log(scale) + (shape-1)*(log(x) - log(scale)) - log(1 + (x/scale)^shape)
    else ret[ind] <- (shape/scale) * (x/scale)^{shape-1} / (1 + (x/scale)^shape)
    ret
}

Hllogis <- function(x, shape=1, scale=1, log = FALSE) 
{
    ret <- - pllogis(x, shape, scale, lower.tail=FALSE, log.p=TRUE)
    if (log) ret <- log(ret)
    ret
}

check.llogis <- function(shape, scale){
    ret <- rep(TRUE, length(shape))
    if (any(shape<=0)) {warning("Non-positive shape parameter"); ret[shape<=0] <- FALSE}
    if (any(scale<=0)) {warning("Non-positive scale parameter"); ret[scale<=0] <- FALSE}
    ret
}

DLdllogis <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    ilt <- tss / (1 + tss)
    res[,1] <- 1 + (1 - 2*ilt)*shape*log(t/scale)
    res[,2] <- - shape + 2*ilt*shape
    res
}

DLSllogis <- function(t, shape, scale){
    res <- matrix(nrow=length(t), ncol=2)
    tss <- (t/scale)^shape
    ilt <- tss / (1 + tss)    
    res[,1] <- ifelse(t==0, 0, -ilt*log(t/scale)*shape)
    res[,2] <- shape*ilt
    res
}

mean.llogis <- function(shape=1, scale=1){
    ifelse(shape > 1,
       {b <- pi/shape
        scale*b / sin(b)},
           NaN)          
}

var.llogis <- function(shape=1, scale=1){
    ifelse(shape > 2,
       {b <- pi/shape
        scale^2 * (2*b/sin(2*b) - b^2/sin(b)^2)},
           NaN)    
}
