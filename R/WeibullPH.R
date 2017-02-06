### Hazard and cumulative hazard functions for R built in
### distributions.  Where possible, use more numerically stable
### formulae than d/(1-p) and -log(1-p)


### Weibull with proportional hazards
## haz(alpha, lambda) for weibull PH is alpha * lam * x^{alpha-1}
## haz(shape, scale) for weibull PH is shape * scale * x^{shape-1}

dweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    dweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

pweibullPH <- function(q, shape, scale = 1,
                       lower.tail=TRUE, log.p=FALSE) {
    pweibull(q, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}

qweibullPH <- function(p, shape, scale = 1, lower.tail=TRUE, log.p=FALSE) {
    qweibull(p, shape=shape, scale=scale^{-1/shape},
             lower.tail=lower.tail, log.p=log.p)
}

hweibullPH <- function(x, shape, scale = 1, log=FALSE) {
    hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

HweibullPH <- function(x, shape, scale=1, log=FALSE) {
    Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

rweibullPH <- function(n, shape, scale=1) {
    rweibull(n, shape=shape, scale=scale^{-1/shape})
}
