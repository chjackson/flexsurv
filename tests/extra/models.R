
### LIBRARY OF NEW MODELS

## Logistic: as in manual

library(eha)
custom.llogis <- list(name="llogis",  pars=c("shape","scale"), location="scale",
                      transforms=c(log, log), inv.transforms=c(exp, exp),
                      inits=function(t){ c(1, median(t)) })
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.llogis)
plot(fs1)
fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
lines(fs2, col="blue")
## Log-logistic / proportional odds model fits better


## Gompertz-Makeham: wrapping an existing distribution
## Need proper vectorisation to work for covariates

library(eha)

dmakeham3 <- function(x, shape1, shape2, scale, ...)  {
    dmakeham(x, shape=c(shape1, shape2), scale=scale, ...)
}
dmakeham3 <- Vectorize(dmakeham3) 

pmakeham3 <- function(q, shape1, shape2, scale, ...)  {
    pmakeham(q, shape=c(shape1, shape2), scale=scale, ...)
}
pmakeham3 <- Vectorize(pmakeham3)

custom.makeham <- list(name="makeham3",
                       pars=c("shape1","shape2","scale"),
                       location="scale",
                       transforms=c(log, log, log),
                       inv.transforms=c(exp, exp, exp),
                       inits=function(t){ c(0.01, 0.02, median(t)) })

x <- rnorm(1000); beta <- 0.1
msim <- rmakeham(1000, shape=c(0.01, 0.02), scale=10*exp(beta*x))
fit <- flexsurvreg(Surv(msim, rep(1, 1000)) ~ x, dist=custom.makeham)
fit # This one needs tight initial values to converge

flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.makeham)
## shape2 and scale parameters appear to be not identifiable
flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="gompertz")

## Weibull proportional hazards parameterisation

hweibullPH <- function(x, shape, scale = 1, log=FALSE){
    hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

HweibullPH <- function(x, shape, scale=1, log=FALSE){
    Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
}

custom.weibullPH <- list(name="weibullPH", pars=c("shape","scale"), location="scale",
                         transforms=c(log, log), inv.transforms=c(exp, exp),
                         inits = function(t){ ## TODO inits
                             lt <- log(t[t>0])
                             c(1, exp(mean(lt) + 0.572))
                         })

bcw <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibull")
bcwph <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.weibullPH)
plot(bcw, type="hazard")
lines(bcwph, type="hazard", col="blue")   # indistinguishable, as expected. same model, different param
## since HR in PH model = exp(AFT effect)^(-shape)




### GG prop haz model

hgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

HgengammaPH <- function(x, dummy, mu=0, sigma=1, Q){
    dummy * Hgengamma(x=x, mu=mu, sigma=sigma, Q=Q)
}

custom.gengammaPH <- list(name="gengammaPH", pars=c("dummy","mu","sigma","Q"), location="dummy",
                          transforms=c(log, identity, log, identity),
                          inv.transforms=c(exp, identity, exp, identity),
                          inits=function(t){
                              lt <- log(t[t>0])
                              c(1, mean(lt), sd(lt), 0)
                          })

bc.gg <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="gengamma")
bc.ggph <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.gengammaPH, fixedpars=1)

plot(bc.gg, type="hazard")
lines(bc.ggph, type="hazard", col="blue")

summ <- summary(bc.ggph, type="hazard", ci=FALSE)
summ[[2]]$est / summ[[1]]$est  #  should equal
exp(bc.ggph$res["groupMedium","est"])
summ[[3]]$est / summ[[1]]$est  #  
exp(bc.ggph$res["groupPoor","est"])




## Log hazard as spline function of log time 
## User supplies knots and inits

hsurvspline.lh <- function(x, gamma, knots=c(-10,10)){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    loghaz <- rowSums(basis(knots, log(x)) * gamma)
    exp(loghaz)
}

hsurvspline.lh3 <- unroll.function(hsurvspline.lh, gamma=0:2)

custom.hsurvspline.lh3 <- list(
    name = "survspline.lh3",
    pars = c(paste0("gamma",0:2)),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms=rep(c(identity), 3)
    )

dtime <- log(bc$recyrs)[bc$censrec==1]
kmin <- min(dtime)
kmax <- max(dtime)
ak <- list(knots=c(kmin,0,kmax))

### Exponential, constant log hazard
flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, inits=0.14, dist="exp") # rate 0.14
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=ak,
                   inits=c(log(0.14), 0, 0), dist=custom.hsurvspline.lh3, fixedpars=2:3)
### Weibull, log hazard is linear function of log time.  right answer but very slow
flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, dist="weibull")
fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=ak,
            inits=c(log(0.14), 0, 0), dist=custom.hsurvspline.lh3, fixedpars=3)

system.time(
fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=ak,
                   inits=c(log(0.14), 0, 0), dist=custom.hsurvspline.lh3, method="L-BFGS-B",
                   lower=c(-Inf,-Inf,-0.5), upper=c(Inf,Inf,0.5), control=list(trace=1,REPORT=1))
) # 146 sec 

system.time(
fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, aux=ak,
                   inits=c(0, 1, 0.3, 0, 0), dist=custom.hsurvspline.lh3, method="L-BFGS-B",
                   lower=c(-Inf,-Inf,-0.5), upper=c(Inf,Inf,0.5), control=list(trace=1,REPORT=1))
)

### without constraint, gammas get too low (-7,-7,-1) and integral diverges 
### log hazards shoot off to minus infinity. get massive values for gamma, both pos and neg
### Runs if constrain parameters: mle of gamma2:  0.2918  (0.1877, 0.3958)
### though takes 5 minutes to run 
fs3 <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, aux=ak,
                   inits=c(0, 1, 0.3, 0, 0), dist=custom.hsurvspline.lh3, method="L-BFGS-B",
                   lower=c(-Inf,-Inf,-0.5), upper=c(Inf,Inf,0.5), control=list(trace=1,REPORT=1))

plot(sp1, type="hazard")
lines(fs3, type="hazard", col="blue")


### any chance of more knots?

hsurvspline.lh5 <- unroll.function(hsurvspline.lh, gamma=0:4)

custom.hsurvspline.lh5 <- list(
    name = "survspline.lh5",
    pars = c(paste0("gamma",0:4)),
    location = c("gamma0"),
    transforms = rep(c(identity), 5), inv.transforms=rep(c(identity), 5)
    )

ak <- list(knots=c(kmin,quantile(dtime, c(0.25, 0.5, 0.75)), kmax))

system.time(
fs5 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=ak,
                   inits=c(log(0.14), 0, 0, 0, 0), dist=custom.hsurvspline.lh5,
                   method="L-BFGS-B", lower=c(-Inf,-Inf,-0.5,-0.5,-0.5), upper=c(Inf,Inf,0.5,0.5, 0.5),
                   control=list(trace=1,REPORT=1))
) # no SEs at convergence.   model probably not much better anyway

plot(fs5, type="hazard")
sp6 <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=6, scale="odds")
## Can we show that models with loads of knots can be made to fit arbitrarily well
plot(sp6, type="cumhaz") # cum haz plot shows well.  bandwidth choice hard in haz plots 
plot(sp6, type="hazard", ylim=c(0,0.5))


### What about log hazard as spline function of (not log) time

hsurvspline.lh <- function(x, gamma, knots=c(-10,10)){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    loghaz <- rowSums(basis(knots, x) * gamma)
    exp(loghaz)
}

hsurvspline.lh3 <- unroll.function(hsurvspline.lh, gamma=0:2)

custom.hsurvspline.lh3 <- list(
    name = "survspline.lh3",
    pars = c(paste0("gamma",0:2)),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms=rep(c(identity), 3)
    )

dtime <- bc$recyrs[bc$censrec==1]
kmin <- min(dtime); kmax <- max(dtime)
ak <- list(knots=c(kmin,median(dtime),kmax))

fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=ak,
                   inits=c(0, 0, 0), dist=custom.hsurvspline.lh3, fixedpars=3,
                   method="L-BFGS-B", lower=c(-Inf,-Inf,-0.5), upper=c(Inf,Inf,0.5),
                   )




## Log hazard as fractional polynomial function of (not log) time

hfp.lh <- function(x, gamma, powers){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    basis <- cbind(1, bfp(x, powers))
    loghaz <- rowSums(basis * gamma)
    exp(loghaz)
}

hfp.lh3 <- unroll.function(hfp.lh, gamma=0:2)

custom.hfp.lh3 <- list(
    name = "fp.lh3",
    pars = c(paste0("gamma",0:2)),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms=rep(c(identity), 3)
    )

aux <- list(powers=c(0,1)) # a log and a linear term
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=aux,
                   inits=c(-2, 0, 0), dist=custom.hfp.lh3)
plot(fs1, type="hazard", ylim=c(0,0.3))  # plot.survival method too slow, but this one seems to fit.
## todo is an extra term worth it?


## Log cumulative hazard as fractional polynomial function of (not log) time
## This can be differentiated, but hazard can't be integrated

Hfp.lh <- function(x, gamma, powers){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    basis <- cbind(1, bfp(x, powers))
    logch <- rowSums(basis * gamma)
    exp(logch)
}

hfp.lh <- function(x, gamma, powers){
    if(!is.matrix(gamma)) gamma <- matrix(gamma, nrow=1)
    lg <- nrow(gamma)
    nret <- max(length(x), lg)
    gamma <- apply(gamma, 2, function(x)rep(x,length=nret))
    x <- rep(x, length=nret)
    basis <- cbind(1, bfp(x, powers))
    dbasis <- cbind(0, dbfp(x, powers))
    loghaz <- rowSums(basis * gamma)
    dloghaz <- rowSums(dbasis * gamma)
    ret <- ifelse(dloghaz <= 0, 0, dloghaz * exp(loghaz))
    ret
}

hfp.lh3 <- unroll.function(hfp.lh, gamma=0:2)
Hfp.lh3 <- unroll.function(Hfp.lh, gamma=0:2)

custom.hfp.lh3 <- list(
    name = "fp.lh3",
    pars = c(paste0("gamma",0:2)),
    location = c("gamma0"),
    transforms = rep(c(identity), 3), inv.transforms=rep(c(identity), 3)
    )

aux <- list(powers=c(1,0))
## null model: constant hazard:  dloghaz=0
fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bc, aux=aux,
                   inits=c(-2, 0.1, 0.01), dist=custom.hfp.lh3)
## non-finite finite difference error
