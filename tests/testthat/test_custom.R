context("Custom distributions in flexsurvreg")

## These previously didn't work when the d or h functions are
## defined in functions or testthat environments.  To find the h or d
## function in form.dp, get() and exists() would need to look in
## parent.frame(2), i.e. the grandparent calling environment, but
## inherits=TRUE searches back through _enclosing_ (definition)
## environments.
## Even with d/h functions defined at top level, these still didn't
## work when called from R CMD check.

## Solve by adding dfns argument to flexsurvreg, passing functions through.

test_that("Custom distributions from another package",{
    if (is.element("eha", installed.packages()[,1])) {
        library(eha)
        custom.llogis <- list(name="llogis",
                              pars=c("shape","scale"),
                              location="scale",
                              transforms=c(log, log),
                              inv.transforms=c(exp, exp),
                              inits=function(t){ c(1, median(t)) })
        fitll <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.llogis)
    }
})
    
test_that("Custom distributions: hazard and cumulative hazard",{
    hfoo <- function(x, rate=1){ rate }
    Hfoo <- function(x, rate=1){ rate*x }
    custom.foo <- list(name="foo", pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp), inits=function(t)1/median(t))
    fitf <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo))
    fite <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
    expect_equal(fitf$loglik, fite$loglik)
})

test_that("Custom distributions: hazard only",{
    hbar <- function(x, rate=1){ rate }
    hbar <- Vectorize(hbar)
    custom.bar <- list(name="bar", pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp), inits=function(t)1/mean(t))
    fitf <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.bar, dfns=list(h=hbar))
    fite <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
    expect_equal(fitf$loglik, fite$loglik)

    ## with covariates.  Approximation is less good.  Many more integrations needed
    fitf <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist=custom.bar, dfns=list(h=hbar))
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    expect_equal(fitf$loglik, fite$loglik, tol=1e-05)

    ## options to integrate()
    flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.bar, dfns=list(h=hbar), integ.opts=list(rel.tol=0.01))
    flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.bar, dfns=list(h=hbar), integ.opts=list(subdivisions=200))

    ## OK to omit inits from custom list if supply it to flexsurvreg
    custom.bar <- list(name="bar", pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp))
    flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.bar, dfns=list(h=hbar), inits=0.001)
})

test_that("Custom distributions: density only",{
    custom.baz <- list(name="baz", pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp), inits=function(t)1/mean(t))
    dbaz <- function(x, rate=1, log=FALSE){
        if (log) {log(rate) - x*rate} else rate*exp(-rate*x)
    }
    dbaz <- Vectorize(dbaz)
    fitf <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.baz, dfns=list(d=dbaz))
    fite <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
    expect_equal(fitf$loglik, fite$loglik)

    ## with covariates
    fitf <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist=custom.baz, dfns=list(d=dbaz))
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    expect_equal(fitf$loglik, fite$loglik, tol=1e-05)
})

## Integration breaks for the Gompertz

if (0) {
    hbar2 <- function(x, a=1, b=1){ b*exp(a*x) }
    hbar2 <- Vectorize(hbar2)
    custom.bar2 <- list(name="bar2", pars=c("a","b"), location="b", transforms=c(identity,log), inv.transforms=c(identity,exp), inits=function(t)c(0.001, 1/mean(t)))
    fitf <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.bar2, control=list(trace=1,REPORT=1,reltol=1e-16,maxit=10000))
    fitgo <- flexsurvreg(formula = Surv(futime, fustat) ~ 1,
                         data = ovarian, dist="gompertz",
                         control=list(trace=1,REPORT=1,reltol=1e-16))
}

## test_that("Custom distributions: hazard only, unknown function, more than one arg",{
## })
## TODO a more complicated one

test_that("Errors in custom distributions",{
    hfoo <- function(x, rate=1){ rate }
    Hfoo <- function(x, rate=1){ rate*x }
    custom.foo <- list(name="foo", pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "\"inits\" not supplied, and no function to estimate them found")
    custom.foo <- list(name="foo", pars=c("rate"), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "transforms not given")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "transforms not given")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=log, inv.transforms=exp, location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "transforms not given")
})
