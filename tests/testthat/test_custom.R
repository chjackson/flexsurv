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
        dllogis2 <- eha::dllogis; pllogis2 <- eha::pllogis       
        custom.llogis <- list(name="llogis2",
                              pars=c("shape","scale"),
                              location="scale",
                              transforms=c(log, log),
                              inv.transforms=c(exp, exp),
                              inits=function(t){ c(1, median(t)) })
        fitll <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.llogis, dfns=list(d=dllogis2, p=pllogis2))

        custom.ev <- list(name="EV",
                          pars=c("shape","scale"),
                          location="scale",
                          transforms=c(log, log),
                          inv.transforms=c(exp, exp),
                          inits=function(t){ c(1, median(t)) })
        fitev <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.ev, dfns=list(d=dEV, p=pEV))
###   h(x) = (b/a)(x/a)^(b-1)exp((x / a)^b)
###   H(x) = exp( (x / a)^b) ) - 1


        detach("package:eha")

        fitll.builtin <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="llogis")
        expect_equal(fitll$loglik, fitll.builtin$loglik)
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

### FIXME q is Inf.  is it Surv object bug?  time2 being used
### it's trying to get p from exp(-H(t)) with t = Inf, should get 


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
    fitf <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist=custom.baz, dfns=list(d=dbaz), inits=c(9e-07, 0.12))
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
    custom.foo <- list(pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "\"name\" element of custom distribution list not given")
    custom.foo <- list(name=0, pars=c("rate"), location="rate", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "\"name\" element of custom distribution list should be a string")
    custom.foo <- list(name="foo", location="rate", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "parameter names \"pars\" not given")
    custom.foo <- list(name="foo", pars=2, location="rate", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "parameter names \"pars\" should be a character vector")
    custom.foo <- list(name="foo", pars="rate", transforms=c(log), inv.transforms=c(exp))
    expect_warning(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo), fixedpars=TRUE, inits=1),
                 "location parameter not given, assuming it is the first one")
    custom.foo <- list(name="foo", pars="rate", location="bar", transforms=c(log), inv.transforms=c(exp))
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo), fixedpars=TRUE, inits=1),
                 "location parameter \"bar\" not in list of parameters")
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
                 "\"transforms\" must be a list")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log), inv.transforms=exp, location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "\"inv.transforms\" must be a list")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=list(2), inv.transforms=c(exp), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "some of \"transforms\" are not functions")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log), inv.transforms=list(2), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "some of \"inv.transforms\" are not functions")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log,log), inv.transforms=c(exp,exp), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "transforms vector of length 2, parameter names of length 1")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log), inv.transforms=c(exp,exp), location="rate")
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "inverse transforms vector of length 2, parameter names of length 1")
    custom.foo <- list(name="foo", pars=c("rate"), transforms=c(log), inv.transforms=c(exp), location="rate", inits=1)
    expect_error(flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist=custom.foo, dfns=list(h=hfoo, H=Hfoo)),
                 "\"inits\" element of custom distribution list must be a function")
})

test_that("Errors in form.dp",{
    expect_error(form.dp(dlist=list(), dfns=list()), "Neither density function .+ nor hazard function")
})
