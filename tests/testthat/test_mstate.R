context("Multi-state modelling and prediction")

library(mstate) # masks flexsurv's generic msfit, so call as msfit.flexsurvreg below

### This would all be better tested through a vignette of worked examples
### e.g. to illustrate more flexible models fitting better

bexp <- flexsurvreg(Surv(years, status) ~ trans, data=bosms3, dist="exp") 
tmat <- rbind(c(NA,1,2),c(NA,NA,3),c(NA,NA,NA))
tgrid <- seq(0,14,by=0.1)
bwei <- flexsurvreg(Surv(years, status) ~ trans + shape(trans), data=bosms3, dist="weibull")
bspl <- flexsurvspline(Surv(years, status) ~ trans + gamma1(trans), data=bosms3, k=3)
bexp.markov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans, data=bosms3, dist="exp")

test_that("msfit.flexsurvreg",{
    mexp <- msfit.flexsurvreg(bexp, t=0.01, trans=tmat, tvar="trans")
    summ <- summary.flexsurvreg(bexp, t=0.01, type="cumhaz", ci=FALSE, newdata=list(trans=1:3))
    summ <- as.numeric(unlist(lapply(summ, function(x)x$est[x$time==0.01])))   
    expect_equal(mexp$Haz$Haz[mexp$Haz$time==0.01], summ)
    mwei <- msfit.flexsurvreg(bwei, t=c(0.01, 0.02), trans=tmat, tvar="trans", B=10)
    mspl <- msfit.flexsurvreg(bspl, t=c(0.01, 0.02), trans=tmat, tvar="trans", B=10)
})

## With covariates
set.seed(1)
bosms3$x <- rnorm(nrow(bosms3))
bexp.cov <- flexsurvreg(Surv(years, status) ~ trans + x, data=bosms3, dist="exp")
bexp.markov.cov <- flexsurvreg(Surv(Tstart, Tstop, status) ~ trans + x, data=bosms3, dist="exp")

test_that("newdata in msfit.flexsurvreg",{
    msfit.flexsurvreg(bexp.cov, newdata=list(x=1), t=c(0,5,10), trans=tmat, variance=FALSE)
    msfit.flexsurvreg(bexp.cov, newdata=list(x=2), t=c(0,5,10), trans=tmat, variance=FALSE)
    msfit.flexsurvreg(bexp.cov, newdata=list(x=c(1,2,3)), t=c(0,5,10), trans=tmat, variance=FALSE)
})

test_that("Errors in msfit.flexsurvreg",{
    expect_error(msfit.flexsurvreg(bexp.cov, t=seq(0,150,10), trans=tmat), "Value.* of covariate.* .+ not supplied")
    expect_error(msfit.flexsurvreg(bexp.cov, t=seq(0,150,10), trans=tmat, tvar="foo"), "variable .* not in model")
    expect_error(msfit.flexsurvreg(bexp.cov, newdata=list(x=c(1,2)), t=c(0,5,10), trans=tmat, variance=FALSE), "length of variables .+ must be")   
})

### TODO test values

test_that("pmatrix.fs",{
    pmatrix.fs(bexp.markov, t=c(5,10), trans=tmat)
    pmatrix.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1))
})

test_that("pmatrix.simfs",{
    pmatrix.simfs(bexp, t=5, trans=tmat, M=100)
    pmatrix.simfs(bwei, t=5, trans=tmat, M=100)
    pmatrix.simfs(bexp.cov, t=5, trans=tmat, newdata=list(x=1), M=100)
})

test_that("totlos.fs",{
    totlos.fs(bexp.markov, t=c(5,10), trans=tmat)
    totlos.fs(bexp.markov.cov, t=c(5,10), trans=tmat, newdata=list(x=1))
})

test_that("totlos.simfs",{
    totlos.simfs(bexp, t=5, trans=tmat, M=100)
    totlos.simfs(bwei, t=5, trans=tmat, M=100)
    totlos.simfs(bexp.cov, t=5, trans=tmat, newdata=list(x=1), M=100)
})
