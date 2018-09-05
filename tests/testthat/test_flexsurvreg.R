context("flexsurvreg model fits")

test_that("Generalized F (p parameter) not identifiable from ovarian data",{
expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf"), "non-finite finite-difference")
expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf.orig"), "non-finite finite-difference")
})

test_that("Generalized gamma fit",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma")
    print(fitg)
    print(fitg, digits=4)
    fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gengamma")
    ovarian2 <- ovarian
    fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian2, dist="gengamma")
    ## GF with "p" fixed at 0
    fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                           fixedpars=4, inits=c(NA,NA,NA,1e-05))
    expect_equal(fitffix$loglik, sum(fitffix$logliki))
    expect_equal(fitffix$res[1:3,"est"], fitg$res[1:3,"est"], tol=1e-03)
    expect_equal(fitffix$res[1:3,2:3], fitg$res[1:3,2:3], tol=1e-02)
})

test_that("Same answers as survreg for Weibull regression",{
    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
    fitws <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
    expect_equal(fitw$loglik, fitws$loglik[1], tol=1e-04)
    expect_equal(fitws$scale, 1 / fitw$res["shape","est"], tol=1e-03)
    expect_equal(as.numeric(coef(fitws)[1]), log(fitw$res["scale","est"]), tol=1e-03)
})

test_that("Exponential",{
    sr <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="exponential")
    fite <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="exp")
    expect_equal(sr$loglik[2], fite$loglik)
    expect_warning(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="exp", sr.control=list(maxiter=2)), "Ran out of iterations")
})

test_that("Log-normal",{
    sr <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="lognormal")
    fitl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="lognormal")
    expect_equal(sr$loglik[2], fitl$loglik)
    expect_warning(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="exp", sr.control=list(maxiter=2)), "Ran out of iterations")
})

test_that("Weighted fits",{
    wt <- rep(1, nrow(ovarian))
    wt[c(1,3,5,7,9)] <- 10
    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", weights=wt)
    fitws <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", weights=wt)
    expect_equal(fitws$loglik[2],fitw$loglik,tol=1e-06)
})

test_that("subset",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, subset=-(1:2), dist="gengamma")
    fitg2 <- flexsurvreg(formula = Surv(ovarian$futime[-(1:2)], ovarian$fustat[-(1:2)]) ~ 1, dist="gengamma")
    expect_equal(fitg$loglik,fitg2$loglik)
})

test_that("na.action",{
    ovarian2 <- ovarian
    ovarian2$futime[1] <- NA
    fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data=ovarian2, na.action=na.omit, dist="gengamma")
    ovarian3 <- ovarian[-1,]
    fitg2 <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data=ovarian3, na.action=na.omit, dist="gengamma")
    expect_equal(fitg$loglik,fitg2$loglik)
    expect_error(flexsurvreg(formula = Surv(futime, fustat) ~ 1, data=ovarian2, na.action=na.fail, dist="gengamma"), "missing values")
})

test_that("Log-normal",{
    fitln <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm")
    expect_equal(fitln$loglik, -97.12174204265681, tol=1e-06)
})

test_that("Gompertz",{
    fitgo <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gompertz", fixedpars=TRUE) # model fit is unstable
    expect_equal(fitgo$loglik, -112.8294446076947, tol=1e-06)
})

test_that("Log-logistic",{
    fitlls <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ age, data = ovarian, dist="loglogistic")
    fitll <- flexsurvreg(formula = Surv(futime, fustat) ~ age, data = ovarian, dist="llogis")
    expect_equal(fitll$loglik, fitlls$loglik[2], tol=1e-06)
})

test_that("Gamma",{
    fitga <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian, dist="gamma")
    expect_equal(fitga$loglik,  -97.86379723453011, tol=1e-06)
})

test_that("Loglikelihoods of flexible distributions reduce to less flexible ones for certain parameters",{
    ## Test distributions reducing to others with fixed pars
    fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                           fixedpars=TRUE, inits=c(0,1,0,1))
    ## GG = GF with p -> 0
    fitffix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf",
                           fixedpars=TRUE, inits=c(0,1,0,1e-08))
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                           fixedpars=TRUE, inits=c(0,1,0))
    expect_equal(fitgfix$loglik, fitffix$loglik, tol=1e-02)
    ## Weib = GG with q=1
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                           fixedpars=TRUE, inits=c(6,0.8,1))
    fitwfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull",
                           fixedpars=TRUE, inits=c(1/0.8,exp(6)))
    expect_equal(fitwfix$loglik, fitgfix$loglik)
    ## Gamma = GG with sig=q
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                           fixedpars=TRUE, inits=c(6,0.5,0.5))
    fitgafix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gamma",
                            fixedpars=TRUE, inits=c(1/0.5^2,exp(-6)/0.5^2))
    expect_equal(fitgafix$loglik,fitgfix$loglik)
    ## Log-normal = GG with q=0
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma",
                           fixedpars=TRUE, inits=c(6,0.8,0))
    fitlfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",
                           fixedpars=TRUE, inits=c(6,0.8))
    expect_equal(fitlfix$loglik,fitgfix$loglik)
    ## Compare with weib/lnorm fit from survreg
    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull",
                        inits=c(1/0.8,exp(6)))
    fitw2 <- survreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
    expect_equal(1 / fitw2$scale, fitw$res["shape","est"], tol=1e-03)
    expect_equal(as.numeric(coef(fitw2)[1]), log(fitw$res["scale","est"]), tol=1e-03)
})

test_that("Fits of flexible distributions reduce to less flexible ones with fixed parameters",{
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma.orig",
                           fixedpars=3, inits=c(NA,NA,1))
    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
    expect_equal(logLik(fitgfix), logLik(fitw), tol=1e-06)
    fitgfix <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gengamma.orig",
                           fixedpars=1, inits=c(1,NA,NA))
    fitga <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="gamma")
    expect_equal(logLik(fitgfix), logLik(fitga), tol=1e-06)

    fite <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="exp")
    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull", inits=c(1, mean(ovarian$futime)), fixedpars=1)

    expect_equal(logLik(fite), logLik(fitw), tol=1e-06)
    expect_equal(fitw$res["scale",1], 1 / fite$res["rate",1], tol=1e-06)
})

fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(rx), data = ovarian, dist="weibull")

test_that("Model fit with covariates",{
    expect_equal(fitg$loglik, -97.3641506645869, tol=1e-06)
    if (interactive()) {
        plot.flexsurvreg(fitg, ci=TRUE)
        plot.flexsurvreg(fitg, X=rbind(c(0), c(1)), ci=TRUE, col="red")
        lines.flexsurvreg(fitg, X=rbind(c(1.1), c(1.2)), ci=TRUE, col="blue")
        plot.flexsurvreg(fitg, type="hazard")
        plot.flexsurvreg(fitg, type="cumhaz")

        fitg1 <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="weibull")
        plot.flexsurvreg(fitg1, type="hazard")

    }
})

test_that("Summary function: alternative ways to supply covariates",{
    expect_equal(summary(fitg, X=c(0), ci=FALSE)[[1]],
                 summary(fitg, newdata=data.frame(rx=1), ci=FALSE)[[1]])
    expect_equal(summary(fitg, X=c(1), ci=FALSE)[[1]],
                 summary(fitg, newdata=data.frame(rx=2), ci=FALSE)[[1]])
    expect_equivalent(summary(fitg, X=matrix(c(0,1),ncol=1), ci=FALSE)[1:2],
                      summary(fitg, newdata=data.frame(rx=c(1,2)), ci=FALSE)[1:2])
    fitg2 <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(rx) + factor(ecog.ps), data = ovarian, dist="weibull")
    expect_equivalent(summary(fitg2, newdata=data.frame(rx=1,ecog.ps=2), ci=FALSE)[[1]][1:2],
                      summary(fitg2, X=c(0,1))[[1]][1:2])
})

test_that("summary with CIs",{
    summ <- summary(fitg, newdata=data.frame(rx=1), B=2, type="survival")
    print(summ)
    expect_true(all(unlist(summ[[1]][,2:4]) <= 1))
    expect_true(all(unlist(summ[[1]][,2:4]) >= 0))
})

test_that("Errors in summary function",{
    expect_error(summary(fitg, newdata=list(foo=1)), "Value of covariate \"rx\" not supplied")
    expect_error(summary(fitg, X=matrix(c(0,1),ncol=2), ci=FALSE), "expected X to be a matrix with 1 column or a vector with 1 element")
    expect_error(summary(fitg, X=matrix(c(0,1),ncol=1), start=1:2, ci=FALSE), "length of \"start\"")
    expect_error(summary(fitg, start=1:2, ci=FALSE), "length of \"start\"")
})


test_that("Model fit with covariates and simulated data",{
    x <- rnorm(500,0,1)
    sim <- rgenf(500, 1.5 - 0.2*x, 1, -0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ x, dist="genf", control=list(maxit=10000))
    if (interactive()) {
        fit # estimate should be -0.2
        summary(fit)
        plot.flexsurvreg(fit)
        lines.flexsurvreg(fit, X=matrix(c(1,2),nrow=2))
        plot(fit)
        plot(fit, type="hazard", min.time=0, max.time=25)
        lines(fit, type="hazard", X=matrix(c(1,2),nrow=2))
        x2 <- factor(rbinom(500, 1, 0.5))
        fit <- flexsurvreg(Surv(simt, dead) ~ x + x2, dist="genf", control=list(maxit=10000))
        plot(fit)
        plot(fit, type="cumhaz")
        plot(fit, type="hazard", min.time=0, max.time=25)
        x3 <- factor(rbinom(500, 1, 0.5))
        fit <- flexsurvreg(Surv(simt, dead) ~ x2 + x3, dist="genf", control=list(maxit=10000))
        fit <- flexsurvreg(Surv(simt, dead) ~ x2, dist="genf", control=list(maxit=10000))
        plot(fit)
        summary(fit, type="hazard", ci=FALSE)
        plot(fit, type="hazard", ci=FALSE)
    }
    x2 <- factor(rbinom(500, 1, 0.5))
    x3 <- rnorm(500,0,1)
    sim <- rgengamma(500, 1.5 + 2*x3, 1, -0.4)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    fit <- flexsurvreg(Surv(simt, dead) ~ x3, dist="gengamma", control=list(maxit=10000))
    fit <- flexsurvreg(Surv(simt, dead) ~ x + x2 + x3, dist="gengamma", control=list(maxit=10000))
    fit <- flexsurvreg(Surv(simt, dead)[1:100] ~ x[1:100] + x2[1:100], dist="gengamma", control=list(maxit=10000), method="BFGS")
    fit <- flexsurvreg(Surv(simt, dead)[1:100] ~ x[1:100], dist="gengamma", control=list(maxit=10000))
    fit
})

test_that("Covariates on ancillary parameters",{
    set.seed(11082012)
    x3 <- rnorm(1500,0,1)
    x4 <- rnorm(1500,0,1)
    x5 <- rnorm(1500,0,1)
    sim <- rgengamma(1500, 1, exp(0.5 + 0.1*x3 + -0.3*x4), -0.4 + 1.2*x5)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)

    fit <- flexsurvreg(Surv(simt, dead) ~ sigma(x3), dist="gengamma", control=list(maxit=10000))
    fit
    cl <- confint(fit)

    ## Cov on ancillary, not on location
    flexsurvreg(Surv(simt, dead) ~ sigma(x3), dist="gengamma")
    flexsurvreg(Surv(simt, dead) ~ 1, anc=list(sigma=~x3), dist="gengamma")

    ## Cov on both location and ancillary
    flexsurvreg(Surv(simt, dead) ~ x3 + sigma(x3), dist="gengamma")
    flexsurvreg(Surv(simt, dead) ~ x3, anc=list(sigma=~x3), dist="gengamma")

    ## More than one covariate on an ancillary parameter
    flexsurvreg(Surv(simt, dead) ~ x3 + sigma(x3) + sigma(x4), dist="gengamma")
    flexsurvreg(Surv(simt, dead) ~ x3, anc=list(sigma=~x3+x4), dist="gengamma")

    ## More than one ancillary parameter with covariates
    flexsurvreg(Surv(simt, dead) ~ x3 + sigma(x3) + sigma(x4) + Q(x5), dist="gengamma")
    x <- flexsurvreg(Surv(simt, dead) ~ x3, anc=list(sigma=~x3+x4, Q=~x5), dist="gengamma")
})

test_that("Various errors",{
    expect_error(flexsurvreg(data = ovarian, dist="genf", inits = c(1,2,3)),"\"formula\" is missing")
    expect_error(flexsurvreg(formula="foo", data = ovarian, dist="genf", inits = c(1,2,3)),"\"formula\" must be a formula")
    expect_error(flexsurvreg(formula= futime ~ fustat, data = ovarian, dist="genf", inits = c(1,2,3)),"Response must be a survival object")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, inits = c(1,2,3)),"Distribution \"dist\" not specified")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist=1, data = ovarian, inits = c(1,2,3)),"\"dist\" should be a")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", data = ovarian, anc=1, inits = c(1,2,3)),"\"anc\" must be a")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,2,3)),"Initial values .+ length")
    expect_warning(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,2,3,4,5), fixedpars=TRUE),"Initial values are a vector length .+ using only the first")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = "foo"),"init.+ must be a numeric vector")
    suppressWarnings({
        expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,2,3,-1)),"Initial value for parameter 4 out of range")
        expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", inits = c(1,-2,3,-1)),"Initial values for parameters 2,4 out of range")
    })
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", fixedpars = c(3,4,5,6,7)), "fixedpars must be TRUE/FALSE")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="genf", fixedpars = "foo"), "fixedpars must be TRUE/FALSE")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=-1), "cl must be a number in")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=1.1), "cl must be a number in")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl=c(1,2)), "cl must be a number in")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="lnorm",cl="foo"), "cl must be a number in")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian, dist="foo"), "\'arg\' should be one of")
    expect_error(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, data = ovarian), "Distribution \"dist\" not specified")
})

test_that("Calling flexsurvreg from within a function",{
    f <- function(){
        ovarian2 <- ovarian
        fitg <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian2, dist="gengamma")
        fitg <- flexsurvreg(formula = Surv(ovarian2$futime, ovarian2$fustat) ~ factor(ovarian2$rx), dist="gengamma",method="Nelder-Mead")
    }
    f()
})

test_that("Calling flexsurvreg from a function within a function",{
    f <- function(){
        ovarian2 <- ovarian
        g <- function(){
            fitw <- flexsurvreg(formula = Surv(futime, fustat) ~ 1, data = ovarian2, dist="weibull")
            fitw <- flexsurvreg(formula = Surv(ovarian2$futime, ovarian2$fustat) ~ factor(ovarian2$rx), dist="weibull")
        }
        g()
    }
    f()
})

## Left-truncation.
## Time passed as arg to initial values is stop - start,
## since, e.g. mean of trunc exponential dist is 1/lam + b, mean par plus trunc point
## time at risk in returned object is currently sum of (stop - start)
## default knot choice for spline - start + quantiles of log dt

test_that("Left-truncation",{
    set.seed(12082012)
    sim <- rgenf(3000, 1.5, 1, -0.4, 0.6)
    dead <- as.numeric(sim<=30)
    simt <- ifelse(sim<=30, sim, 30)
    obs <- simt>3; simt <- simt[obs]; dead <- dead[obs]
    fit <- flexsurvreg(Surv(simt, dead) ~ 1, dist="gengamma")
    summ <- summary(fit, ci=FALSE)
    expect_true(all(summ$time>3))
    if (interactive()) plot(fit, ci=FALSE, xlim=c(0,10))
    fit <- flexsurvreg(Surv(rep(3, length(simt)), simt, dead) ~ 1, dist="gengamma")
    if (interactive()) lines(fit, ci=FALSE, col="blue") # truncated model fits truncated data better.
})

test_that("Interval censoring",{
    set.seed(1)
    simt <- rweibull(1000, 2, 0.5)
    tmin <- simt
    status <- ifelse(simt>0.6, 0, 1)
    simt[status==0] <- 0.6
    tmin <- simt
    tmax <- ifelse(status==1, simt, Inf)
    tmax.sr <- ifelse(status==1, simt, NA) # need -Inf instead of Inf, don't understand why. bug?
## Currently flexsurvreg fails with Inf, just as survreg fails, with Weibull. Consistent with survreg for the moment
    
    sr1 <- survreg(Surv(tmin, tmax.sr, type="interval2") ~ 1, dist="weibull")
    sr2 <- survreg(Surv(tmin, status) ~ 1, dist="weibull")
    expect_equal(sr1$loglik[2], sr2$loglik[2])
    
    fs1 <- flexsurvreg(Surv(tmin, tmax.sr, type="interval2") ~ 1, dist="weibull")
    fs2 <- flexsurvreg(Surv(simt, status) ~ 1, dist="weibull")
    expect_equal(fs1$loglik, fs2$loglik)
    expect_equal(fs1$loglik, sr1$loglik[2])

    ## put an upper bound on censored times
    tmax <- ifelse(status==1, simt, 0.7)
    fs1 <- flexsurvreg(Surv(tmin, tmax, type="interval2") ~ 1, dist="weibull")
    fs2 <- flexsurvreg(Surv(simt, status) ~ 1, dist="weibull")
    expect_true(fs1$loglik != fs2$loglik)

    ### FIXME n events wrong for interval2. 
    
    ## using type="interval"
    status[status==0] <- 3
    fs3 <- flexsurvreg(Surv(tmin, tmax, status, type="interval") ~ 1, dist="weibull")
    expect_equal(fs1$loglik, fs3$loglik)

})

test_that("inits",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", inits=c(6,1,-1))
    fitg2 <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma")
    expect_equal(fitg$loglik, fitg2$loglik, tol=1e-05)
})

test_that("fixedpars",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", fixedpars=1:2, inits=c(6,1,-1))
    expect_equivalent(fitg$res[1:2,"est"], c(6,1))
    expect_equivalent(fitg$res[1:2,"L95%"], c(NA_real_,NA_real_))
})

test_that("aux is ignored if it contains parameters",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma")
    fitg2 <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", aux=list(sigma=1))
    expect_equal(fitg$loglik, fitg2$loglik)
})

test_that("cl",{
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", cl=0.99)
    fitg2 <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gengamma", cl=0.999)
    expect_true(fitg2$res[1,2] < fitg$res[1,2])
    expect_true(fitg2$res[1,3] > fitg$res[1,2])
})

test_that("Relative survival", { 
    bc$bh <- rep(0.01, nrow(bc))

    ## Compare with stata stgenreg, using Weibull PH model

    hweibullPH <- function(x, shape, scale = 1, log=FALSE){
        hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
    }
    HweibullPH <- function(x, shape, scale=1, log=FALSE){
        Hweibull(x, shape=shape, scale=scale^{-1/shape}, log=log)
    }
    custom.weibullPH <- list(name="weibullPH", 
                             pars=c("shape","scale"), location="scale",
                             transforms=c(log, log), inv.transforms=c(exp, exp),
                             inits = function(t){
                                 c(1, median(t[t>0]) / log(2))
                             })

#     fs6b <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.weibullPH, bhazard=bh, dfns=list(h=hweibullPH, H=HweibullPH), fixedpars=TRUE)

    fs6b <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist=custom.weibullPH, bhazard=bh, dfns=list(h=hweibullPH, H=HweibullPH))
    fs6bd <- flexsurvreg(Surv(recyrs, censrec) ~ group, data=bc, dist="weibullPH", bhazard=bh)
    
    expect_equal(log(fs6b$res[1,"est"]), 0.3268327417773233, tol=1e-05)
    expect_equal(log(fs6b$res[2,"est"]), -3.5308925743338038, tol=1e-05)
    expect_equal(fs6b$res["groupMedium","est"], 0.9343799681269026, tol=1e-05)
    expect_equal(fs6b$res["groupPoor","est"], 1.799204192587765, tol=1e-05)

    ## same results as 
    ## cd /home/chris/flexsurv/stata
    ## use stpm/bc
    ## gen rec = 1 - censrec
    ## gen recyrs = rectime / 365
    ## gen bh = 0.01
    ## stset recyrs, failure(censrec)
    ## stgenreg, loghazard([ln_lambda] :+ [ln_gamma] :+ (exp([ln_gamma]) :- 1) :* log(#t)) nodes(100) ln_lambda(group2 group3) bhazard(bh)
### don't know what scale Stata returns the loglik on

})

test_that("warning with strata", { 
    ## need double backslash to escape $
    expect_warning(flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ strata(ovarian$resid.ds), dist="gengamma", inits=c(6,1,-1,0)),
                   "Ignoring \"strata\" function: interpreting \"ovarian\\$resid.ds\" as a covariate on \"mu\"")
})

test_that("Distribution names are case insensitive",{
  fs1 = flexsurvreg(Surv(rectime, censrec)~group,dist="weibull",data=bc)
  fs2 = flexsurvreg(Surv(rectime, censrec)~group,dist="Weibull",data=bc)
 expect_equal(fs1$loglik, fs2$loglik, tol=1e-06)
})

test_that("Weibull hazards from summary are reliable",{
  fs1 = flexsurvreg(Surv(rectime, censrec)~group ,dist="weibull",data=bc)
  output = summary(fs1, t=seq(from=0,to=30000,length.out=100), ci=F, tidy=T)
  expect_true(all(is.finite(output$est)))
})

test_that("RMST/Mean/Median calculations are working",{
  
  fs1 = flexsurvreg(Surv(rectime, censrec)~group ,dist="weibull",data=bc)
  fs2 = flexsurvreg(Surv(rectime, censrec)~group ,dist="exp",data=bc)
  fs3 = flexsurvreg(Surv(rectime, censrec)~group ,dist="llogis",data=bc)
  fs4 = flexsurvreg(Surv(rectime, censrec)~group ,dist="lnorm",data=bc)
  fs5 = flexsurvreg(Surv(rectime, censrec)~group ,dist="gamma",data=bc)
  
  res1 = summary(fs1,t=c(Inf),start=0,type="rmst")
  res2 = summary(fs1,type="mean")
  
  res1_len = length(res1)
  for(i in seq_len(res1_len)){
    expect_equal(
      res1[[i]]$est,
      res2[[i]]$est,
      tolerance=1e-3
    )
  }
  
  # Exponential analytical RMST should be consistent w/ analytical
  # mean.
  expect_equal(summary(fs2,type="mean",tidy=T)$est,summary(fs2,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  
  # Analytical mean should closely match result from integration
  expect_equal(summary(fs1,type="mean",tidy=T)$est,summary(fs1,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs3,type="mean",tidy=T)$est,summary(fs3,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs4,type="mean",tidy=T)$est,summary(fs4,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  expect_equal(summary(fs5,type="mean",tidy=T)$est,summary(fs5,t=Inf,type="rmst",tidy=T)$est, tolerance=1e-3)
  
  
  # RMST of exponential to 100 starting at 0 should be the same
  # as RMST to 200 starting at 100.
  res3 = summary(fs2,t=c(100,200),start=c(0,100),type="rmst")
  res3_len = length(res3)
  for(i in seq_len(res1_len)){
    expect_equal(
      res3[[i]]$est[1],
      res3[[i]]$est[2],
      tolerance=1e-3
    )
  }
  
  expect_warning(
    summary(fs1,t=10,type="mean"),
    "Mean selected, but time specified.  For restricted mean, set type to 'rmst'."
  )
  
  expect_warning(
    summary(fs1,t=10,type="median"),
    "Median selected, but time specified."
  )
  
})

test_that("No events in the data",{
    set.seed(1)
    tmin <- rexp(100, 1) 
    tmax <- tmin + 0.1
    mod <- flexsurvreg(Surv(tmin, tmax, type="interval2") ~ 1, dist="exponential")
    expect_equal(mod$loglik, -337.9815, tol=1e-03)
})

test_that("No censoring in the data",{
    bcev <- bc[bc$censrec==1,]
    mod <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data=bcev, dist="weibull")
    expect_equal(mod$loglik, -477.2455, tol=1e-03)
})

test_that("summary type=quantile is consistent",{
  expect_equal(summary(fitg, type='quantile', quantiles=.5)[[1]][1,2]
               ,summary(fitg, type='median')[[1]][1,1])
  
  expect_equal(summary(fitg, type='quantile', quantiles=.5, start = 50)[[1]][1,2]
               ,summary(fitg, type='median', start = 50)[[1]][1,1])
})

test_that("Errors in summary type=quantile",{
  expect_error(summary(fitg, type='quantile', quantiles=1.5), "Quantiles should not be less than 0 or greater than 1")
  expect_error(summary(fitg, type='quantile', quantiles=-.5), "Quantiles should not be less than 0 or greater than 1")
})

test_that("SEs in summary function",{
    expect_true(is.numeric(summary(fitg, se=TRUE)[[1]]$se))
})

test_that("summary type `link`",{
    expect_equal(summary(fitg, type="link")[["factor(rx)=1"]]$est, 
                 fitg$res["scale","est"])
    expect_equal(summary(fitg, type="link")[["factor(rx)=2"]]$est, 
                 exp(fitg$res.t["scale","est"] + fitg$res.t["factor(rx)2","est"]))
})
