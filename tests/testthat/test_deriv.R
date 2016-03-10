context("Accuracy of analytic derivatives")

if (require("numDeriv")) { 

deriv_error <- function(object){
    if (!isTRUE(getOption("flexsurv.test.analytic.derivatives")))
        stop("flexsurv.test.analytic.derivatives option not set")
    object$deriv.test$error
}

options(flexsurv.test.analytic.derivatives=TRUE)
err <- 1e-04

test_that("Analytic derivatives match numeric",{
    fite <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="exp", inits=0.002)
    expect_lt(deriv_error(fite), err)
    fite <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(ovarian$rx) + factor(ovarian$resid.ds), dist="exp", inits=c(0.002,0.01,-0.001))
    expect_lt(deriv_error(fite), err)

    ov2 <- ovarian[ovarian$futime>200,]
    fite <- flexsurvreg(formula = Surv(rep(200,nrow(ov2)), ov2$futime, ov2$fustat) ~ factor(ov2$rx) + factor(ov2$resid.ds), dist="exp", inits=c(0.002,0.01,-0.001))
    expect_lt(deriv_error(fite), err)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="weibull")
    expect_lt(deriv_error(fitw), err)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(ovarian$rx) + factor(ovarian$resid.ds), dist="weibull")
    expect_lt(deriv_error(fitw), err)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ shape(factor(ovarian$rx)), dist="weibull")    
    expect_lt(deriv_error(fitw), err)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(ovarian$rx), dist="weibull")    

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(ovarian$resid.ds) + shape(factor(ovarian$rx)), inits=c(0.1,0.1,0.2,0.3), dist="weibull")    
    expect_lt(deriv_error(fitw), err)

    ov2 <- ovarian[ovarian$futime>200,]
    fitw <- flexsurvreg(formula = Surv(rep(200,nrow(ov2)), ov2$futime, ov2$fustat) ~ factor(ov2$rx), dist="weibull")
    expect_lt(deriv_error(fitw), err)

    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gompertz", inits=c(0.001, 0.0005), fixedpars=TRUE)
    expect_lt(deriv_error(fitg), err)
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gompertz")
    expect_lt(deriv_error(fitg), err)
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="gompertz", method="Nelder-Mead")
    expect_lt(deriv_error(fitg), err)
    fitg # With derivs, finds wrong MLE if initialize shape at zero, so default to 0.001 instead.
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ factor(ovarian$rx) + factor(ovarian$resid.ds), dist="gompertz")
    expect_lt(deriv_error(fitg), err)
    ov2 <- ovarian[ovarian$futime>200,]
    fitg <- flexsurvreg(formula = Surv(rep(200,nrow(ov2)), ov2$futime, ov2$fustat) ~ factor(ov2$rx) + factor(ov2$resid.ds), dist="gompertz") # truncation
    expect_lt(deriv_error(fitg), err)
    fitg <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ shape(factor(ovarian$resid.ds)), dist="gompertz")
    expect_lt(deriv_error(fitg), err)

    fitl <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ 1, dist="llogis")
    expect_lt(deriv_error(fitl), err)

    
## fixedpars
    fitg <- flexsurvreg(formula = Surv(rep(200,nrow(ov2)), ov2$futime, ov2$fustat) ~ factor(ov2$rx) + factor(ov2$resid.ds), dist="gompertz", fixedpars=2)
    expect_lt(deriv_error(fitg), err)

    fitl <- flexsurvreg(formula = Surv(rep(200,nrow(ov2)), ov2$futime, ov2$fustat) ~ factor(ov2$rx) + factor(ov2$resid.ds), dist="llogis", fixedpars=2)
    expect_lt(deriv_error(fitl), err)
    
    ## spline, hazard. covs, no trunc.
    bc$foo <- factor(sample(1:3, nrow(bc), replace=TRUE))
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=0, scale="odds")
    expect_lt(deriv_error(spl), err)
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group, data=bc, k=1, scale="odds")
    expect_lt(deriv_error(spl), err)
    spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + gamma1(foo), data=bc, k=1, scale="odds")
    expect_lt(deriv_error(spl), err)

    ## truncation
    bc <- bc[bc$recyrs>2,]
    (spl <- flexsurvspline(Surv(recyrs, censrec) ~ 1, data=bc, k=0))
    expect_lt(deriv_error(spl), err)
    (spl <- flexsurvspline(Surv(rep(0, nrow(bc)), recyrs, censrec) ~ 1, data=bc, k=0))
    expect_lt(deriv_error(spl), err)
    
    ## performance increase with derivatives
    time.deriv <- system.time(spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=1, scale="odds"))["elapsed"]
    time.noderiv <- system.time(spl <- flexsurvspline(Surv(recyrs, censrec) ~ group + foo, data=bc, k=1, scale="odds", method="Nelder-Mead"))["elapsed"]
    expect_true(time.deriv < time.noderiv)
})

options(flexsurv.test.analytic.derivatives=FALSE)

}
