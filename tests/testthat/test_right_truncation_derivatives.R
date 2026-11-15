test_that("right-truncated fits use the conditional likelihood and curvature", {
    old <- options(flexsurv.test.analytic.derivatives=TRUE)
    on.exit(options(old))
    set.seed(20261113)
    t <- rexp(1000, .2)
    cases <- list(list(t=c(.2,.4,.5,.7), cutoff=1), list(t=t[t<=5], cutoff=5))
    for (case in cases) {
        t <- case$t
        cutoff <- case$cutoff
        fit <- flexsurvreg(Surv(t)~1, dist="exp", rtrunc=rep(cutoff,length(t)),
                          weights=rep(1,length(t)))
        nll <- function(p) -sum(dexp(t,exp(p),log=TRUE) -
                               pexp(cutoff,exp(p),log.p=TRUE))
        ref <- optimize(nll, log(c(.001,10)), tol=1e-10)
        x <- exp(ref$minimum)*cutoff
        logse <- 1/sqrt(length(t)*(1-x^2*exp(x)/expm1(x)^2))
        expect_equal(unname(fit$res["rate","est"]), exp(ref$minimum), tolerance=1e-4)
        expect_equal(unname(sqrt(fit$cov[1,1])), logse, tolerance=1e-3)
        expect_equal(fit$opt$convergence, 0L)
        expect_null(fit$cov.robust)
        expect_null(fit$deriv.test)
        expect_null(fit$hess.test)
    }
})
