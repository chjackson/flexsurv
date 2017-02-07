context("Distribution functions and utilities")

## note - standard q fns in R return zero for p=0 for positive dists, but -Inf for real dists.

tol <- 1e-06


test_that("Exponential hazards",{
    expect_equal(hexp(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,2,3)), c(0,NaN,3,NA,0,3,1))
    expect_equal(hexp(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,2,3), log=TRUE), log(c(0,NaN,3,NA,0,3,1)))
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3)), c(0,NaN,Inf, NA,0,0, 1, 4, 12))
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3), log=TRUE), log(c(0,NaN,Inf, NA,0,0, 1, 4, 12)))
    expect_equal(dexp(c(1,2,3), c(2,3,4)) / (1 - pexp(c(1,2,3), c(2,3,4))), hexp(c(1,2,3), c(2,3,4)))
    expect_equal(-log(1 - pexp(c(1,2,3), c(2,3,4))), Hexp(c(1,2,3), c(2,3,4)))
    expect_equal(log(-log(1 - pexp(c(1,2,3), c(2,3,4)))), Hexp(c(1,2,3), c(2,3,4), log=TRUE))
    expect_warning(hexp(c(1,1,1), c(-1, 0, 2)), "Negative rate")
    expect_warning(dexp(c(1,1,1), c(-1, 0, 2)), "NaNs produced")
})

test_that("Weibull hazards",{
    expect_equal(hweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)), 
         dweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)) / (1 - pweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3))))
    expect_equal(Hweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3)),
         -pweibull(c(1,1,1,1), c(1,1,2,2), c(1,3,1,3), lower.tail=FALSE, log.p=TRUE))
    ## exponential reduction
    expect_equal(hweibull(c(-Inf, NaN, Inf, NA, -1, 0, 1), c(1,1,1), c(1,2,3)), c(0,NaN,1/3,NA,0,1/3,1))
### positive shape, increasing 
    hweibull(c(-Inf, -1, 0, 1, 2, Inf), 2, 2.5)
### neg shape, decreasing 
    hweibull(c(-Inf, -1, 0, 1, 2, Inf), 0.5, 1)
    expect_equal(Hexp(c(-Inf, NaN, Inf,  NA, -1, 0,  1, 2, 4), c(1,2,3)), c(0,NaN,Inf, NA,0,0, 1, 4, 12))
})



test_that("WeibullPH",{
    a <- 0.1; m <- 2
    b <- m^(-1/a)
    x <- c(-Inf, NaN, NA, -1, 0, 1, 2, Inf)
    expect_equal(dweibullPH(x, a, m), dweibull(x, a, b))
    expect_equal(dweibullPH(x, a, m), dweibull(x, a, b), log=TRUE)
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b))
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b), log.p=TRUE)
    expect_equal(pweibullPH(x, a, m), pweibull(x, a, b), lower.tail=FALSE)
    qq <- c(0, 0.5, 0.7, 1)
    expect_equal(qweibullPH(qq, a, m), qweibull(qq, a, b))
    expect_equal(hweibullPH(x, a, m), hweibull(x, a, b))
    expect_equal(hweibullPH(x, a, m), hweibull(x, a, b), log=TRUE)
    expect_equal(HweibullPH(x, a, m), Hweibull(x, a, b))
    expect_equal(HweibullPH(x, a, m), Hweibull(x, a, b), log=TRUE)
    set.seed(1); x1 <- rweibull(10, a, b)
    set.seed(1); x2 <- rweibullPH(10, a, m)
    expect_equal(x1, x2)

    fitw <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ rx, data = ovarian, dist="weibull")
    fitwp <- flexsurvreg(formula = Surv(ovarian$futime, ovarian$fustat) ~ rx, data = ovarian, dist="weibullPH")
    expect_equal(fitw$res["shape","est"], fitwp$res["shape","est"], tol=1e-06)
    expect_equal(fitw$res["scale","est"], fitwp$res["scale","est"]^(-1/fitwp$res["shape","est"]), tol=1e-05)
    expect_equal(coef(fitw)["rx"], -coef(fitwp)["rx"] / fitwp$res["shape","est"])
    
})
