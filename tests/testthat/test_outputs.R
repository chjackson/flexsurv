test_that("normboot.flexsurvreg",{
    fite <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="exp")
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=0))
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=50))
    set.seed(1); b2 <- normboot.flexsurvreg(fite, B=10, X=matrix(50,nrow=1))
    expect_equivalent(b1, b2)
    set.seed(1); b1 <- normboot.flexsurvreg(fite, B=10, newdata=list(age=c(0,50)))
    set.seed(1); b2 <- normboot.flexsurvreg(fite, B=10, X=matrix(c(0,50),nrow=2))
    expect_equivalent(b1, b2)   

    ## return cov effs, not adjusted
    set.seed(1)
    normboot.flexsurvreg(fite, B=5, raw=TRUE)
    set.seed(1)
    normboot.flexsurvreg(fite, B=5, raw=TRUE, transform=TRUE)

    ## no covs
    fite <- flexsurvreg(Surv(futime, fustat) ~ 1, data = ovarian, dist="exp")
    normboot.flexsurvreg(fite, B=5)
    normboot.flexsurvreg(fite, B=5, transform=TRUE)   
})


test_that("custom function in summary.flexsurvreg",{
    fitw <- flexsurvreg(Surv(futime, fustat) ~ age, data = ovarian, dist="weibull")     

    median.weibull <- function(t, start, shape, scale) { qweibull(0.5, shape=shape, scale=scale) }
    summ <- summary(fitw, newdata=list(age=50), fn=median.weibull, t=1, B=10)
    expect_equal(summ[[1]][1,"est"], 1575.803185910278, tol=1e-04)

    summ <- summary(fitw, newdata=list(age=50), fn=median.weibull, t=c(1,2,3), B=10)
    expect_equal(summ[[1]][1,"est"], summ[[1]][2,"est"])

    mean.weibull <- function(shape, scale=1) { scale * gamma(1 + 1/shape) }   
    median.weibull <- function(t, start, shape, scale) { scale * log(2)^(1/shape) }

})
